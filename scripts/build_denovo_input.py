#!/usr/bin/env python
"""Build the de novo VHH design input folder for HER2.

Writes ``data/denovo-vhh/`` with the three files a design spec needs:

* ``target_her2.cif``    — HER2 domain IV, chains C (489-580) and D (591-607),
  taken from the tutorial's 1N8Z crop. Trastuzumab's own chains A/B are dropped,
  so the design sees only the antigen.
* ``framework_vhh.cif``  — chain H of the caplacizumab VHH (86 residues, CDRs
  already excised, gaps at 23/38/76). Its PMVK antigen is dropped.
* ``design.json``        — copied from ``specs/denovo-vhh/`` (the tracked copy).

The epitope is pinned by the ``C_HOT`` atoms in the spec: the HER2 atoms nearest
trastuzumab's CDRH3 in 1N8Z, which is what makes the VHH compete with
trastuzumab in the competitive-binding assay.

    python scripts/build_denovo_input.py
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import gemmi

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "data" / "docs" / "Data"
SPEC = ROOT / "specs" / "denovo-vhh" / "design.json"
OUT = ROOT / "data" / "denovo-vhh"

HER2_SRC = DOCS / "12-her2" / "her2_fv_1n8z.cif"
VHH_SRC = DOCS / "03-vhh" / "pmvk_caplacizumab.cif.zst"

HER2_CHAINS = {"C", "D"}  # HER2 domain IV; A/B are trastuzumab VL/VH
VHH_CHAINS = {"H"}  # the VHH framework; A is PMVK, its original antigen


def _read(path: Path) -> gemmi.Structure:
    """Read an mmCIF, transparently decompressing .zst with the zstd CLI."""
    if path.suffix == ".zst":
        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tmp:
            subprocess.run(["zstd", "-dqc", str(path)], stdout=tmp, check=True)
            path = Path(tmp.name)
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    return st


def _is_standard_aa(name: str) -> bool:
    info = gemmi.find_tabulated_residue(name)
    return bool(info and info.is_amino_acid())


def keep_chains(
    src: Path, keep: set[str], dest: Path, *, drop_placeholders: bool = False
) -> dict[str, int]:
    """Write ``src`` to ``dest`` with only the named chains.

    ``drop_placeholders`` removes residues that are not standard amino acids.
    The tutorial's VHH file encodes its CDRs as residues *named* for their
    allowed length (``13-13``, ``9-10``, ``12-21`` at positions 23, 38 and 76),
    which is how it works as a direct ``--input``. Our spec states the lengths
    in its own contig instead, so those placeholders have to go — otherwise the
    module tries to look them up in the CCD and rejects the input with
    "Failed to compute leaving groups". Removing them leaves real numbering
    gaps, exactly like the shipped adalimumab framework.
    """
    src_st = _read(src)
    st = gemmi.Structure()
    st.name = src_st.name
    st.cell, st.spacegroup_hm = src_st.cell, src_st.spacegroup_hm
    model = gemmi.Model('1')
    for chain in src_st[0]:
        if chain.name not in keep:
            continue
        new_chain = gemmi.Chain(chain.name)
        dropped = []
        for res in chain:
            if drop_placeholders and not _is_standard_aa(res.name):
                dropped.append(f"{res.seqid.num}:{res.name}")
                continue
            new_chain.add_residue(res)
        if dropped:
            print(f"  chain {chain.name}: dropped CDR placeholders {dropped}")
        model.add_chain(new_chain)
    st.add_model(model)
    if not st[0]:
        raise SystemExit(f"{src}: none of the chains {sorted(keep)} were found")
    st.setup_entities()
    dest.parent.mkdir(parents=True, exist_ok=True)
    st.make_mmcif_document().write_file(str(dest))
    return {ch.name: len(ch) for ch in st[0]}


def main() -> int:
    if not HER2_SRC.exists() or not VHH_SRC.exists():
        raise SystemExit(
            "tutorial data missing — run: bash scripts/fetch_data.sh small"
        )

    target = keep_chains(HER2_SRC, HER2_CHAINS, OUT / "target_her2.cif")
    print(f"target_her2.cif    chains {target}")

    framework = keep_chains(
        VHH_SRC, VHH_CHAINS, OUT / "framework_vhh.cif", drop_placeholders=True
    )
    print(f"framework_vhh.cif  chains {framework}")

    shutil.copy(SPEC, OUT / "design.json")
    print(f"design.json        copied from {SPEC.relative_to(ROOT)}")

    print(
        f"\nready: {OUT.relative_to(ROOT)}\n"
        "next:  WEIGHT=antibody N=100 SEED=1 "
        "INPUT=data/denovo-vhh/design.json OUT=runs/denovo-vhh-s1 bash scripts/design.sh"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
