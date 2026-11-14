"""Check that a design binds where trastuzumab binds.

Epitope preservation is one of the three judged axes, measured by competitive
binding against trastuzumab. A de novo binder that grips HER2 somewhere else
scores zero on it, however good its affinity — so designs are checked against
the trastuzumab contact patch before they are ranked.

The reference patch is taken from the design spec's ``C_HOT`` atoms (the HER2
atoms nearest trastuzumab's CDRH3 in 1N8Z) and widened to every HER2 residue
within :data:`EPITOPE_RADIUS` of trastuzumab in that structure.

    python -m abhack.epitope runs/denovo-vhh-s1
"""

from __future__ import annotations

import argparse
from pathlib import Path

import gemmi
import polars as pl

from abhack.constants import DOCS

#: Hotspot residues named in specs/denovo-vhh/design.json, as (chain, seqid).
HOTSPOTS: tuple[tuple[str, int], ...] = (("C", 572), ("C", 573), ("D", 591), ("D", 593))

#: A design residue this close to a HER2 residue counts as contacting it (A).
CONTACT_CUTOFF = 4.5
#: How far from the hotspots a HER2 residue may be and still count as epitope (A).
EPITOPE_RADIUS = 10.0

REFERENCE_CIF = DOCS / "Data" / "12-her2" / "her2_fv_1n8z.cif"
ANTIGEN_CHAINS = ("C", "D")


def _atoms(chain) -> list:
    return [(res, atom) for res in chain for atom in res]


def epitope_residues(
    cif: Path | str = REFERENCE_CIF, radius: float = EPITOPE_RADIUS
) -> set[tuple[str, int]]:
    """HER2 residues within ``radius`` of any hotspot atom in the reference."""
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    model = st[0]

    anchors = []
    for chain_name, seqid in HOTSPOTS:
        if chain_name not in [c.name for c in model]:
            continue
        for res in model[chain_name]:
            if res.seqid.num == seqid:
                anchors.extend(atom.pos for atom in res)

    if not anchors:
        raise ValueError(f"no hotspot residues found in {cif}")

    out: set[tuple[str, int]] = set()
    for chain_name in ANTIGEN_CHAINS:
        if chain_name not in [c.name for c in model]:
            continue
        for res in model[chain_name]:
            if any(atom.pos.dist(a) <= radius for atom in res for a in anchors):
                out.add((chain_name, res.seqid.num))
    return out


def _first_residue_numbers(cif: Path | str) -> dict[str, int]:
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    return {ch.name: ch[0].seqid.num for ch in st[0] if len(ch)}


def contacts(cif_path: Path | str, binder_chain: str = "H") -> dict:
    """Which HER2 residues the designed binder touches, and the epitope overlap.

    RFD4 renumbers every chain of a generated structure from 1, so antigen
    residues are shifted back to the reference (1N8Z) numbering before the
    hotspots and epitope set are applied.
    """
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    model = st[0]
    names = [c.name for c in model]
    if binder_chain not in names:
        return {"error": f"chain {binder_chain} not in {names}"}

    ref_first = _first_residue_numbers(REFERENCE_CIF)
    design_first = _first_residue_numbers(cif_path)
    offset = {
        c: ref_first.get(c, 1) - design_first.get(c, 1) for c in ANTIGEN_CHAINS if c in names
    }

    binder = _atoms(model[binder_chain])
    touched: set[tuple[str, int]] = set()
    min_hotspot = float("inf")

    for chain_name in ANTIGEN_CHAINS:
        if chain_name not in names:
            continue
        for res in model[chain_name]:
            key = (chain_name, res.seqid.num + offset[chain_name])
            for atom in res:
                d = min(atom.pos.dist(b.pos) for _, b in binder)
                if d <= CONTACT_CUTOFF:
                    touched.add(key)
                if key in HOTSPOTS:
                    min_hotspot = min(min_hotspot, d)

    epitope = epitope_residues()
    hit = touched & epitope
    return {
        "n_contacts": len(touched),
        "n_epitope_contacts": len(hit),
        "epitope_fraction": len(hit) / len(touched) if touched else 0.0,
        "min_distance_to_hotspot": None if min_hotspot == float("inf") else round(min_hotspot, 2),
        "contacts_hotspot": min_hotspot <= CONTACT_CUTOFF,
    }


def score_run(run_dir: Path | str, binder_chain: str = "H") -> pl.DataFrame:
    """Epitope check for every generated design in a run."""
    run_dir = Path(run_dir)
    rows = []
    for cif in sorted(run_dir.glob("**/inference_samples/*/*_generated.cif")):
        rows.append({"design_id": cif.parent.name, **contacts(cif, binder_chain)})
    if not rows:
        raise FileNotFoundError(f"no generated structures under {run_dir}")
    return pl.DataFrame(rows)


def main() -> None:  # pragma: no cover - thin CLI
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir")
    ap.add_argument("-o", "--output", default=None)
    ap.add_argument("--chain", default="H", help="designed binder chain")
    args = ap.parse_args()

    df = score_run(args.run_dir, args.chain)
    print(f"epitope: {len(epitope_residues())} HER2 residues within "
          f"{EPITOPE_RADIUS} A of the trastuzumab hotspots")
    print(df.sort("epitope_fraction", descending=True))
    on_target = df.filter(pl.col("contacts_hotspot")).height
    print(f"{on_target}/{df.height} designs touch a hotspot residue directly")
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.write_csv(args.output)
        print(f"-> {args.output}")


if __name__ == "__main__":  # pragma: no cover
    main()
