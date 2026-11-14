"""pH-selectivity of a designed complex from propka, by proton linkage.

For every titratable group i the bound-state pKa (``pKb``, in the HER2 complex) and
the free-state pKa (``pKf``, the Fv alone) give its contribution to the observed
association constant at a given pH:

    log10 K_obs(pH) = log10 K_0 + sum_i [ log10(1 + 10^(pKb_i - pH))
                                         - log10(1 + 10^(pKf_i - pH)) ]

so the predicted preference for tumour pH over plasma pH is

    selectivity = log10 K(pH_low) - log10 K(pH_high)          (default 6.5 vs 7.4)

Positive means the complex binds better at the lower pH. A trastuzumab control gives
the noise floor. Propka is an empirical, static-structure predictor (errors of about
one pKa unit), so this ranks candidates; it does not measure a switch.

Variants are threaded onto the 1N8Z crystal (``PDBFixer.applyMutations``: side chains
are rebuilt, the backbone is not), which is a triage model. Real structures from a
refold should replace it for the shortlist.
"""

from __future__ import annotations

import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import gemmi
from pdbfixer import PDBFixer
from openmm.app import PDBFile

from abhack.constants import TRASTUZUMAB_VH, TRASTUZUMAB_VL
from abhack.epitope import REFERENCE_CIF

PH_LOW, PH_HIGH = 6.5, 7.4  # tumour microenvironment vs plasma/heart
THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
PROPKA = str(Path(sys.executable).with_name("propka3"))
_ROW = re.compile(r"^\s*(ASP|GLU|HIS|CYS|TYR|LYS|ARG|C-|N\+)\s+(\d+)\s+(\S)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")


def mutations(vh: str, vl: str) -> list[tuple[str, int, str, str]]:
    """(chain, 1-based pos, wt, mut) differences from trastuzumab; lengths must match."""
    out = []
    for chain, new, ref in (("B", vh, TRASTUZUMAB_VH), ("A", vl, TRASTUZUMAB_VL)):
        if len(new) != len(ref):
            raise ValueError(f"chain {chain} length {len(new)} != {len(ref)} (indels unsupported)")
        out += [(chain, i + 1, w, m) for i, (w, m) in enumerate(zip(ref, new, strict=True)) if w != m]
    return out


def thread(vh: str, vl: str, out_pdb: Path) -> Path:
    """Write the HER2-domain-IV complex carrying the given variable domains (no hydrogens)."""
    st = gemmi.read_structure(str(REFERENCE_CIF))
    st.remove_ligands_and_waters()
    tmp = out_pdb.with_suffix(".in.pdb")
    st.write_pdb(str(tmp))
    fx = PDBFixer(filename=str(tmp))
    by_chain: dict[str, list[str]] = {}
    for chain, pos, wt, mut in mutations(vh, vl):
        by_chain.setdefault(chain, []).append(f"{THREE[wt]}-{pos}-{THREE[mut]}")
    for chain, muts in by_chain.items():
        fx.applyMutations(muts, chain)
    fx.findMissingResidues()
    fx.missingResidues = {}
    fx.findMissingAtoms()
    fx.addMissingAtoms()
    with open(out_pdb, "w") as f:
        PDBFile.writeFile(fx.topology, fx.positions, f, keepIds=True)
    tmp.unlink(missing_ok=True)
    return out_pdb


def _propka(pdb: Path, workdir: Path) -> dict[tuple[str, int, str], tuple[float, float]]:
    subprocess.run([PROPKA, "-q", str(pdb.name)], cwd=workdir, check=True, capture_output=True)
    text = (workdir / pdb.with_suffix(".pka").name).read_text()
    block = text.split("SUMMARY OF THIS PREDICTION", 1)[1]
    rows = {}
    for line in block.splitlines():
        m = _ROW.match(line)
        if m:
            rows[(m.group(3), int(m.group(2)), m.group(1))] = (float(m.group(4)), float(m.group(5)))
    return rows


def _apo(pdb: Path, out: Path, antigen=("C", "D")) -> Path:
    keep = [ln for ln in pdb.read_text().splitlines()
            if ln.startswith(("ATOM", "HETATM", "TER")) and ln[21] not in antigen]
    out.write_text("\n".join(keep) + "\nEND\n")
    return out


def _log_term(pka: float, ph: float) -> float:
    return math.log10(1.0 + 10.0 ** (pka - ph))


def selectivity(bound: dict, free: dict, ph_low: float = PH_LOW, ph_high: float = PH_HIGH) -> dict:
    """Proton-linkage preference for ph_low over ph_high, summed over every shared group."""
    total, parts = 0.0, []
    for key, (pk_b, _) in bound.items():
        if key not in free:
            continue  # antigen groups have no free-Fv counterpart
        pk_f = free[key][0]
        d = (_log_term(pk_b, ph_low) - _log_term(pk_f, ph_low)) - (
            _log_term(pk_b, ph_high) - _log_term(pk_f, ph_high)
        )
        total += d
        if abs(d) > 0.01:
            parts.append((key, round(pk_f, 2), round(pk_b, 2), round(d, 3)))
    parts.sort(key=lambda x: -abs(x[3]))
    return {"log10_sel": total, "fold_sel": 10.0**total, "top_groups": parts[:6]}


def score_variant(
    vh: str,
    vl: str,
    ph_low: float = PH_LOW,
    ph_high: float = PH_HIGH,
    pdb: Path | None = None,
    antigen: tuple[str, ...] = ("C", "D"),
) -> dict:
    """Thread (or use ``pdb``) -> propka complex and free Fv -> selectivity."""
    with tempfile.TemporaryDirectory(prefix="phscore_") as d:
        wd = Path(d)
        cplx = pdb if pdb is not None else thread(vh, vl, wd / "cplx.pdb")
        if pdb is not None:
            (wd / "cplx.pdb").write_text(Path(pdb).read_text())
            cplx = wd / "cplx.pdb"
        apo = _apo(cplx, wd / "apo.pdb", antigen)
        res = selectivity(_propka(cplx, wd), _propka(apo, wd), ph_low, ph_high)
    his = [(k, pf, pb) for k, pf, pb, _ in res["top_groups"] if k[2] == "HIS"]
    res["his_summary"] = ";".join(f"{k[0]}{k[1]}:{pf}->{pb}" for k, pf, pb in his)
    res["top_groups"] = ";".join(f"{k[2]}{k[1]}{k[0]}:{pf}->{pb}({d:+.2f})" for k, pf, pb, d in res["top_groups"])
    return res
