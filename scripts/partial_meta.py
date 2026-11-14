#!/usr/bin/env python3
"""Write metadata for MD runs that were stopped before finishing, so they can be analysed.

The run scripts write their metadata only at the end. For a killed run this reads the production time
actually reached from the StateDataReporter log (production starts at frame 0 of the DCD) and writes a
metadata file marked "partial": true, so no partial result is mistaken for a finished one.

    python scripts/partial_meta.py runs/thermal-poc runs/phmd
"""
import json
import re
import sys
from pathlib import Path

import mdtraj as md


def prod_ns(prefix: Path) -> float:
    """Production nanoseconds actually written: frames in the DCD x the 10 ps frame interval."""
    top = f"{prefix}_top.pdb"
    n = md.load_frame(f"{prefix}.dcd", 0, top=top) and sum(1 for _ in md.iterload(f"{prefix}.dcd", top=top, chunk=1))
    return round(n * 10 / 1000.0, 3)


for d in map(Path, sys.argv[1:]):
    for dcd in sorted(d.glob("*.dcd")):
        prefix = dcd.with_suffix("")
        m = re.search(r"(.+)_T(\d+)_r(\d+)$", prefix.name)
        p = re.search(r"(.+)_pH([\d.]+)_r(\d+)$", prefix.name)
        thermal = bool(m)
        meta_path = Path(f"{prefix}_meta.json") if thermal else Path(f"{prefix}.json")
        if meta_path.exists():
            continue
        ns = prod_ns(prefix)
        if thermal:
            meta = {"design_id": m.group(1), "temp_K": float(m.group(2)), "seed": int(m.group(3)), "ns": ns,
                    "binder_chains": ["A", "B"], "format": "fv", "partial": True}
        else:
            meta = {"ph": float(p.group(2)), "seed": int(p.group(3)), "ns": ns, "temp_K": 310.0,
                    "his_states": {}, "antigen": ["C", "D"], "partial": True}
        meta_path.write_text(json.dumps(meta, indent=1))
        print(f"{prefix.name}: {ns} ns production, partial -> {meta_path.name}")
