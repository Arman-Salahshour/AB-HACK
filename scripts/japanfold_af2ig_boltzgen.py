#!/usr/bin/env python
"""Score BoltzGen de novo designs with AF2-IG on JapanFold.

BoltzGen's own ranking disagrees with its interface confidence (its rank-1
design is not its best-ipTM design), and its refold is the same family of
evidence that mis-scored the RFD4 VHHs. AF2-IG, starting from the designed
coordinates, is the orthogonal check that decided the RFD4 question, so the
BoltzGen candidates get the same treatment before any of them takes a
submission slot.

BoltzGen CIFs are whole complexes: chain A is HER2 domain IV (one chain
already — nothing to merge, unlike the RFD4 spec's C/D split) and chain B
is the designed binder.

    python scripts/japanfold_af2ig_boltzgen.py --submit
    python scripts/japanfold_af2ig_boltzgen.py --poll <job_id>
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import uuid
from pathlib import Path

import gemmi

BASE = os.environ.get("JAPANFOLD_BASE_URL", "https://api.japanfold.aiand.com")
KEY_FILE = Path.home() / "japanfold" / "japanfold.key"

MANIFEST = Path(__file__).resolve().parents[1] / "runs" / "boltzgen" / "manifest.json"


def _key() -> str:
    if os.environ.get("JAPANFOLD_API_KEY"):
        return os.environ["JAPANFOLD_API_KEY"]
    if KEY_FILE.exists():
        return KEY_FILE.read_text().strip()
    raise SystemExit(f"no key: set JAPANFOLD_API_KEY or write {KEY_FILE}")


def _curl(method: str, path: str, payload: dict | None = None) -> dict:
    cmd = ["curl", "-s", "--max-time", "120",
           "-H", f"Authorization: Bearer {_key()}",
           "-H", "X-JapanFold-Client: skill"]
    tmp_name: str | None = None
    if payload is not None:
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as tmp:
            json.dump(payload, tmp)
        tmp_name = tmp.name
        cmd += ["-X", method, "-H", "Content-Type: application/json",
                "-H", f"Idempotency-Key: {uuid.uuid4()}", "-d", f"@{tmp_name}"]
    cmd.append(f"{BASE}{path}")
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    finally:
        if tmp_name is not None:
            os.unlink(tmp_name)
    try:
        return json.loads(out)
    except json.JSONDecodeError:
        raise SystemExit(f"unexpected reply from {path}: {out[:400]}")


def target_from_complex(cif_path: Path) -> dict:
    """AF2-IG target: HER2 chain A + binder chain B, straight from the CIF."""
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    model = st[0]
    out = gemmi.Structure()
    out.cell, out.spacegroup_hm = st.cell, st.spacegroup_hm
    m = gemmi.Model("1")
    for name in ("A", "B"):
        if name not in [c.name for c in model]:
            raise ValueError(f"chain {name} not in {cif_path}")
        ch = gemmi.Chain(name)
        for i, res in enumerate(model[name], start=1):
            ch.add_residue(res)
            ch[len(ch) - 1].seqid = gemmi.SeqId(i, " ")
        m.add_chain(ch)
    out.add_model(m)
    out.setup_entities()
    seq = gemmi.one_letter_code([r.name for r in model["B"]]).upper()
    pdb = out.make_pdb_string()
    content = (
        "target:\n  structure: |\n"
        + "".join(f"    {line}\n" for line in pdb.splitlines())
        + "  chain: A\nbinder:\n"
        + f"  sequence: {seq}\n  chain: B\n"
    )
    return {"name": cif_path.stem, "content": content, "binder_seq": seq}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--poll")
    args = ap.parse_args()

    if args.poll:
        job = _curl("GET", f"/v1/jobs/{args.poll}?wait=60")
        print(json.dumps({k: job.get(k) for k in ("id", "status", "queue", "error")}, indent=1))
        if job.get("status") == "succeeded":
            res = _curl("GET", f"/v1/jobs/{args.poll}/results")
            for row in res.get("rows", []):
                print(row)
        return 0

    manifest = json.loads(MANIFEST.read_text())
    targets = []
    for d in manifest:
        t = target_from_complex(Path(d["file"]))
        targets.append({"name": t["name"], "content": t["content"]})
        print(f"  {t['name']}: binder {len(t['binder_seq'])} aa, "
              f"boltzgen iptm {d['boltzgen_iptm']}")
    if not targets:
        raise SystemExit("no designs in manifest")
    if not args.submit:
        print("dry run; pass --submit to send")
        return 0

    job = _curl("POST", "/v1/predictions",
                {"model": "af2ig", "name": "af2ig-boltzgen", "targets": targets})
    if "id" not in job:
        raise SystemExit(f"submit failed: {job}")
    print(f"job {job['id']} status={job.get('status')}")
    print(f"poll: python scripts/japanfold_af2ig_boltzgen.py --poll {job['id']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
