#!/usr/bin/env python
"""Score designed complexes with AF2-IG on JapanFold.

AlphaFold2 initial-guess starts from the designed coordinates and asks whether
the complex is self-consistent, which is the field-standard check for a binder.
It is a fairer test than refolding from sequence alone: our local Boltz-2 run
had no MSA on any chain and had to rediscover the pose from scratch.

AF2-IG wants one target chain plus one binder chain, so the HER2 chains C and D
(one protein, split by the organisers at its unresolved 581-590 loop) are merged
into a single chain A, and the designed VHH becomes chain B.

    python scripts/japanfold_af2ig.py runs/denovo-vhh-s1 --limit 20 --submit

The key is read from a file, never printed. Jobs are async: this submits and
prints the job id; poll it with --poll <job_id>.
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

TARGET_CHAINS = ("C", "D")  # HER2, merged into one chain for AF2-IG
BINDER_CHAIN = "H"


def _key() -> str:
    if os.environ.get("JAPANFOLD_API_KEY"):
        return os.environ["JAPANFOLD_API_KEY"]
    if KEY_FILE.exists():
        return KEY_FILE.read_text().strip()
    raise SystemExit(f"no key: set JAPANFOLD_API_KEY or write {KEY_FILE}")


def _curl(method: str, path: str, payload: dict | None = None) -> dict:
    cmd = [
        "curl", "-s", "--max-time", "120",
        "-H", f"Authorization: Bearer {_key()}",
        "-H", "X-JapanFold-Client: skill",
    ]
    tmp = None
    if payload is not None:
        # A batch of structures is far past the shell's argv limit, so the body
        # goes through a file rather than the command line.
        tmp = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
        json.dump(payload, tmp)
        tmp.close()
        cmd += ["-X", method, "-H", "Content-Type: application/json",
                "-H", f"Idempotency-Key: {uuid.uuid4()}", "-d", f"@{tmp.name}"]
    cmd.append(f"{BASE}{path}")
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    finally:
        if tmp is not None:
            os.unlink(tmp.name)
    try:
        return json.loads(out)
    except json.JSONDecodeError:
        raise SystemExit(f"unexpected reply from {path}: {out[:400]}")


def build_target(cif_path: Path) -> tuple[str, str]:
    """Return (PDB text with merged target chain A + binder chain B, binder sequence)."""
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    model = st[0]

    out = gemmi.Structure()
    out.cell, out.spacegroup_hm = st.cell, st.spacegroup_hm
    m = gemmi.Model("1")

    target = gemmi.Chain("A")
    n = 0
    for name in TARGET_CHAINS:
        if name not in [c.name for c in model]:
            continue
        for res in model[name]:
            n += 1
            target.add_residue(res)
            target[len(target) - 1].seqid = gemmi.SeqId(n, " ")
    m.add_chain(target)

    binder = gemmi.Chain("B")
    for i, res in enumerate(model[BINDER_CHAIN], start=1):
        binder.add_residue(res)
        binder[len(binder) - 1].seqid = gemmi.SeqId(i, " ")
    m.add_chain(binder)

    out.add_model(m)
    out.setup_entities()
    seq = gemmi.one_letter_code([r.name for r in binder]).upper()
    return out.make_pdb_string(), seq


def targets_from_run(run_dir: Path, limit: int | None, only: set[str] | None) -> list[dict]:
    targets = []
    for cif in sorted(run_dir.glob("**/inference_samples/*/*_generated.cif")):
        design_id = cif.parent.name
        if only and design_id not in only:
            continue
        pdb_text, seq = build_target(cif)
        content = (
            "target:\n  structure: |\n"
            + "".join(f"    {line}\n" for line in pdb_text.splitlines())
            + "  chain: A\nbinder:\n"
            + f"  sequence: {seq}\n  chain: B\n"
        )
        targets.append({"name": design_id, "content": content})
        if limit and len(targets) >= limit:
            break
    return targets


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir", nargs="?")
    ap.add_argument("--limit", type=int, default=20)
    ap.add_argument("--only", help="comma-separated design_ids to score")
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--poll", help="job id to check instead of submitting")
    ap.add_argument("--name", default="af2ig")
    args = ap.parse_args()

    if args.poll:
        job = _curl("GET", f"/v1/jobs/{args.poll}?wait=60")
        print(json.dumps({k: job.get(k) for k in ("id", "status", "queue", "error")}, indent=1))
        if job.get("status") == "succeeded":
            res = _curl("GET", f"/v1/jobs/{args.poll}/results")
            for row in res.get("rows", []):
                print(row)
        return 0

    if not args.run_dir:
        raise SystemExit("run_dir required unless --poll")
    only = set(args.only.split(",")) if args.only else None
    targets = targets_from_run(Path(args.run_dir), args.limit, only)
    print(f"{len(targets)} targets prepared from {args.run_dir}")
    if not targets:
        return 1
    print(f"  example: {targets[0]['name']}, binder "
          f"{len(targets[0]['content'].split('sequence: ')[1].split()[0])} aa")

    if not args.submit:
        print("dry run; pass --submit to send")
        return 0

    job = _curl("POST", "/v1/predictions",
                {"model": "af2ig", "name": args.name, "targets": targets})
    if "id" not in job:
        raise SystemExit(f"submit failed: {job}")
    print(f"job {job['id']} status={job.get('status')}")
    print(f"poll: python scripts/japanfold_af2ig.py --poll {job['id']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
