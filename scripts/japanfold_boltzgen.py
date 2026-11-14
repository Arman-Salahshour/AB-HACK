#!/usr/bin/env python
"""Design HER2 nanobodies with BoltzGen on JapanFold, aimed at trastuzumab's epitope.

A second, independent shot at the de novo track: a different model family from
RFD4-Proteina, running off-box so it does not compete with the local GPU.
BoltzGen generates, inverse-folds, refolds, scores and filters in one job, and
returns a ranked set with confidence metrics.

The target is HER2 domain IV (1N8Z, residues 489-607 as one continuous chain).
`binding_types` aims the binder at the trastuzumab contact face: HER2 572/573
and 591/593 are spec positions 84/85 and 103/105.

    python scripts/japanfold_boltzgen.py --num-designs 50 --budget 20 --submit
    python scripts/japanfold_boltzgen.py --poll <job_id>
    python scripts/japanfold_boltzgen.py --credits
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

BASE = os.environ.get("JAPANFOLD_BASE_URL", "https://api.japanfold.aiand.com")
KEY_FILE = Path.home() / "japanfold" / "japanfold.key"

HER2_DOMAIN_IV = (
    "CHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTCFGPE"
    "ADQCVACAHYKDPPFCVARCPSGVKPDLSYMPIWKFPDEEGACQPCPIN"
)
#: Trastuzumab contact face, 1-based in the sequence above (HER2 572-593).
EPITOPE = "84..105"
#: VHH length range to design.
BINDER_LENGTH = "110..130"


def _key() -> str:
    if os.environ.get("JAPANFOLD_API_KEY"):
        return os.environ["JAPANFOLD_API_KEY"]
    if KEY_FILE.exists():
        return KEY_FILE.read_text().strip()
    raise SystemExit(f"no key: set JAPANFOLD_API_KEY or write {KEY_FILE}")


def _call(method: str, path: str, payload: dict | None = None) -> dict:
    cmd = ["curl", "-s", "--max-time", "120",
           "-H", f"Authorization: Bearer {_key()}", "-H", "X-JapanFold-Client: skill"]
    tmp = None
    if payload is not None:
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


def credits() -> dict:
    return _call("GET", "/v1/credits")


def report_credits(label: str) -> float:
    c = credits()
    bal = c["balance"]
    print(f"[credits] {label}: ${bal:,.2f} at {c['rate']['description']}")
    return bal


def spec() -> str:
    return (
        "entities:\n"
        f"  - protein: {{id: A, sequence: {HER2_DOMAIN_IV}, "
        f"binding_types: {{binding: '{EPITOPE}'}}}}\n"
        f"  - protein: {{id: B, sequence: {BINDER_LENGTH}}}\n"
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--num-designs", type=int, default=50)
    ap.add_argument("--budget", type=int, default=20)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--protocol", default="nanobody-anything")
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--poll")
    ap.add_argument("--credits", action="store_true")
    args = ap.parse_args()

    if args.credits:
        report_credits("balance")
        return 0

    if args.poll:
        job = _call("GET", f"/v1/jobs/{args.poll}?wait=60")
        print(json.dumps({k: job.get(k) for k in ("id", "status", "queue", "error", "partial")},
                         indent=1))
        if job.get("status") == "succeeded":
            res = _call("GET", f"/v1/jobs/{args.poll}/results")
            for d in res.get("designs", []):
                print(d)
        report_credits("after poll")
        return 0

    print("spec:\n" + spec())
    if not args.submit:
        print("dry run; pass --submit to send")
        return 0

    before = report_credits("before submit")
    params = {"num_designs": args.num_designs, "budget": args.budget}
    if args.seed is not None:
        params["seed"] = args.seed
    job = _call("POST", "/v1/designs", {
        "model": "boltzgen", "protocol": args.protocol, "params": params, "spec": spec(),
    })
    if "id" not in job:
        raise SystemExit(f"submit refused: {json.dumps(job, indent=1)[:900]}")
    print(f"job {job['id']} status={job.get('status')}")
    print(f"poll: python scripts/japanfold_boltzgen.py --poll {job['id']}")
    after = report_credits("after submit")
    print(f"[credits] reserved so far: ${before - after:,.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
