#!/usr/bin/env python
"""Validate antibody designs against HER2 with OpenDDE-abag on JapanFold.

A second opinion on the optimisation track, from a model built for
antibody-antigen complexes, running off-box so it does not queue behind the
local GPU. Run **unguided**: no epitope constraint, so the prediction is
evidence about where the antibody binds rather than a restatement of what we
asked for.

Each design goes in as its own target: HER2 domain IV plus the design's VH and
VL. MSA is on by default for this model, which the local Boltz-2 refold lacked.

    python scripts/japanfold_abag.py designs/opt-s11.ranked.csv --limit 20 --submit
    python scripts/japanfold_abag.py --poll <job_id>
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

import polars as pl

BASE = os.environ.get("JAPANFOLD_BASE_URL", "https://api.japanfold.aiand.com")
KEY_FILE = Path.home() / "japanfold" / "japanfold.key"

#: HER2 domain IV, 1N8Z residues 489-607 as one continuous chain.
HER2 = (
    "CHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTCFGPE"
    "ADQCVACAHYKDPPFCVARCPSGVKPDLSYMPIWKFPDEEGACQPCPIN"
)


def _key() -> str:
    if os.environ.get("JAPANFOLD_API_KEY"):
        return os.environ["JAPANFOLD_API_KEY"]
    if KEY_FILE.exists():
        return KEY_FILE.read_text().strip()
    raise SystemExit(f"no key: set JAPANFOLD_API_KEY or write {KEY_FILE}")


def _call(method: str, path: str, payload: dict | None = None) -> dict:
    cmd = ["curl", "-s", "--max-time", "180",
           "-H", f"Authorization: Bearer {_key()}", "-H", "X-JapanFold-Client: skill"]
    tmp_name = None
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


def report_credits(label: str) -> float:
    c = _call("GET", "/v1/credits")
    print(f"[credits] {label}: ${c['balance']:,.2f}")
    return c["balance"]


def target_yaml(vh: str, vl: str | None) -> str:
    chains = [f"  - protein: {{id: A, sequence: {HER2}}}",
              f"  - protein: {{id: H, sequence: {vh}}}"]
    if vl:
        chains.append(f"  - protein: {{id: L, sequence: {vl}}}")
    return "sequences:\n" + "\n".join(chains) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("table", nargs="?")
    ap.add_argument("--limit", type=int, default=20)
    ap.add_argument("--rank-by", default="p_high")
    ap.add_argument("--model", default="opendde-abag")
    ap.add_argument("--name", default="abag")
    ap.add_argument("--submit", action="store_true")
    ap.add_argument("--poll")
    ap.add_argument("--params", default=None, help='JSON, e.g. {"write_pae": true, "output_format": "pdb"}')
    ap.add_argument("--fetch", default=None, help="with --poll: download the job archive into this directory")
    args = ap.parse_args()

    if args.poll:
        job = _call("GET", f"/v1/jobs/{args.poll}?wait=60")
        print(json.dumps({k: job.get(k) for k in ("id", "status", "queue", "error")}, indent=1))
        if job.get("status") == "succeeded":
            res = _call("GET", f"/v1/jobs/{args.poll}/results")
            rows = res.get("rows", [])
            rows.sort(key=lambda r: -(r.get("iptm") or 0))
            print(f"{'design':18} {'iptm':>6} {'plddt':>7}")
            for r in rows:
                print(f"{r.get('id',''):18} {r.get('iptm',0):6.3f} {r.get('complex_plddt',0):7.3f}")
        if args.fetch and job.get("status") == "succeeded":
            import zipfile
            out = Path(args.fetch)
            out.mkdir(parents=True, exist_ok=True)
            subprocess.run(["curl", "-fsL", "-o", str(out / "archive.zip"), "-H", f"Authorization: Bearer {_key()}",
                            "-H", "X-JapanFold-Client: skill", f"{BASE}/v1/jobs/{args.poll}/archive"], check=True)
            zipfile.ZipFile(out / "archive.zip").extractall(out)
            print(f"fetched -> {out}")
        report_credits("after poll")
        return 0

    if not args.table:
        raise SystemExit("table required unless --poll")
    df = pl.read_csv(args.table, infer_schema_length=10_000)
    if args.rank_by in df.columns:
        df = df.sort(args.rank_by, descending=True, nulls_last=True)
    df = df.head(args.limit)

    targets = [
        {"name": r["design_id"], "content": target_yaml(r["vh_sequence"], r.get("vl_sequence"))}
        for r in df.to_dicts()
    ]
    print(f"{len(targets)} targets from {args.table} (ranked by {args.rank_by})")
    if not args.submit:
        print("dry run; pass --submit to send")
        return 0

    report_credits("before submit")
    payload = {"model": args.model, "name": args.name, "targets": targets}
    if args.params:
        payload["params"] = json.loads(args.params)
    job = _call("POST", "/v1/predictions", payload)
    if "id" not in job:
        raise SystemExit(f"submit refused: {json.dumps(job, indent=1)[:900]}")
    print(f"job {job['id']} status={job.get('status')}")
    print(f"poll: python scripts/japanfold_abag.py --poll {job['id']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
