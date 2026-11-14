"""One command per pipeline stage, so an agent runs a tested path.

    python -m abhack.pipeline generate --spec data/opt-epitope/design.json --seed 21
    python -m abhack.pipeline validate --run runs/opt-epitope-s21
    python -m abhack.pipeline finalise --runs runs/opt-novel-s11 runs/opt-novel-s12 --limit 20

Every stage is idempotent: re-running `validate` on a finished run recomputes
the tables without touching the GPU. `status` reports what exists so an agent
can decide what still needs doing.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import polars as pl

from abhack.constants import ROOT

#: HER2 epitope residues contacting trastuzumab, in 1N8Z numbering.
#: Chain C is 489-580, chain D is 591-607 (split at the unresolved 581-590 loop).
#: 14 are confirmed by both 1N8Z and 7MN8; the rest by one structure.
EPITOPE = [("C", n) for n in (542, 557, 558, 560, 561, 569, 570, 571, 572, 573, 579, 580)] + \
          [("D", n) for n in (591, 593, 596, 598, 600, 601, 602, 603, 604, 605)]
CONTACT_CUTOFF = 4.5
REF_FIRST = {"C": 489, "D": 591}

#: A refold is "strong" when an independent model reproduces the design.
STRONG = {"min_ipae": 0.5, "iptm": 0.75, "backbone_rmsd": 2.0}


def sh(cmd: list[str] | str, **kw) -> int:
    """Run a shell step, streaming output; returns the exit code."""
    printable = cmd if isinstance(cmd, str) else " ".join(cmd)
    print(f"$ {printable}", flush=True)
    return subprocess.run(cmd, shell=isinstance(cmd, str), cwd=ROOT, check=False, **kw).returncode


# --------------------------------------------------------------------------- generate
def generate(spec: str, seed: int, n: int, weight: str, out: str | None) -> Path:
    out_dir = Path(out or f"runs/{Path(spec).parent.name}-s{seed}")
    env = f"WEIGHT={weight} N={n} SEED={seed} INPUT={spec} OUT={out_dir}"
    if sh(f"{env} bash scripts/design.sh"):
        raise SystemExit(f"generation failed; see {out_dir}/design/generation.log")
    return out_dir


def refold(run: str, top_n: int, fold: str, weight: str) -> None:
    if sh(f"OUT={run} TOPN={top_n} FOLD={fold} WEIGHT={weight} bash scripts/refold.sh"):
        raise SystemExit("refolding failed")


# --------------------------------------------------------------------------- epitope
def epitope_coverage(run: str, binder_chains=("A", "B")) -> pl.DataFrame:
    """Which epitope residues each design's antibody actually contacts."""
    import gemmi

    want = set(EPITOPE)
    rows = []
    for cif in sorted(Path(run).glob("**/inference_samples/*/*_generated.cif")):
        st = gemmi.read_structure(str(cif))
        st.setup_entities()
        m = st[0]
        names = [c.name for c in m]
        first = {c.name: c[0].seqid.num for c in m}
        # generated structures renumber every chain from 1
        off = {c: REF_FIRST[c] - first.get(c, 1) for c in REF_FIRST if c in names}
        ab = [a for cn in binder_chains if cn in names for r in m[cn] for a in r]
        if not ab:
            continue
        hit, tot = set(), set()
        for cn in REF_FIRST:
            if cn not in names:
                continue
            for res in m[cn]:
                key = (cn, res.seqid.num + off[cn])
                for atom in res:
                    if min(atom.pos.dist(b.pos) for b in ab) <= CONTACT_CUTOFF:
                        tot.add(key)
                        if key in want:
                            hit.add(key)
                        break
        rows.append({
            "design_id": cif.parent.name,
            "n_contacts": len(tot),
            "n_epitope": len(hit),
            "epitope_coverage": round(len(hit) / len(want), 3),
            "epitope_precision": round(len(hit) / len(tot), 3) if tot else 0.0,
        })
    if not rows:
        raise FileNotFoundError(f"no generated structures under {run}")
    return pl.DataFrame(rows)


# --------------------------------------------------------------------------- validate
def validate(run: str, out: str | None = None) -> dict:
    """Collect, score the epitope, filter novelty, and run abguard. CPU only."""
    run_path = Path(run)
    tag = run_path.name
    table = Path(out or f"designs/{tag}.csv")

    if sh([sys.executable, "-m", "abhack.designs", "collect", str(run_path), "-o", str(table)]):
        raise SystemExit("collect failed")

    df = pl.read_csv(table, infer_schema_length=10_000)
    epi = epitope_coverage(str(run_path))
    df = df.join(epi, on="design_id", how="left")
    df.write_csv(table)

    passing = df.filter(pl.col("passes") == True)
    unique = passing.unique(subset=["cdrh3"], keep="first")

    from abhack.submission import to_submission

    sub_path = Path(f"designs/{tag}-submission.csv")
    to_submission(unique, kind="paired").write_csv(sub_path)
    verdict = subprocess.run(["bash", "scripts/validate.sh", str(sub_path), "--drop",
                              "-o", f"designs/{tag}-valid.csv"],
                             cwd=ROOT, capture_output=True, text=True, check=False)
    valid = 0
    for line in (verdict.stdout + verdict.stderr).splitlines():
        if "design(s) valid" in line:
            print(line.strip())
            valid = int(line.split(":")[-1].strip().split("/")[0])

    refolded = df.filter(pl.col("min_ipae_boltz2").is_not_null()) \
        if "min_ipae_boltz2" in df.columns else df.head(0)
    strong = refolded.filter(
        (pl.col("min_ipae_boltz2") < STRONG["min_ipae"])
        & (pl.col("iptm_boltz2") > STRONG["iptm"])
        & (pl.col("backbone_rmsd_boltz2") < STRONG["backbone_rmsd"])
    ) if refolded.height else refolded

    summary = {
        "run": tag,
        "designs": df.height,
        "unique_cdrh3": df["cdrh3"].n_unique(),
        "novelty_passing": passing.height,
        "abguard_valid": valid,
        "refolded": refolded.height,
        "strong_refold": strong.height,
        "mean_epitope_coverage": round(df["epitope_coverage"].mean() or 0, 3),
        "max_epitope_contacts": int(df["n_epitope"].max() or 0),
        "table": str(table),
    }
    for k, v in summary.items():
        print(f"  {k:24} {v}")
    return summary


# --------------------------------------------------------------------------- finalise
def finalise(runs: list[str], limit: int, out: str) -> pl.DataFrame:
    """Merge validated runs into one ranked submission table."""
    frames = []
    for run in runs:
        tag = Path(run).name
        t = Path(f"designs/{tag}.csv")
        if not t.exists():
            print(f"  skipping {tag}: no table, run `validate` first")
            continue
        frames.append(pl.read_csv(t, infer_schema_length=10_000).with_columns(run=pl.lit(tag)))
    if not frames:
        raise SystemExit("nothing to finalise")

    df = pl.concat(frames, how="diagonal_relaxed")
    df = df.filter(pl.col("passes") == True)

    if "min_ipae_boltz2" in df.columns:
        df = df.filter(pl.col("min_ipae_boltz2").is_not_null())

    # Rank first, then dedupe: the same loop can be sampled more than once with
    # different poses, and the best-posed sample is the one to keep.
    df = df.with_columns(
        reproduced=(pl.col("backbone_rmsd_boltz2") < STRONG["backbone_rmsd"]).fill_null(False)
    ).sort(
        ["reproduced", "n_epitope", "min_ipae_boltz2"],
        descending=[True, True, False], nulls_last=True,
    ).unique(subset=["cdrh3"], keep="first", maintain_order=True).head(limit)

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(out)
    cols = [c for c in ("run", "design_id", "cdrh3", "min_ipae_boltz2", "iptm_boltz2",
                        "backbone_rmsd_boltz2", "n_epitope", "p_high") if c in df.columns]
    print(df.select(cols).to_pandas().to_string(index=False))
    print(f"\n{df.height} designs -> {out}")
    return df


def status() -> None:
    print(f"{'run':28} {'designs':>8} {'refolded':>9} {'table':>6}")
    for run in sorted(Path("runs").glob("*")):
        if not run.is_dir():
            continue
        n = len(list(run.glob("**/inference_samples/*/*_generated.cif")))
        r = len(list(run.glob("**/refold_metrics.csv")))
        t = "yes" if Path(f"designs/{run.name}.csv").exists() else "-"
        print(f"{run.name:28} {n:8} {r:9} {t:>6}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("generate")
    g.add_argument("--spec", required=True)
    g.add_argument("--seed", type=int, default=1)
    g.add_argument("--n", type=int, default=100)
    g.add_argument("--weight", default="her2-novel")
    g.add_argument("--out")
    g.add_argument("--refold", action="store_true", help="refold immediately after")
    g.add_argument("--top-n", type=int, default=30)

    v = sub.add_parser("validate")
    v.add_argument("--run", required=True)
    v.add_argument("-o", "--output")

    f = sub.add_parser("finalise")
    f.add_argument("--runs", nargs="+", required=True)
    f.add_argument("--limit", type=int, default=20)
    f.add_argument("-o", "--output", default="designs/submission-final.csv")

    sub.add_parser("status")
    args = ap.parse_args()

    if args.cmd == "generate":
        out = generate(args.spec, args.seed, args.n, args.weight, args.out)
        if args.refold:
            refold(str(out), args.top_n, "boltz2", args.weight)
        print(f"-> {out}")
    elif args.cmd == "validate":
        validate(args.run, args.output)
    elif args.cmd == "finalise":
        finalise(args.runs, args.limit, args.output)
    else:
        status()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
