"""Loaders for the hackathon tables, with a parquet cache under data/processed/.

    python -m abhack.data cache     # build the caches, print a summary
"""

from __future__ import annotations

import zipfile
from functools import lru_cache
from pathlib import Path

import polars as pl

from abhack.constants import (
    DATA_PROCESSED,
    FACS_CSV,
    SPR_CSV,
    TRASTUZUMAB_CDRH3_WINDOW,
)

FACS_PARQUET = DATA_PROCESSED / "facs.parquet"
SPR_PARQUET = DATA_PROCESSED / "spr.parquet"
PROVIDED_PARQUET = DATA_PROCESSED / "provided_cdrh3.parquet"

#: CDRH3 of the one non-trastuzumab antibody in the SPR table; drop it when
#: modelling trastuzumab variants (different VH and VL, no KD).
OTHER_ANTIBODY_CDRH3 = "ARDRIVDDYYYYGMDV"


def _read_csv_maybe_zipped(path: Path) -> pl.DataFrame:
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as zf:
            names = [
                n for n in zf.namelist() if n.endswith(".csv") and not n.startswith("__MACOSX")
            ]
            if len(names) != 1:
                raise ValueError(f"expected one CSV in {path}, found {names}")
            (name,) = names
            with zf.open(name) as fh:
                return pl.read_csv(fh.read())
    return pl.read_csv(path)


def load_facs(*, use_cache: bool = True) -> pl.DataFrame:
    """524,346 trastuzumab CDRH3 variants with a FACS class (high/medium/low)."""
    if use_cache and FACS_PARQUET.exists():
        return pl.read_parquet(FACS_PARQUET)
    df = _read_csv_maybe_zipped(FACS_CSV)
    DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
    df.write_parquet(FACS_PARQUET)
    return df


def load_spr(*, use_cache: bool = True, trastuzumab_only: bool = False) -> pl.DataFrame:
    """768 SPR-measured designs; ``kd_nM`` is set for 190 of them.

    With ``trastuzumab_only`` the rows of the unrelated antibody are dropped.
    """
    if use_cache and SPR_PARQUET.exists():
        df = pl.read_parquet(SPR_PARQUET)
    else:
        df = _read_csv_maybe_zipped(SPR_CSV)
        DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
        df.write_parquet(SPR_PARQUET)
    if trastuzumab_only:
        df = df.filter(pl.col("cdrh3") != OTHER_ANTIBODY_CDRH3)
    return df


def load_spr_quantified(**kwargs) -> pl.DataFrame:
    """Only the rows with a measured KD, with ``pkd = -log10(KD in M)`` added."""
    df = load_spr(trastuzumab_only=True, **kwargs)
    return df.filter(pl.col("kd_nM").is_not_null()).with_columns(
        pkd=-(pl.col("kd_nM") * 1e-9).log10()
    )


def trastuzumab_baseline_kd() -> dict[str, float]:
    """KD statistics of the trastuzumab reference rows in the SPR table."""
    ref = load_spr().filter(
        (pl.col("cdrh3") == TRASTUZUMAB_CDRH3_WINDOW) & pl.col("kd_nM").is_not_null()
    )
    kd = ref["kd_nM"]
    return {"n": len(kd), "min": kd.min(), "median": kd.median(), "max": kd.max()}


@lru_cache(maxsize=1)
def provided_cdrh3() -> frozenset[str]:
    """Every CDRH3 in the provided data — what the novelty gate compares against.

    The union of the FACS and SPR tables (about 525k unique sequences).
    """
    if PROVIDED_PARQUET.exists():
        return frozenset(pl.read_parquet(PROVIDED_PARQUET)["cdrh3"].to_list())
    seqs = set(load_facs()["cdrh3"].to_list()) | set(load_spr()["cdrh3"].to_list())
    DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
    pl.DataFrame({"cdrh3": sorted(seqs)}).write_parquet(PROVIDED_PARQUET)
    return frozenset(seqs)


def known_binders(level: str = "high") -> frozenset[str]:
    """CDRH3s the FACS screen put in one class. ``level`` is high/medium/low."""
    df = load_facs().filter(pl.col("binding_class") == level)
    return frozenset(df["cdrh3"].to_list())


def main() -> None:  # pragma: no cover - thin CLI
    import sys

    if len(sys.argv) > 1 and sys.argv[1] != "cache":
        print(__doc__)
        raise SystemExit(1)
    facs, spr = load_facs(), load_spr()
    classes = dict(facs["binding_class"].value_counts().iter_rows())
    print(f"FACS   {facs.height:>7,} rows  classes={classes}")
    print(f"SPR    {spr.height:>7,} rows  with KD={spr['kd_nM'].is_not_null().sum()}")
    print(f"provided unique CDRH3  {len(provided_cdrh3()):,}")
    print(f"trastuzumab baseline KD (nM): {trastuzumab_baseline_kd()}")
    print(f"caches in {DATA_PROCESSED}")


if __name__ == "__main__":  # pragma: no cover
    main()
