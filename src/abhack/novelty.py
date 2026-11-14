"""The novelty gate, mirrored locally so designs can be filtered before `abguard`.

`abguard` is the authority — always run `scripts/validate.sh` before submitting.
This module exists so a design run can be filtered in-process, cheaply.

What `abguard` actually does (verified on a probe table, 2026-11-14):

* it numbers every chain with `immunum` and compares the **IMGT CDR3**
  (``SRWGGDGFYAMDY`` for trastuzumab, the 10-residue window plus ``SR``/``Y`` flanks);
* ``identity_to_trastuzumab`` is computed over the whole domain, not the loop. A
  CDRH3-only redesign on the trastuzumab framework scores ~0.97, far above 70 %,
  so it lands in the **novel** (optimisation) class;
* that class needs **at least 3 mutations from every known CDRH3** in the
  reference, which bundles the 524k-row FACS table, the SPR table and
  trastuzumab itself. The **de_novo** class (identity < 70 % to any known
  antibody) needs 2.

The practical consequence: being 4 mutations away from trastuzumab is not enough
— the design must also be 3 mutations away from all 525,063 provided CDRH3s.
"""

from __future__ import annotations

from dataclasses import dataclass

from abhack.constants import (
    MIN_CDRH3_MUTATIONS_DE_NOVO,
    MIN_CDRH3_MUTATIONS_OPTIMISATION,
    NOVELTY_IDENTITY_THRESHOLD,
    TRASTUZUMAB_CDRH3_WINDOW,
)


def hamming(a: str, b: str) -> int:
    """Mismatches between two equal-length strings; differing lengths count fully."""
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(x != y for x, y in zip(a, b))


def levenshtein(a: str, b: str) -> int:
    """Edit distance, for comparing CDRH3s of different length."""
    if a == b:
        return 0
    if not a or not b:
        return len(a) or len(b)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        cur = [i]
        for j, cb in enumerate(b, start=1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def identity(a: str, b: str) -> float:
    """Fraction of identical positions, normalised by the longer sequence."""
    if not a or not b:
        return 0.0
    return 1.0 - levenshtein(a, b) / max(len(a), len(b))


def distance(a: str, b: str) -> int:
    """Hamming for equal lengths, edit distance otherwise."""
    return hamming(a, b) if len(a) == len(b) else levenshtein(a, b)


def mutations_from_trastuzumab(cdrh3: str) -> int:
    """Positions changed relative to ``WGGDGFYAMD``.

    Pass the 10-residue window; the IMGT flanks are constant on this framework,
    so the count is the same either way.
    """
    return distance(cdrh3, TRASTUZUMAB_CDRH3_WINDOW)


def min_mutations_to_provided(cdrh3: str, provided: frozenset[str] | set[str]) -> int:
    """Distance to the closest provided CDRH3 (exact match -> 0).

    ``provided`` holds the 10-residue windows from the FACS and SPR tables, so
    pass a window here too.
    """
    if cdrh3 in provided:
        return 0
    best = len(cdrh3)
    for known in provided:
        d = distance(cdrh3, known)
        if d < best:
            best = d
            if best == 1:
                break
    return best


@dataclass(frozen=True)
class NoveltyVerdict:
    cdrh3: str
    track: str  # "optimisation" (abguard calls it "novel") | "de_novo"
    identity_to_trastuzumab_cdrh3: float
    mutations_from_trastuzumab: int
    min_mutations_to_provided: int
    required_mutations: int
    passes: bool
    reason: str

    def as_dict(self) -> dict:
        return self.__dict__.copy()


def classify(
    cdrh3: str,
    provided: frozenset[str] | set[str] | None = None,
    *,
    domain_identity_to_known: float = 1.0,
    enforce_cdrh3_identity: bool = True,
) -> NoveltyVerdict:
    """Classify one CDRH3 window and say whether it clears its track's gate.

    ``domain_identity_to_known`` is identity of the whole variable domain to the
    nearest known antibody. It defaults to 1.0 — the right assumption for a
    CDRH3 redesign on trastuzumab's framework, which always lands in the
    optimisation track. Pass the measured value for a genuinely new scaffold.

    ``enforce_cdrh3_identity`` adds the task description's "< 70 % CDRH3 identity
    to trastuzumab" rule, which needs 4 of the 10 window positions changed. The
    shipped `abguard` (checked 2026-11-14) enforces only the mutation counts, so
    set it to False to see exactly what `abguard` would accept.
    """
    provided = provided if provided is not None else frozenset()
    muts_tz = mutations_from_trastuzumab(cdrh3)
    cdr_identity = identity(cdrh3, TRASTUZUMAB_CDRH3_WINDOW)

    if domain_identity_to_known < NOVELTY_IDENTITY_THRESHOLD:
        track, required = "de_novo", MIN_CDRH3_MUTATIONS_DE_NOVO
    else:
        track, required = "optimisation", MIN_CDRH3_MUTATIONS_OPTIMISATION

    # Trastuzumab is part of the provided data; count it explicitly so the
    # verdict is right even when `provided` has not been loaded.
    min_provided = min(
        muts_tz,
        min_mutations_to_provided(cdrh3, provided) if provided else muts_tz,
    )

    reasons = []
    if (
        enforce_cdrh3_identity
        and track == "optimisation"
        and cdr_identity >= NOVELTY_IDENTITY_THRESHOLD
    ):
        reasons.append(f"CDRH3 identity to trastuzumab {cdr_identity:.0%} >= 70%")
    if min_provided < required:
        reasons.append(
            f"only {min_provided} mutation(s) from the nearest provided CDRH3, need {required}"
        )

    return NoveltyVerdict(
        cdrh3=cdrh3,
        track=track,
        identity_to_trastuzumab_cdrh3=cdr_identity,
        mutations_from_trastuzumab=muts_tz,
        min_mutations_to_provided=min_provided,
        required_mutations=required,
        passes=not reasons,
        reason="; ".join(reasons) or "ok",
    )


def unpaired_cysteines(sequence: str) -> bool:
    """True if the chain has an odd number of cysteines (a hard-filter failure)."""
    return sequence.count("C") % 2 == 1
