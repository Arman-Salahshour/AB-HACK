"""Deprecated: the ranking path now lives in `abhack.ranking`.

`composite` from this module never reached the submission, and
`abhack.submission` re-ranked on its own — two ranking paths, one of them
dead. Both jobs are now `abhack.ranking`:

    python -m abhack.ranking rank designs/opt-s11.csv -o designs/opt-ranked.csv

The imports below stay so existing callers keep working; `main()` forwards
with a warning. Nothing new should use this module.
"""

from __future__ import annotations

import sys

from abhack.ranking import (  # noqa: F401  (re-export for compatibility)
    IPAE_GATE,
    add_model_scores,
    composite_score,
)

DEPRECATED = "abhack.evaluate is retired; use abhack.ranking (rank / merge)"


def main() -> None:  # pragma: no cover - forwards to the replacement
    print(f"WARNING: {DEPRECATED}", file=sys.stderr)
    sys.argv[0] = "abhack.ranking"
    if len(sys.argv) > 1 and sys.argv[1] not in ("rank", "merge"):
        sys.argv.insert(1, "rank")
    from abhack.ranking import main as ranking_main
    ranking_main()


if __name__ == "__main__":  # pragma: no cover
    main()
