# Repository Guidelines

## Project Structure & Module Organization

`src/abhack/` contains data, novelty, ranking, and submission logic. Shell
entry points live in `scripts/`. Candidate tables and structures are in
`designs/`; disposable run outputs belong in `runs/`. Tests are in `tests/`.

## Build, Test, and Development Commands

Run from the repository root:

```bash
uv sync
uv run pytest
uv run ruff check src tests
bash scripts/smoke_test.sh
```

For optional thermal MD: `uv sync --extra md`, then run
`bash scripts/run_all_thermal.sh designs/structures/*_boltz2.pdb`.

## Coding Style & Naming Conventions

Use Python 3.12, four-space indentation, type hints on public functions, and
clear `snake_case` names. Ruff enforces 100 columns. Keep columns explicit and
stable: join data on `design_id`, suffix model metrics (for example,
`min_ipae_boltz2`), and never use CSV row order as a ranking.

## Testing Guidelines

Add focused `pytest` coverage for behavior changes. Use small `polars` tables
for rank and join tests. Normal tests must not require a GPU. GPU execution
must fail clearly without CUDA; it must never silently fall back to CPU.

## Commit & Pull Request Guidelines

Use concise imperative commits, e.g. `Add thermal shortlist gate`. Do not
commit trajectories, checkpoints, credentials, or generated run output. In a
PR, state the ranking or biological decision affected, commands run, and any
candidate-count or metric change. `abguard` is authoritative for eligibility.
