# Generator agent

**Role:** turn a target structure, an epitope and a scaffold into a set of
candidate antibody designs, with everything the evaluator needs to judge them.

**Non-role:** deciding whether a design is good. The generator produces
candidates and provenance. Every quality claim belongs to the evaluator, which
uses a model that never saw the designed structure. See
[`evaluator.md`](evaluator.md).

> **Why the roles are split — this is the measured part.** 100 de novo designs scored
> **100/100** on this generator's own hotspot-obedience metric, and **0 of 20 survived**
> independent refolding. A generator that also grades its output has no structural reason
> to reject it. The prohibition below — you may not claim a design binds — is what that
> failure bought, and it is the architectural contribution of this project.

---

## Input

```jsonc
{
  "target":   {"structure": "<path>", "chains": ["C"], "crop": [489, 607]},
  "epitope":  {"hotspot_atoms": [{"chain": "C", "seqid": 593, "atom": "NZ"}, ...],
               "source": "PDB 1N8Z, atoms nearest trastuzumab CDRH3"},
  "scaffold": {"format": "vhh", "structure": "<path>", "chain": "H",
               "cdr_lengths": {"cdr1": [13, 13], "cdr2": [9, 10], "cdr3": [12, 21]}},
  "budget":   {"weight": "antibody-2026-10-30", "n_designs": 100, "seeds": [1, 2, 3]}
}
```

Any structure format is accepted on the way in; the agent converts. **The
epitope is mandatory** — it is the only thing aiming the binder, and a design
that binds the wrong patch cannot be graded on epitope preservation.

---

## Procedure

### 1. Prepare structures

```bash
python scripts/prepare_structure.py <in.pdb|in.cif> <out.cif> --chain C --range 489-607
```

Produces mmCIF with hydrogens, waters, ligands and altlocs removed, and the
polymer bookkeeping restored. **A PDB carries no `label_*` records**, so a naive
conversion yields rows the module silently rejects:

```
rejected:  ATOM 1 N N . CYS Cxp C  .   ... 489 C 1
accepted:  ATOM 1 N N . CYS C   3  489 ... 489 C 1
```

The scaffold needs one extra step: **strip residues that are not standard amino
acids**. A framework file may encode its CDRs as residues *named* for their
length (`13-13`, `9-10`, `12-21`). Those work as a direct `--input` but collide
with a contig that states lengths itself, and the run dies with `Failed to
compute leaving groups`.

### 2. Check the numbering frame

If the hotspots came from a different file than the target, confirm the frames
agree before trusting them — compare the residue sequence over the shared range,
and confirm each hotspot atom exists at its chain and seqid. A hotspot typo
produces a *silently wrong* answer: the design binds elsewhere and every
downstream number looks fine.

### 3. Write the spec

```
contig    C489-607, /0, H1-22, 13-13, H24-37, 9-10, H39-75, 12-21, H77-86
C_CRD     target backbone pinned in absolute coordinates
C_CRDGRP  scaffold rigid internally, free to dock
C_SEQ     keep target and framework sequence; the rest is invented
C_HOT     the epitope atoms
C_CATH    2.60.40 on the antibody chain
```

Chain-lettered ranges are kept; bare ranges are designed loops of that length.
Keep the spec in version control, separate from the generated input folder.

### 4. Probe before scaling

```bash
WEIGHT=antibody N=2 SEED=<s> INPUT=<spec> OUT=runs/probe bash scripts/design.sh
```

Spec errors surface in ~90 s and the only message is `error[input-data]:
Invalid RFD4-Proteina design input`. **Never queue a long run against an
unprobed spec.** If the probe fails, isolate by changing one variable at a time
against a configuration known to work — that is how the `label_seq_id` problem
above was found.

### 5. Generate

```bash
WEIGHT=antibody N=100 SEED=1 INPUT=<spec> OUT=runs/<name>-s1 bash scripts/design.sh
```

~10.5 s per VHH design on an H100. `--num-replicates` caps at 200, so larger
batches are several runs with different seeds. Run seeds sequentially: they
contend for the same GPU.

### 6. Parse and self-check (CPU, run during generation)

```bash
python -m abhack.designs collect runs/<name>-s1 -o designs/<name>-s1.csv
python -m abhack.epitope  runs/<name>-s1 -o designs/epitope-s1.csv
```

Report, per batch: unique sequences, length range, CDR3 length distribution,
unpaired cysteines, identity to the reference antibody (which decides the
track), and the epitope contact fraction.

---

## Output

A candidate set the evaluator can consume without re-reading structures:

| Field | Notes |
| --- | --- |
| `design_id`, `run_id`, `seed`, `weight` | provenance |
| `format`, `heavy_sequence`, `light_sequence` | the molecule |
| `cdr1`, `cdr2`, `cdr3_imgt` | IMGT loops, read with `immunum` |
| `total_length`, `unpaired_cys` | hard-filter inputs |
| `identity_to_reference`, `track` | derived, never chosen |
| `n_contacts`, `epitope_fraction`, `contacts_hotspot` | self-reported geometry |
| `structure_path` | the docked complex |

Plus the run's own `per_example_metrics.csv`: `hotspot_obedience`,
`clash_metrics`, `pre_oracle_pass`, interface contact counts.

---

## The one thing to be honest about

**Every metric the generator emits is self-assessment.** `hotspot_obedience`
answers "did I build what I was told to build", not "will this bind". The model
always emits something plausible.

Measured on this pipeline (2026-11-14): 100 de novo VHH designs, 100/100
touching a hotspot, median epitope fraction 1.00, every one of them clean. Of
the top 20 by the model's own score, **0 survived independent refolding** —
binder pLDDT 89–93 (they fold fine) but 15–50 Å from the designed pose (they
dock somewhere else entirely).

So the generator must report its numbers as *geometry of the proposal*, never as
evidence of binding, and must not rank its output as though it were validated.

---

## Failure modes

| Symptom | Cause | Fix |
| --- | --- | --- |
| `Invalid RFD4-Proteina design input` | missing `label_seq_id`/entity after PDB conversion | `prepare_structure.py` |
| `Failed to compute leaving groups for '12-21'` | CDR length placeholders left in scaffold | strip non-standard residues |
| designs bind the wrong patch | hotspot atom wrong, or numbering frame mismatch | verify atoms exist before running |
| `inference_metrics/` missing | run not finished — it is written at the end | wait |
| sequences come back null | CDR search used framework anchors the scaffold lacks | number with `immunum` |

---

## Handoff

Give the evaluator: the candidate table, the run folders, the spec, and the
epitope definition. State the batch's self-reported geometry and say plainly
that it is unvalidated.
