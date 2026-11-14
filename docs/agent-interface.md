# Designing an antibody-design agent: inputs and outputs

Notes for a future agent that drives this pipeline end to end, written from what
actually broke while building it by hand (Oslo hackathon, November 2026).

The job: **given a target structure and an epitope, return a ranked list of
antibody designs that a wet lab can order**, each with evidence for why it is on
the list. Everything below is the contract that makes that automatable.

---

## 1. The shape of the contract

```
                    ┌───────────────────────────────┐
  design_request →  │  generate → refold → validate │  → design_bundle
     (JSON)         └───────────────────────────────┘      (table + structures
                                                             + report)
```

One request in, one bundle out. Every intermediate (`runs/…`) is cache, not
product: an agent must be able to delete it and rebuild from the request alone.

---

## 2. Input: the design request

```jsonc
{
  "target": {
    "structure": "data/denovo-vhh/target_her2.cif",  // mmCIF ONLY — .pdb is rejected
    "chains": ["C", "D"],                             // what to keep as antigen
    "numbering": "author"                             // the frame hotspots are given in
  },
  "epitope": {
    // REQUIRED. Without it the binder lands wherever it likes, and a design that
    // binds the wrong patch scores zero on competitive binding.
    "hotspot_atoms": [
      {"chain": "D", "seqid": 593, "atom": "NZ"},
      {"chain": "C", "seqid": 573, "atom": "CD1"},
      {"chain": "D", "seqid": 591, "atom": "CG2"},
      {"chain": "C", "seqid": 572, "atom": "CD"}
    ],
    "reference_complex": "PDB 1N8Z",   // where the hotspots came from, for provenance
    "radius": 10.0                      // widen to a patch for scoring contacts
  },
  "scaffold": {
    "format": "vhh",                   // vhh | fv | existing_antibody
    "structure": "data/denovo-vhh/framework_vhh.cif",
    "chain": "H",
    "cdr_lengths": {"cdr1": [13, 13], "cdr2": [9, 10], "cdr3": [12, 21]}
  },
  "budget": {
    "weight": "antibody-2026-10-30",
    "n_designs": 100,
    "seeds": [1, 2, 3],
    "fold_models": ["boltz2"],
    "refold_top_n": 20
  },
  "constraints": {
    "max_total_length": 250,
    "require_paired_cysteines": true,
    "novelty_reference": "abguard-bundled"
  }
}
```

### What each input is really for

| Input | Why the agent cannot skip it |
| --- | --- |
| `target.structure` | the antigen surface. Must be mmCIF; the module rejects `.pdb` outright |
| `epitope.hotspot_atoms` | the **only** thing aiming the binder. Epitope preservation is a judged property, so a run without hotspots is a run whose output cannot be graded |
| `scaffold.structure` | the framework whose fold is trusted. For a de novo design it also decides the *track*: identity to the reference antibody must fall under 70 % |
| `scaffold.cdr_lengths` | the search space. Wider ranges = more diversity, less convergence |
| `budget` | generation is cheap (~10 s/design), refolding is not (~45 s). The split between them is the main cost decision |

### Preconditions the agent must check before spending GPU

1. **mmCIF, not PDB.** Convert if needed; do not pass `.pdb` through.
2. **Strip non-standard residues from the scaffold.** A framework file may encode
   its CDRs as residues *named* for their length (`13-13`, `9-10`, `12-21`).
   That works as a direct `--input`, but collides with a contig that states
   lengths itself, and the run dies with `Failed to compute leaving groups`.
3. **Hotspot atoms must exist** in the target file, by chain and seqid, in the
   numbering the request declares. A typo here produces a silent wrong answer,
   not an error — the design simply binds elsewhere.
4. **Chain names must not collide** between target and scaffold files.
5. **Run a 2-design probe first.** Spec errors surface in ~90 s. Never queue
   hundreds of designs against an unvalidated spec.

---

## 3. Output: the design bundle

### 3.1 Per-design record

The unit of output. One row per design, stable across formats so tables from
different runs concatenate:

| Field | Meaning |
| --- | --- |
| `design_id` | unique within a run (`her2_vhh_denovo_42`) |
| `run_id`, `seed`, `weight` | provenance — which configuration produced it |
| `format` | `vhh` / `paired` / `scfv` |
| `heavy_sequence`, `light_sequence` | the molecule itself |
| `cdr1`, `cdr2`, `cdr3_imgt` | IMGT loops, read with `immunum` |
| `total_length`, `unpaired_cys` | hard-filter inputs |
| `identity_to_reference` | decides the track |
| `track` | `de_novo` / `optimisation` — **derived, never chosen** |
| `min_ipae`, `iptm`, `binder_plddt`, `rmsd_to_design`, `pass` | refold evidence, suffixed per folding model |
| `n_contacts`, `epitope_fraction`, `contacts_hotspot` | did it bind the right patch |
| `validator_verdict`, `validator_reason` | the authority's answer, not ours |
| `rank`, `rank_rationale` | why it is where it is — **prose, not just a number** |
| `structure_path` | the docked complex |

`rank_rationale` is not decoration. A human has to defend the list; a bare score
is not a defence.

### 3.2 Structures

Each design is a full **docked complex** (antigen chains + binder), one model,
mmCIF. Note that generated structures **renumber every chain from 1** — HER2
C489 becomes residue 1. Any comparison against the reference numbering must
shift, or it silently finds zero epitope contacts.

### 3.3 Report

Prose, for the humans deciding what to order: the request, how many designs
survived each stage, the ranked table, what the agent is uncertain about, and
what it would do with more budget. Attrition per stage is the most informative
thing in it — "300 generated → 280 numbered → 27 refolded → 24 validated".

---

## 4. Stage contracts

| Stage | In | Out | Fails when |
| --- | --- | --- | --- |
| **build inputs** | request | target cif, scaffold cif, design spec | chain missing, placeholders left in, hotspot not found |
| **probe** | spec | 2 designs | spec invalid — cheap, so always run it |
| **generate** | spec, seed, N | `*_generated.cif` per design + `per_example_metrics.csv` (written only at the end) | OOM, GPU contention |
| **parse** | structures | per-design records | chain is not a numberable antibody domain |
| **epitope check** | structures | contacts, epitope fraction | CPU-only, so run it *during* generation for a free early signal |
| **refold** | top-N designs | folded structures + `refold_metrics.csv` (keyed by `query_id`, not `design_id`) | the real filter: expect ~90 % loss |
| **validate** | submission table | verdict per design | the organisers' tool is the authority; never substitute a local reimplementation |
| **rank + report** | everything | bundle | — |

---

## 5. Invariants

- **The validator is the authority.** A local novelty implementation exists only
  to filter cheaply before spending compute. When they disagree, the validator
  is right and the local copy is the bug.
- **Rank within one folding model.** Cut-offs were calibrated per model; scores
  from different models are not comparable.
- **Confidence is a filter, not a score.** Interface confidence happily rates
  non-binding pairs highly — calibrate the threshold on known negatives, then
  use it to discard, not to order.
- **The track is derived, never chosen.** A design that misses its track's
  criteria is rejected, not reclassified.
- **Generation is blind.** Every design looks plausible; the model always emits
  something. Nothing is a result until it has survived refolding.

---

## 6. Failure modes seen in practice

| Symptom | Cause | Fix |
| --- | --- | --- |
| `Failed to compute leaving groups for '12-21'` | CDR length placeholders left in the scaffold | strip non-standard residues |
| 0 epitope contacts on every design | generated structures renumbered from 1 | shift to reference numbering |
| novelty table all null | CDRH3 found by framework anchors that a VHH lacks | number with `immunum` instead |
| refold metrics don't join | `refold_metrics.csv` keys on `query_id`; `design` is a path | join on `query_id` |
| designs rejected despite being far from the reference | the bar is distance from **every** known loop, not just the reference | check against the full reference set |

---

## 7. What the agent should decide, and what it must ask

**Decide alone:** seeds, batch sizes, which folding model, retry on transient
failure, how to spend a fixed GPU budget across generate/refold, filtering
thresholds that have a calibration set.

**Ask a human:** the epitope (a scientific claim about mechanism), the scaffold
format, how to trade diversity against validation depth when the deadline binds,
and anything that publishes or submits — the final list is a commitment of
someone else's lab budget.

---

## 8. Worth building next

- An **epitope-from-complex** helper: given a reference complex, emit hotspot
  atoms automatically instead of copying them by hand.
- **Calibrated confidence thresholds** per folding model from a known-negative
  set, so the filter is evidence-based rather than a guessed cut-off.
- **Developability scoring** (aggregation, solubility) before ranking —
  a design that does not express cannot be measured at all.
- A **cost model**: seconds per design per stage, so the agent can plan a budget
  rather than discover it has run out of time.
