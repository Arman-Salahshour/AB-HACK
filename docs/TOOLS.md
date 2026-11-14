# Tools used

Every tool that produced a number in this project, what it was used for, where it ran,
and what it cannot tell you. Tools we tried and could not use are listed too, because the
gaps matter when reading the results.

---

## Design

| Tool | Version / weights | Used for | Where |
| --- | --- | --- | --- |
| **RFD4-Proteina** | `her2-novel-2026-11-09` | CDRH3 redesign on trastuzumab's framework — every submitted design | local H100, Apheris Foundry Docker |
| **RFD4-Proteina** | `antibody-2026-10-30` | de novo VHH design (track dropped) | local H100 |
| **RFD4-Proteina** | `her2-2026-11-07`, `--sequence-tsr-k 1` | exploratory sampling (not submitted) | local H100 |
| **BoltzGen** | `nanobody-anything` | second, independent de novo attempt | JapanFold |
| **RFD4-Proteina inverse folding** | `--inverse-folding true` | sequence decode on a fixed backbone, benchmarked at 21.1 s/seq | local H100 |

Design inputs came from PDB **1N8Z** (trastuzumab–HER2) and **7MN8** (HER2–HER3–neuregulin
with trastuzumab Fab), the latter supplied by the team.

## Structure validation

| Tool | Used for | Where | Limitation |
| --- | --- | --- | --- |
| **Boltz-2** | refolding designs from **sequence alone** and comparing to the designed pose — the primary filter | local H100 | ran without MSAs; `pass` column undefined for this setup, so thresholds are ours |
| **AF2-IG** | independent scoring of designed binders; killed the de novo track at ipTM 0.11–0.34 | JapanFold | no MSA, no sampling — deterministic |
| **OpenDDE-abag** | antibody–antigen specific second opinion, **with MSAs**; confirmed the optimisation track at ipTM 0.63–0.96 | JapanFold | run unguided on purpose; guided docking cannot validate |

## Filtering and compliance

| Tool | Used for |
| --- | --- |
| **abguard** | the organisers' hard filters: format, length, unpaired cysteines, IMGT numberability, and the novelty gate against a 610k CDRH3 reference. **The authority** — it rejected designs our local filter passed |
| **immunum 2.0.0** | IMGT numbering and CDR extraction, the same library abguard uses |
| **abguard liability-report** | LAP Antibody Liability Reference (PLOS Comput Biol 2024;20(3):e1011881) — deamidation, isomerisation, oxidation, glycosylation, fragmentation |

## Scoring and ranking

| Tool | Trained on | Used for | Validation |
| --- | --- | --- | --- |
| **FACS classifier** (ours) | 524,346 measured trastuzumab CDRH3 variants | P(high binder) | 0.91 held-out accuracy; scores the parent at **0.995** |
| **SPR regressor** (ours) | 190 measured K_D values | predicted K_D | CV R² 0.56; predicts the parent at **0.22 nM** vs measured 0.029–0.178 nM |
| **`abhack.epitope`** (ours) | — | contacts within 4.5 Å against the mapped epitope | geometry only |
| **`scripts/structural_stability.py`** (ours) | — | CDRH3 vs framework pLDDT, buried SASA, interface contacts | descriptive |
| **Biopython ProtParam** | — | pI, charge at pH, GRAVY, instability index, aromaticity | sequence heuristics |
| **gemmi** | — | structure parsing, chain surgery, PDB↔mmCIF conversion | — |

Both trained models pass the obvious control: run on trastuzumab's own CDRH3 they return
0.995 and 0.22 nM. That is why their verdict on the designs is worth reading — and their
verdict is that none approaches the parent.

## Literature and regulatory evidence

| Tool | Used for |
| --- | --- |
| **AMASS BiomedCore** | published determinants of antibody Tm, CDR tryptophan oxidation, deamidation and aggregation |
| **AMASS RegulatoryCore** | FDA/EMA labels: formulation pH 6, measured anti-drug antibody rate (1/903), cardiomyopathy incidence (CHF 0.4–3.2 %) |
| **Adaptyv EGFR competition reports** | hit rates by design approach — optimised 18–28 %, best de novo 9 % — which justified dropping the de novo track |

---

## Tried and could not use

Listing these because their absence shapes what the benchmark can claim.

| Tool | Intended for | Why not |
| --- | --- | --- |
| **ProteinMPNN (soluble weights)** | solubility-optimised sequences, recommended by the organisers | cloned with weights, torch installed and CUDA working — **sandbox blocks executing downloaded code**, twice |
| **Aggrescan3D** | aggregation propensity | Python 2 package; build fails on Python 3.12 |
| **NetSolP** | solubility | not on PyPI; DTU download plus ESM weights |
| **SoluProt** | solubility | not on PyPI; web service |
| **FoldX / ThermoMPNN / DeepSTABp** | folding ΔΔG and Tm | licensed or blocked from execution |
| **NetMHCIIpan** | immunogenicity (T-cell epitopes) | licensed |
| **Nebius `stab-ddg`** | binding ΔΔG — a third affinity estimate | endpoint located (`POST /models/stab-ddg:submit`, returns 422 not 404) but the batch protocol needs an artifact upload whose tooling (`~/japanfold/agent/forge.py`) is missing from the provided folder |

**Consequences, stated plainly:**

- **Thermostability is not predicted.** It is a judged axis. What we report are structural
  correlates (VH–VL interface contacts, CDRH3 glycine content) drawn from the literature,
  not a temperature or a free energy.
- **Toxicity is not predicted.** Nebius's `ctoxpred2` and `admet-ai` are small-molecule
  tools. We report the parent's *measured* safety from its FDA label, plus absence-of-red-flag
  checks (no glycosylation sequons, no RGD, no free cysteine, CDR charge matching the parent).
- **Immunogenicity is not predicted**, deliberately. Computational immunogenicity tools in
  this class have limited validated predictive power
  ([FASEB J 2008;22:563](https://faseb.onlinelibrary.wiley.com/doi/10.1096/fasebj.22.2_supplement.563)),
  so the measured 1/903 anti-drug antibody rate for the parent is reported instead.
- **Solubility and aggregation are not predicted**, since none of the four recommended
  tools could be run.

---

## Compute

| | |
| --- | --- |
| Local | 1× NVIDIA H100 80 GB, ~5 h total |
| JapanFold | **$0.20** of $1,100 (AF2-IG, BoltzGen, OpenDDE-abag) |
| Nebius | model catalogue queried; no jobs run |
| AMASS | literature and regulatory search |
