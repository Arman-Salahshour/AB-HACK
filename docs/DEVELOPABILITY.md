# Developability benchmark — seed 21 candidates

Twelve top designs from the epitope-steered run, against trastuzumab as the reference.
Computed metrics from `abguard` (LAP liability reference) and Biopython ProtParam;
interpretation grounded in the literature and FDA/EMA labels via AMASS.

**What is and is not measured here.** Sequence liabilities, isoelectric point, charge,
hydrophobicity and the instability index are computed. **Toxicity, immunogenicity,
aggregation rate and thermostability are not** — those need assays or predictors
(Aggrescan3D, NetSolP, TAP) that are not installed. AMASS searches literature and
regulatory records; it does not compute properties from a sequence.

---

## 1. Sequence liabilities

`abguard team-validate --all-liabilities`, LAP Antibody Liability Reference
(PLOS Comput Biol 2024;20(3):e1011881).

| Design | CDRH3 | total | high | in CDRH3 | CDRH3 motifs |
| --- | --- | --- | --- | --- | --- |
| **D04** | YASSYDYYFD | 5 | 2 | **0** | – |
| **D06** | ATRAVDYAFD | 5 | 2 | **0** | – |
| **D09** | STLSVSYAFD | 5 | 2 | **0** | – |
| D01 | WSFSVDYAFD | 6 | 2 | 1 | TrpOx:W |
| D02 | WASLYEYYFD | 6 | 2 | 1 | TrpOx:W |
| D03 | WASSYYYYFD | 6 | 2 | 1 | TrpOx:W |
| D05 | WASHYDYYFD | 6 | 2 | 1 | TrpOx:W |
| D07 | STSAVDYAFD | 6 | 2 | 1 | FragM:TS |
| D10 | WHSLYYYYFD | 6 | 2 | 1 | TrpOx:W |
| D08 | WNYAYDYAMD | 7 | 2 | 2 | MetOx:M, TrpOx:W |
| D11 | WTSSVDYAHD | 7 | 2 | 2 | FragM:TS, TrpOx:W |
| D12 | STSSVGWAFD | 7 | 2 | 2 | FragM:TS, TrpOx:W |
| *trastuzumab* | WGGDGFYAMD | **8** | **3** | 7 in CDRs | Isom ×2, DeAmd ×3, TrpOx, MetOx |

**Every design carries fewer liabilities than trastuzumab** (5–7 vs 8 total, 2 vs 3 high).
The two remaining high-severity flags are in the framework and CDR1/CDR2, inherited
unchanged from trastuzumab — so they are present in an approved drug.

**Three designs have no liability in CDRH3 at all**: D04, D06, D09.

---

## 2. Why CDRH3 tryptophan matters

Nine of twelve designs begin CDRH3 with tryptophan, flagged `TrpOx`. This is not
cosmetic:

- **Trp oxidation within or near a CDR increases conformational flexibility and
  disrupts antigen binding** (Hageman et al., *Pharm Res* 2018, [PMID 30536043](https://pubmed.ncbi.nlm.nih.gov/30536043/)).
- **Trp oxidation and Asn deamidation promote aggregation across a wide range of
  solution conditions**, by a different mechanism from Met oxidation, which instead
  reduces conformational stability (Alam et al., *J Pharm Sci* 2019, [PMID 31678251](https://pubmed.ncbi.nlm.nih.gov/31678251/)).
- A CDR-H3 Trp is a recognised critical quality attribute requiring analytical control
  during manufacture (Pavon et al., *Anal Chem* 2019, [PMID 30608647](https://pubmed.ncbi.nlm.nih.gov/30608647/)).

**Trastuzumab itself carries a CDR Trp oxidation site**, so this is not disqualifying —
but it is a real risk our three Trp-free designs avoid.

**A mitigation exists that does not touch the Trp.** Oxidation susceptibility is governed
by local electrostatics as well as solvent accessibility; distal charge-altering
mutations reduced oxidation ~50 % in four of five re-engineered antibodies with binding
preserved in two (Lenka et al., *mAbs* 2026, [PMID 42590841](https://pubmed.ncbi.nlm.nih.gov/42590841/)).
Since the hackathon's novelty gate inspects only CDRH3, such framework changes are free
against the filter.

---

## 3. Biophysical panel

| Design | pI | charge pH 7.4 | charge pH 5.5 | GRAVY | instability | Cys | Met | Trp |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| *trastuzumab* | 8.61 | +2.6 | +5.4 | −0.334 | 40.7 | 4 | 3 | 5 |
| D01 WSFSVDYAFD | 8.61 | +2.6 | +5.4 | −0.314 | 39.1 | 4 | 2 | 5 |
| D02 WASLYEYYFD | 8.60 | +2.6 | +5.4 | −0.336 | 38.3 | 4 | 2 | 5 |
| D03 WASSYYYYFD | 8.79 | +3.6 | +6.3 | −0.346 | 40.5 | 4 | 2 | 5 |
| D04 YASSYDYYFD | 8.59 | +2.6 | +5.4 | −0.358 | 39.3 | 4 | 2 | 4 |
| D06 ATRAVDYAFD | 8.82 | +3.6 | +6.4 | −0.322 | **36.6** | 4 | 2 | 4 |
| D09 STLSVSYAFD | 8.82 | +3.6 | +6.3 | **−0.296** | 39.2 | 4 | 2 | 4 |
| D10 WHSLYYYYFD | 8.79 | +3.6 | **+7.1** | −0.348 | 41.4 | 4 | 2 | 5 |

(Full table: `designs/s21-biophysics.csv`.)

**Formulation context.** Every approved trastuzumab product is formulated at **pH 6** in
L-histidine with polysorbate 20 (Ontruzant, Kanjinti, Trazimera, Herzuma FDA labels);
Kadcyla at pH 5.0. At pH 6 our designs carry roughly +5 to +7 net charge — comfortably
distant from their pI of ~8.6–8.8, which favours colloidal stability. No design sits
near its pI at formulation pH.

**Deamidation is pH-dependent**: it compromises colloidal stability and raises
aggregation most at low pH (3.8), which is reached during purification rather than
formulation (Alam et al., *Mol Pharm* 2019, [PMID 30916563](https://pubmed.ncbi.nlm.nih.gov/30916563/)).
Relevant to process design, not to the final product.

**Instability index** (Guruprasad 1990; < 40 predicts stability in vitro): designs span
36.6–41.8 against trastuzumab's 40.7. Five score below trastuzumab, D06 best at 36.6.
The index is a crude sequence heuristic and should not be read as a thermostability
prediction — which is one of the three judged axes and remains unmeasured.

**Cysteines**: 4 in every design, all paired — two canonical intradomain disulfides, no
free thiol. This is a hard filter and all designs pass.

---

## 4. Predicted binding

From the two models trained on this exact system. Full table:
`designs/s21-full-panel.csv`.

| Design | CDRH3 | P(high binder) | pred K_D | epitope | min_ipae | RMSD |
| --- | --- | --- | --- | --- | --- | --- |
| ***trastuzumab*** | WGGDGFYAMD | **0.995** | **0.22 nM** | – | – | – |
| D12 | STSSVGWAFD | 0.689 | 143.8 nM | 18 | 0.335 | 0.89 Å |
| D02 | WASLYEYYFD | 0.583 | 68.9 nM | **20** | 0.301 | 0.99 Å |
| D11 | WTSSVDYAHD | 0.565 | 28.1 nM | 18 | 0.316 | 0.98 Å |
| D01 | WSFSVDYAFD | 0.485 | 95.2 nM | **20** | 0.337 | 0.84 Å |
| **D03** | **WASSYYYYFD** | 0.454 | **5.9 nM** | **20** | 0.314 | **0.70 Å** |
| D04 | YASSYDYYFD | 0.419 | 51.7 nM | 19 | 0.328 | 0.78 Å |
| D07 | STSAVDYAFD | 0.248 | 196.1 nM | 18 | 0.329 | 0.70 Å |
| D10 | WHSLYYYYFD | 0.121 | 13.9 nM | 18 | 0.318 | 0.78 Å |
| D05 | WASHYDYYFD | 0.113 | 85.6 nM | 19 | 0.310 | 0.88 Å |
| D06 | ATRAVDYAFD | 0.108 | 213.1 nM | 19 | 0.313 | 0.91 Å |
| D08 | WNYAYDYAMD | 0.063 | 23.2 nM | 18 | 0.314 | 0.68 Å |
| D09 | STLSVSYAFD | 0.034 | 255.1 nM | 18 | 0.295 | 0.85 Å |

**The parent is the positive control, and the models pass it.** Run on trastuzumab's own
CDRH3 they return P(high) = 0.995 and K_D 0.22 nM against a measured 0.029–0.178 nM —
accurate on the one sequence whose answer is known. That makes the rest of the column
worth reading, and the reading is sobering: **no design comes close.** The best predicted
affinity is 27× weaker than the parent, and every P(high) is below 0.70 against 0.995.

Note also that the models and the structure disagree. `min_ipae` and RMSD are tightly
clustered across all twelve (0.29–0.34, 0.68–0.99 Å) while P(high) spans 0.03 to 0.69.
Structure prediction cannot separate these designs; the sequence models can, and they
were trained on measured binding. Where they conflict, prefer the measured-data model.

### A data-integrity note

Ten of the refolded designs share a CDRH3 with another sample. The same loop can come
back with different poses — `WASSYYYYFD` appears as both `_30` (14 epitope contacts,
1.16 Å) and `_94` (20 contacts, 0.70 Å). The table above keeps the best-posed sample per
loop. An earlier version of this panel joined on CDRH3 and silently picked an arbitrary
one, which understated D03.

---

## 5. Toxicity

**Not predictable from sequence, and not claimed here.** What can be checked, was:

| Flag | Result |
| --- | --- |
| N-linked glycosylation sequons (N-X-S/T) | **0** in all designs and in trastuzumab |
| RGD integrin-binding motif | **absent** in all |
| Free (unpaired) cysteine | **none**; 4 Cys, 2 disulfides, in all |
| CDR positive charge (polyreactivity / fast-clearance proxy) | **2–3**, same as trastuzumab's 2 |
| CDR net charge | −1 to −2, same as trastuzumab's −2 |

Low positive charge across the CDRs is the useful one: high positive CDR charge is
associated with polyreactivity and rapid clearance, and these designs sit exactly where
the parent does.

**The structural argument.** Every design is trastuzumab's framework with ten CDRH3
positions changed and the light chain unchanged — the scaffold of a drug approved by FDA
and EMA since 1998 and marketed as at least eight products.

Two caveats follow from that same fact:

- **Trastuzumab's known class toxicity is cardiac** (left-ventricular dysfunction),
  target-mediated through HER2 in cardiomyocytes. Designs binding the same domain IV
  epitope would be expected to share it.
- **A changed CDRH3 is a new antigenic determinant.** Immunogenicity cannot be inferred
  from the parent. Assessing it needs MHC-II binding prediction (NetMHCIIpan) or an
  ADA assay — neither available here.

---

## 6. Head to head: seed 11 vs seed 21 best

Every test applied to the same three molecules, with trastuzumab as reference.
Source: `designs/head-to-head.csv`.

| Test | trastuzumab | D53 (seed 11) | D03 (seed 21) | D04 (seed 21) |
| --- | --- | --- | --- | --- |
| CDRH3 | WGGDGFYAMD | WTSSVGYAHD | WASSYYYYFD | YASSYDYYFD |
| **Structure** | | | | |
| `min_ipae` | – | 0.475 | **0.314** | 0.328 |
| `iptm` | – | 0.777 | **0.967** | 0.962 |
| backbone RMSD | – | 0.824 Å | **0.704 Å** | 0.781 Å |
| binder pLDDT | – | 97.6 | **98.2** | 97.9 |
| **Epitope** | | | | |
| contacts | – | 15/17 | **20** | 19 |
| precision | – | 0.88 | **1.00** | **1.00** |
| **Binding (measured-data models)** | | | | |
| P(high binder) | **0.995** | **0.780** | 0.454 | 0.419 |
| predicted K_D | **0.22 nM** | 22.2 nM | **5.9 nM** | 51.7 nM |
| **Developability** | | | | |
| liabilities, total | 8 | 7 | 6 | **5** |
| high severity | 3 | 2 | 2 | **2** |
| in CDRH3 | 3 | 2 | 1 | **0** |
| instability index | 40.7 | **38.3** | 40.5 | 39.3 |
| pI | 8.61 | 8.82 | 8.79 | 8.59 |
| charge at pH 5.5 | +5.4 | +7.1 | +6.3 | +5.4 |
| GRAVY | −0.334 | −0.342 | −0.346 | **−0.358** |
| CDR positive charge | 2 | 2 | 2 | 2 |
| glycosylation sequons | 0 | 0 | 0 | 0 |

### Reading it

**Steering won on structure and epitope, decisively.** D03 beats D53 on every structural
measure — `iptm` +24 %, RMSD −15 %, `min_ipae` −34 % — and reaches perfect epitope
precision. That improvement came entirely from replacing 4 hotspot atoms with the
17-residue measured epitope.

**D53 still wins on the binding model that matters most.** P(high) 0.780 against D03's
0.454 — and that model is trained on 524,000 measured variants, with the parent as a
passing positive control (0.995). The two models disagree with each other too: D03 has
the better predicted K_D (5.9 vs 22.2 nM) but the worse P(high).

**Developability improves monotonically across the three**: 8 → 7 → 6 → 5 liabilities,
3 → 2 → 1 → 0 in CDRH3. D04 is the cleanest molecule we have and the only one with no
CDRH3 liability at all, at the cost of ~9× predicted affinity against D03.

**Nothing approaches trastuzumab on affinity.** The parent is 27× tighter than our best
predicted design and has P(high) 0.995. Every design is also *cleaner* than the parent on
liabilities — the parent carries 3 high-severity flags and 3 in CDRH3 — so the trade we
are making is affinity for developability, in the wrong direction for S-tier.

### Which to submit

No single design dominates. The three are different bets:

- **D03** — best structure and epitope, best predicted affinity, one Trp flag
- **D53** — best binding probability, independent OpenDDE-abag confirmation at ipTM 0.951
- **D04** — cleanest developability, no CDRH3 liability, weakest affinity

Submit all three near the top. The disagreement between them is information, not noise:
if the measured-data model is right, D53 binds; if structure is right, D03 does.

## 7. Recommendation

**If ranking on developability alone**, the picks change:

| Rank | Design | CDRH3 | Developability case | but P(high) |
| --- | --- | --- | --- | --- |
| 1 | **D06** | ATRAVDYAFD | no CDRH3 liability, best instability index (36.6) | 0.108 |
| 2 | **D09** | STLSVSYAFD | no CDRH3 liability, least hydrophobic (−0.296) | **0.034** |
| 3 | **D04** | YASSYDYYFD | no CDRH3 liability, 4 Trp vs trastuzumab's 5 | 0.419 |

**The three cleanest designs are also among the weakest predicted binders.** D09 has no
liability in CDRH3 and the lowest P(high) of all twelve, 0.034. A clean molecule that
does not bind is worth nothing, so developability cannot lead the ranking.

**The best overall compromise is D03 (`WASSYYYYFD`)**: 20 epitope contacts at precision
1.00, tightest structure agreement (0.70 Å), best predicted affinity by 2.4× (5.9 nM),
mid-range P(high) 0.454 — carrying one CDRH3 Trp oxidation flag, as trastuzumab does.

**D04 (`YASSYDYYFD`) is the one design that is strong on both**: no CDRH3 liability at
all, 19 epitope contacts, 0.78 Å, P(high) 0.419, 51.7 nM. It gives up ~9× predicted
affinity against D03 to remove the Trp risk entirely.

Suggested split of the 20 slots: lead with D03 and D04, carry D02/D01/D11 on binding
evidence, and include D06/D09 as developability-first entries. The Trp risk in the
first group is mitigable through distal framework charge changes, which the novelty
gate does not inspect.

---

## Sources

- Hageman et al., *Pharm Res* 2018 — [Trp oxidation in CDRs disrupts antigen binding](https://pubmed.ncbi.nlm.nih.gov/30536043/)
- Alam et al., *J Pharm Sci* 2019 — [Met/Trp oxidation and deamidation on stability and aggregation](https://pubmed.ncbi.nlm.nih.gov/31678251/)
- Alam et al., *Mol Pharm* 2019 — [Deamidation and colloidal stability, pH-dependent](https://pubmed.ncbi.nlm.nih.gov/30916563/)
- Pavon et al., *Anal Chem* 2019 — [Selective Trp oxidation of mAbs](https://pubmed.ncbi.nlm.nih.gov/30608647/)
- Lenka et al., *mAbs* 2026 — [Local electrostatics governs Trp oxidation; distal mutations reduce it](https://pubmed.ncbi.nlm.nih.gov/42590841/)
- FDA/EMA labels via AMASS RegulatoryCore: Ontruzant, Kanjinti, Trazimera, Herzuma, Kadcyla
