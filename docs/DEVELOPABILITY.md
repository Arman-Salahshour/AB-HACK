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

## 4. Toxicity

**Not predictable from sequence, and not claimed here.** The strongest available argument
is structural: every design is trastuzumab's framework with ten CDRH3 positions changed,
and the light chain is trastuzumab's unchanged. The scaffold is that of a drug approved
by FDA and EMA since 1998 and marketed as at least eight products.

Two caveats that follow from the same fact:

- **Trastuzumab's known class toxicity is cardiac** (left-ventricular dysfunction),
  target-mediated through HER2 in cardiomyocytes. Designs binding the same domain IV
  epitope would be expected to share it.
- **A changed CDRH3 is a new antigenic determinant.** Immunogenicity risk cannot be
  inferred from the parent and needs its own assessment.

---

## 5. Recommendation

**If ranking on developability alone**, the picks change:

| Rank | Design | CDRH3 | Case |
| --- | --- | --- | --- |
| 1 | **D06** | ATRAVDYAFD | no CDRH3 liability, best instability index (36.6), fewest total liabilities |
| 2 | **D09** | STLSVSYAFD | no CDRH3 liability, least hydrophobic (GRAVY −0.296) |
| 3 | **D04** | YASSYDYYFD | no CDRH3 liability, 4 Trp vs trastuzumab's 5 |

**But these are not the structurally best designs.** D03 (`WASSYYYYFD`) leads on
interface and predicted affinity (20/22 epitope contacts, precision 1.00, predicted
K_D 5.9 nM) and carries a CDRH3 Trp oxidation site.

That is a genuine trade-off, not a ranking error: the binding evidence favours D03, the
developability evidence favours D06/D09. With 20 submission slots, submit both groups —
and note that the Trp risk is mitigable through distal framework charge changes without
touching CDRH3 or the novelty gate.

---

## Sources

- Hageman et al., *Pharm Res* 2018 — [Trp oxidation in CDRs disrupts antigen binding](https://pubmed.ncbi.nlm.nih.gov/30536043/)
- Alam et al., *J Pharm Sci* 2019 — [Met/Trp oxidation and deamidation on stability and aggregation](https://pubmed.ncbi.nlm.nih.gov/31678251/)
- Alam et al., *Mol Pharm* 2019 — [Deamidation and colloidal stability, pH-dependent](https://pubmed.ncbi.nlm.nih.gov/30916563/)
- Pavon et al., *Anal Chem* 2019 — [Selective Trp oxidation of mAbs](https://pubmed.ncbi.nlm.nih.gov/30608647/)
- Lenka et al., *mAbs* 2026 — [Local electrostatics governs Trp oxidation; distal mutations reduce it](https://pubmed.ncbi.nlm.nih.gov/42590841/)
- FDA/EMA labels via AMASS RegulatoryCore: Ontruzant, Kanjinti, Trazimera, Herzuma, Kadcyla
