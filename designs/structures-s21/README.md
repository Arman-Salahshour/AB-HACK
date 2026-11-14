# HER2 candidate structures

Trastuzumab variants with a redesigned CDRH3. Every design passes all
`abguard` filters. Files are named by submission rank.

## Two files per candidate

| Suffix | What it is |
| --- | --- |
| `_designed.pdb` | what RFD4-Proteina drew — the designed pose on HER2 |
| `_boltz2.pdb` | what Boltz-2 rebuilt **from the sequence alone** |

Backbone RMSD between the two is the evidence that matters: under 1 Å means an
independent model reproduced the design. A high RMSD with a good `min_ipae`
means Boltz-2 is confident about a *different* interface than the one designed.

## Chains

| Chain | Contents |
| --- | --- |
| A | trastuzumab VL (unchanged) |
| B | VH, CDRH3 redesigned |
| C | HER2 domain IV |
| D | HER2 domain IV (second segment) |

## Candidates

| Rank | Design | CDRH3 | min_ipae | iptm | pLDDT | RMSD (Å) | P(high) | pred K_D (nM) | muts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | Depitope_72 | `STRSVDYAFD` | 0.292 | 0.969 | 98.5 | 0.70 | – | – | 7 |
| 2 | Depitope_11 | `STLSVSYAFD` | 0.294 | 0.969 | 98.6 | 0.85 | – | – | 7 |
| 3 | Depitope_37 | `WASLYDYYFD` | 0.300 | 0.970 | 98.6 | 0.87 | – | – | 7 |
| 4 | Depitope_76 | `WASLYDYYFD` | 0.300 | 0.970 | 98.6 | 0.84 | – | – | 7 |
| 5 | Depitope_78 | `WASLYEYYFD` | 0.300 | 0.971 | 98.7 | 0.80 | – | – | 7 |
| 6 | Depitope_85 | `WASLYEYYFD` | 0.300 | 0.971 | 98.7 | 0.99 | – | – | 7 |
| 7 | Depitope_88 | `FPWPAIQAFD` | 0.304 | 0.964 | 98.1 | 0.73 | – | – | 8 |
| 8 | Depitope_95 | `SRRAVDYAFD` | 0.306 | 0.966 | 98.6 | 0.86 | – | – | 7 |
| 9 | Depitope_14 | `YASLHDYYFD` | 0.309 | 0.969 | 98.6 | 0.76 | – | – | 8 |
| 10 | Depitope_28 | `WASHYDYYFD` | 0.310 | 0.967 | 98.4 | 0.88 | – | – | 7 |
| 11 | Depitope_73 | `ATRAVDYAFD` | 0.313 | 0.961 | 98.3 | 0.91 | – | – | 7 |
| 12 | Depitope_30 | `WASSYYYYFD` | 0.314 | 0.967 | 98.2 | 1.16 | – | – | 7 |

## Reading the numbers

- `min_ipae` — interface predicted aligned error, lower is better; under 0.5 is good.
- `iptm` — interface confidence, higher is better.
- `binder_plddt` — how well-formed the antibody itself is.
- `RMSD` — backbone agreement between designed and refolded structure.
- `P(high)` — FACS classifier (524k measured variants) on whether the loop binds.
- `pred K_D` — SPR regressor (190 measured affinities). Trastuzumab measures
  0.029–0.178 nM, and no design is predicted to beat it.

Compare scores only within one folding model; cut-offs are calibrated per model.
