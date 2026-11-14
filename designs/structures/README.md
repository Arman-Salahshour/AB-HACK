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
| 1 | D53 | `WTSSVGYAHD` | 0.475 | 0.777 | 97.6 | 0.82 | 0.78 | 22.2 | 6 |
| 2 | D50 | `WASSYYYYFD` | 0.434 | 0.794 | 97.2 | 0.89 | 0.45 | 5.9 | 7 |
| 3 | D55 | `YASLYDYYFD` | 0.493 | 0.792 | 97.4 | 0.93 | 0.54 | 87.7 | 8 |
| 4 | D4 | `WHSLYYYYFD` | 0.388 | 0.837 | 95.9 | 0.81 | 0.12 | 13.9 | 7 |
| 5 | D20 | `STSAYDYAFD` | 0.376 | 0.794 | 96.4 | 0.82 | 0.17 | 180.9 | 7 |
| 6 | D64 | `YYSTSMYYFD` | 0.330 | 0.883 | 94.6 | 8.82 | 0.66 | 41.2 | 8 |
| 7 | D10 | `DLSLADRAFD` | 0.378 | 0.868 | 96.5 | 14.17 | 0.66 | 210.5 | 8 |
| 8 | D23 | `DLLPNHAVFD` | 0.414 | 0.805 | 96.6 | 4.33 | 0.05 | 1097.2 | 9 |
| 9 | D0 | `SLSIYDYSFD` | 0.423 | 0.749 | 97.4 | 11.74 | 0.39 | 227.0 | 8 |
| 10 | D83 | `DMRDYPGYFD` | 0.424 | 0.757 | 96.8 | 20.27 | 0.13 | 343.5 | 8 |
| 11 | D41 | `SLLLNHDVFD` | 0.437 | 0.790 | 96.9 | 3.60 | 0.04 | 501.7 | 9 |
| 12 | D81 | `DYSTYSYAFD` | 0.444 | 0.724 | 97.3 | 8.65 | 0.13 | 110.2 | 7 |
| 13 | D80 | `DLTLSTSDFD` | 0.458 | 0.772 | 96.7 | 6.99 | 0.11 | 307.9 | 9 |
| 14 | D9 | `DRSLYYYEFD` | 0.463 | 0.757 | 97.5 | 23.19 | 0.16 | 140.0 | 8 |
| 15 | D21 | `YVPFSVSAFD` | 0.482 | 0.696 | 96.6 | 20.54 | 0.76 | 278.0 | 8 |
| 16 | D68 | `WKSLSEYYFD` | 0.523 | 0.747 | 96.8 | 0.97 | 0.29 | 29.3 | 7 |
| 17 | D77 | `SLSLSSYAYD` | 0.530 | 0.671 | 97.2 | 7.82 | 0.14 | 52.5 | 7 |

## Reading the numbers

- `min_ipae` — interface predicted aligned error, lower is better; under 0.5 is good.
- `iptm` — interface confidence, higher is better.
- `binder_plddt` — how well-formed the antibody itself is.
- `RMSD` — backbone agreement between designed and refolded structure.
- `P(high)` — FACS classifier (524k measured variants) on whether the loop binds.
- `pred K_D` — SPR regressor (190 measured affinities). Trastuzumab measures
  0.029–0.178 nM, and no design is predicted to beat it.

Compare scores only within one folding model; cut-offs are calibrated per model.
