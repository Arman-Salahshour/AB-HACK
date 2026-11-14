# Submission candidates

> **The delivered file is [`odyssey.csv`](../odyssey.csv)** — 30 designs (24 optimisation,
> 6 de novo), 30/30 `abguard`-valid, in the organisers' format with the `track` column.
> It is ranked by the scheme described in `docs/EFFICIENCY.md`: structure as a gate,
> then P(high) and predicted K_D. The 20-design table below is the earlier
> structure-first ranking, kept because the disagreement between the two orderings is
> itself a result.

# Earlier ranking — 20 designs, structure-first

> **What is new here.** The ranking disagreement is published rather than resolved. This
> order weights structural evidence first, but the FACS classifier — trained on 524,346
> measured variants and calibrated against the parent at 0.995 — puts designs 20 and 16 on
> top instead. Both orderings are defensible and they select different molecules. Stating
> that lets whoever presents this choose deliberately, and lets the experimental result
> adjudicate between two recorded predictions rather than one.

All 20 pass every `abguard` filter (format, length, cysteines, IMGT numbering, novelty
against a 610k CDRH3 reference). Ranked on structure reproduction, then epitope
engagement, then interface confidence.

Structures: `designs/submission-structures/` — designed and refolded PDB per design.
Machine-readable: `designs/SUBMISSION-full-metrics.csv` (27 columns, with sequences).

## Binding and structure

| # | CDRH3 | epitope | prec | min_ipae | iptm | pLDDT | RMSD | P(high) | pred K_D |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `WASLYEYYFD` | 20 | 1.00 | 0.300 | 0.971 | 98.7 | 0.99 Å | 0.583 | 68.9 nM |
| 2 | `WASSYYYYFD` | 20 | 1.00 | 0.314 | 0.967 | 98.2 | 0.70 Å | 0.454 | 5.9 nM |
| 3 | `YKSLYEYYFD` | 20 | 1.00 | 0.332 | 0.956 | 97.8 | 0.79 Å | 0.411 | 55.8 nM |
| 4 | `WSFSVDYAFD` | 20 | 1.00 | 0.337 | 0.956 | 98.2 | 0.84 Å | 0.485 | 95.2 nM |
| 5 | `WASHYDYYFD` | 19 | 1.00 | 0.310 | 0.967 | 98.4 | 0.88 Å | 0.113 | 85.6 nM |
| 6 | `ATRAVDYAFD` | 19 | 0.95 | 0.313 | 0.961 | 98.3 | 0.91 Å | 0.108 | 213.1 nM |
| 7 | `YASSYDYYFD` | 19 | 1.00 | 0.328 | 0.962 | 97.3 | 0.78 Å | 0.419 | 51.7 nM |
| 8 | `STRAVDYAFD` | 19 | 0.95 | 0.344 | 0.952 | 98.5 | 0.91 Å | 0.283 | 211.2 nM |
| 9 | `STLSVSYAFD` | 18 | 1.00 | 0.294 | 0.969 | 98.6 | 0.85 Å | 0.034 | 255.1 nM |
| 10 | `WNYAYDYAMD` | 18 | 1.00 | 0.314 | 0.962 | 98.6 | 0.68 Å | 0.063 | 23.2 nM |
| 11 | `YASLYEYYFD` | 18 | 1.00 | 0.314 | 0.965 | 98.7 | 0.89 Å | 0.603 | 65.7 nM |
| 12 | `WTSSVDYAHD` | 18 | 1.00 | 0.316 | 0.960 | 98.5 | 0.98 Å | 0.565 | 28.1 nM |
| 13 | `WHSLYYYYFD` | 18 | 1.00 | 0.318 | 0.964 | 98.0 | 0.78 Å | 0.121 | 13.9 nM |
| 14 | `STSAVDYAFD` | 18 | 1.00 | 0.329 | 0.956 | 98.6 | 0.70 Å | 0.248 | 196.1 nM |
| 15 | `STSSVGWAFD` | 18 | 1.00 | 0.335 | 0.956 | 98.2 | 0.89 Å | 0.689 | 143.8 nM |
| 16 | `WVYAGSYAHD` | 18 | 1.00 | 0.341 | 0.954 | 98.0 | 1.04 Å | 0.759 | 4.9 nM |
| 17 | `STRSVDYAFD` | 17 | 0.90 | 0.292 | 0.969 | 98.5 | 0.70 Å | 0.363 | 206.0 nM |
| 18 | `FPWPAIQAFD` | 17 | 0.94 | 0.304 | 0.964 | 98.1 | 0.73 Å | 0.665 | 121.0 nM |
| 19 | `SRRAVDYAFD` | 17 | 0.94 | 0.306 | 0.966 | 98.6 | 0.86 Å | 0.569 | 202.8 nM |
| 20 | `YASLHDYYFD` | 17 | 1.00 | 0.309 | 0.969 | 98.6 | 0.76 Å | 0.883 | 87.7 nM |
| — | *trastuzumab* | — | — | — | — | — | — | **0.995** | **0.22 nM** |

## Novelty and developability

| # | CDRH3 | muts from Tz | nearest known | liabilities | high | in CDRH3 | pI | charge pH 6 | GRAVY | instability | len |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `WASLYEYYFD` | 7 | +3 | 6 | 2 | 1 | 8.60 | +4.3 | -0.336 | 38.3 | 227 |
| 2 | `WASSYYYYFD` | 7 | +3 | 6 | 2 | 1 | 8.79 | +5.2 | -0.346 | 40.5 | 227 |
| 3 | `YKSLYEYYFD` | 8 | +4 | 5 | 2 | 0 | 8.78 | +5.3 | -0.363 | 36.1 | 227 |
| 4 | `WSFSVDYAFD` | 6 | +3 | 6 | 2 | 1 | 8.61 | +4.2 | -0.314 | 39.1 | 227 |
| 5 | `WASHYDYYFD` | 7 | +3 | 6 | 2 | 1 | 8.60 | +4.7 | -0.367 | 41.6 | 227 |
| 6 | `ATRAVDYAFD` | 7 | +3 | 5 | 2 | 0 | 8.82 | +5.2 | -0.322 | 36.6 | 227 |
| 7 | `YASSYDYYFD` | 8 | +3 | 5 | 2 | 0 | 8.59 | +4.2 | -0.358 | 39.3 | 227 |
| 8 | `STRAVDYAFD` | 7 | +3 | 5 | 2 | 0 | 8.82 | +5.2 | -0.333 | 38.5 | 227 |
| 9 | `STLSVSYAFD` | 7 | +3 | 5 | 2 | 0 | 8.82 | +5.2 | -0.296 | 39.2 | 227 |
| 10 | `WNYAYDYAMD` | 5 | +4 | 7 | 2 | 2 | 8.60 | +4.2 | -0.360 | 41.8 | 227 |
| 11 | `YASLYEYYFD` | 8 | +3 | 5 | 2 | 0 | 8.59 | +4.3 | -0.337 | 37.1 | 227 |
| 12 | `WTSSVDYAHD` | 6 | +3 | 7 | 2 | 2 | 8.61 | +4.7 | -0.356 | 38.4 | 227 |
| 13 | `WHSLYYYYFD` | 7 | +3 | 6 | 2 | 1 | 8.79 | +5.7 | -0.348 | 41.4 | 227 |
| 14 | `STSAVDYAFD` | 7 | +3 | 6 | 2 | 1 | 8.61 | +4.2 | -0.317 | 38.5 | 227 |
| 15 | `STSSVGWAFD` | 8 | +4 | 7 | 2 | 2 | 8.83 | +5.2 | -0.313 | 38.5 | 227 |
| 16 | `WVYAGSYAHD` | 5 | +3 | 6 | 2 | 1 | 8.81 | +5.7 | -0.333 | 39.2 | 227 |
| 17 | `STRSVDYAFD` | 7 | +3 | 5 | 2 | 0 | 8.82 | +5.2 | -0.345 | 40.4 | 227 |
| 18 | `FPWPAIQAFD` | 8 | +4 | 6 | 2 | 1 | 8.83 | +5.2 | -0.306 | 37.7 | 227 |
| 19 | `SRRAVDYAFD` | 7 | +3 | 5 | 2 | 0 | 8.98 | +6.2 | -0.350 | 41.9 | 227 |
| 20 | `YASLHDYYFD` | 8 | +3 | 5 | 2 | 0 | 8.60 | +4.7 | -0.346 | 37.4 | 227 |
| — | *trastuzumab* | 0 | 0 | **8** | **3** | **3** | 8.61 | +4.3 | −0.334 | 40.7 | 227 |

## Sequences

Light chain is trastuzumab's, identical in all 20:

```
DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
```

| # | CDRH3 | heavy chain |
| --- | --- | --- |
| 1 | `WASLYEYYFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWASLYEYYFDYWGQGTLVTVSS` |
| 2 | `WASSYYYYFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWASSYYYYFDYWGQGTLVTVSS` |
| 3 | `YKSLYEYYFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRYKSLYEYYFDYWGQGTLVTVSS` |
| 4 | `WSFSVDYAFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWSFSVDYAFDYWGQGTLVTVSS` |
| 5 | `WASHYDYYFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWASHYDYYFDYWGQGTLVTVSS` |
| 6 | `ATRAVDYAFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRATRAVDYAFDYWGQGTLVTVSS` |
| 7 | `YASSYDYYFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRYASSYDYYFDYWGQGTLVTVSS` |
| 8 | `STRAVDYAFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRSTRAVDYAFDYWGQGTLVTVSS` |
| 9 | `STLSVSYAFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRSTLSVSYAFDYWGQGTLVTVSS` |
| 10 | `WNYAYDYAMD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWNYAYDYAMDYWGQGTLVTVSS` |
| 11 | `YASLYEYYFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRYASLYEYYFDYWGQGTLVTVSS` |
| 12 | `WTSSVDYAHD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWTSSVDYAHDYWGQGTLVTVSS` |
| 13 | `WHSLYYYYFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWHSLYYYYFDYWGQGTLVTVSS` |
| 14 | `STSAVDYAFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRSTSAVDYAFDYWGQGTLVTVSS` |
| 15 | `STSSVGWAFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRSTSSVGWAFDYWGQGTLVTVSS` |
| 16 | `WVYAGSYAHD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWVYAGSYAHDYWGQGTLVTVSS` |
| 17 | `STRSVDYAFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRSTRSVDYAFDYWGQGTLVTVSS` |
| 18 | `FPWPAIQAFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRFPWPAIQAFDYWGQGTLVTVSS` |
| 19 | `SRRAVDYAFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRSRRAVDYAFDYWGQGTLVTVSS` |
| 20 | `YASLHDYYFD` | `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRYASLHDYYFDYWGQGTLVTVSS` |

## How to read this

- **epitope / prec** — HER2 residues contacted within 4.5 Å of the mapped trastuzumab
  epitope, and the fraction of all contacts falling inside it. Precision 1.00 means no
  off-target contact anywhere on the receptor.
- **min_ipae / iptm** — interface confidence from refolding with Boltz-2 given the
  **sequence alone**. Lower ipae and higher iptm are better.
- **RMSD** — backbone agreement between the design and that independent rebuild. Under
  1 Å means a model that never saw the design reconstructed it.
- **P(high)** — FACS classifier over 524,346 measured variants. It scores trastuzumab
  at 0.995, so the column is calibrated; none of our designs approaches that.
- **pred K_D** — SPR regressor over 190 measured affinities. Predicts the parent at
  0.22 nM against a measured 0.029–0.178 nM.
- **liabilities** — LAP reference via `abguard`. Every design carries fewer than
  trastuzumab's 8, and several have none in CDRH3.

## What these numbers do not cover

Thermostability, toxicity, immunogenicity, solubility and aggregation are **not
predicted** — no usable tool was reachable. See `docs/TOOLS.md` for what was tried.
The affinity columns are predictions from models trained on this system, not
measurements, and no design is predicted to beat trastuzumab.

## A ranking caveat

This order weights structural evidence first. The FACS classifier disagrees: design 20
(`YASLHDYYFD`) has the highest P(high) of all twenty at 0.883, and design 16
(`WVYAGSYAHD`) pairs 0.759 with the best predicted K_D at 4.9 nM. If you trust the
measured-data model over structure prediction, those two belong at the top.

