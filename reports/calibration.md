# ipTM calibration against known cognate / shuffled pairs

Source: Smorodina et al., `iptm_confidence_scores.zip`. Diagonal = real cognate pair, off-diagonal = shuffled pairing (non-binder by construction).

Negatives are **computational** shuffled pairings, not measured non-binders, and are an easier negative than an optimised design — so every FPR below is a *lower bound* on what we would see on our own candidates.

## Separation

| predictor | n real | n shuffled | median real | median shuffled | AUROC |
|---|---|---|---|---|---|
| af3 | 91 | 8137 | 0.500 | 0.180 | **0.754** |
| boltz2 | 91 | 8190 | 0.911 | 0.884 | **0.645** |
| chai1 | 89 | 7952 | 0.208 | 0.191 | **0.568** |

## Cutoff at a 5% false-positive rate

| predictor | cutoff | real kept (TPR) | shuffled kept (FPR) | enrichment |
|---|---|---|---|---|
| af3 | ipTM >= 0.560 | 49.5% | 5.2% | 9.45x |
| boltz2 | ipTM >= 0.951 | 17.6% | 5.0% | 3.51x |
| chai1 | ipTM >= 0.393 | 16.9% | 5.0% | 3.37x |

## What a 0.8 cutoff buys, the number people actually use

| predictor | real kept | shuffled kept | enrichment |
|---|---|---|---|
| af3 | 41.8% | 1.5% | 28.80x |
| boltz2 | 90.1% | 81.7% | 1.10x |
| chai1 | 1.1% | 0.0% | 89.35x |

## Conclusion

- **Boltz-2 ipTM must not be used to rank designs.** AUROC 0.645; 82% of known non-binders score >= 0.8, against 90% of real pairs — an enrichment of 1.10x. `scripts/refold.sh` folds with Boltz-2, so its `iptm` column carries almost no evidence of binding.
- AF3 ipTM does separate (AUROC 0.754), but only 42% of real pairs reach 0.8, so a 0.8 cutoff on AF3 would also discard most true binders. AF3 is not available locally.
- Therefore rank on **interface PAE** against the organisers' vendored antibody criterion (`min_ipae <= 2.0` A, `data/docs/06-evaluate.md`) rather than on ipTM, and treat an orthogonal predictor (AF2-IG) as the confirmatory check. Note this benchmark contains ipTM only, so the PAE threshold is vendored, not calibrated here.
