# Making the pipeline more efficient

Measured on the hackathon VM (1× H100), 14 November 2026. Every number below comes from a
run in `runs/`, not from documentation.

---

## 1. The main inefficiency is not speed

**The filters run in the wrong order.** We spend 14 s per design on a measurement that does
not discriminate, then rank with a model that costs milliseconds and does.

Across the final 20 designs:

| Signal | Cost per design | Range across the 20 | Discriminates? |
| --- | --- | --- | --- |
| `min_ipae` (Boltz-2 refold) | **14 s** | 0.292 – 0.344 | barely |
| `iptm` (Boltz-2 refold) | (same run) | 0.952 – 0.971 | barely |
| backbone RMSD | (same run) | 0.68 – 1.04 Å | barely |
| **P(high binder)** (FACS classifier) | **~1 ms** | **0.034 – 0.883** | **26× spread** |
| predicted K_D (SPR regressor) | ~1 ms | 4.9 – 255 nM | 52× spread |

And the refold's own pass rate says the same thing. Seed 21 returned **30 strong refolds out
of 30**. A filter that passes everything carries no information: ~9 minutes of GPU bought
confirmation that the batch was healthy, and nothing at all about which designs to rank.

**Reorder:**

```
now       generate ──► refold top-30 by self-score ──► rank by structure + models
                       └─ 9 min GPU, ~0 bits of ranking information

better    generate ──► abguard ──► FACS/SPR models ──► refold top-8 by P(high)
                       └─ ms ─────┴─ ms ─────────────┴─ 2 min, as a sanity check
```

The refold still belongs in the pipeline — it is what killed the de novo track, where the
pass rate was 0/40 and the information content was total. But **its value collapses as the
pass rate approaches 100 %**, and on the optimisation track it has.

**Rule of thumb:** spend GPU on a filter in proportion to how much it rejects. Track its
pass rate and cut the batch size when it stops rejecting.

---

## 2. Cost model, per 100 designs

Measured:

| Stage | Cost | Yield |
| --- | --- | --- |
| generate | 100 × 18.0 s = **30 min** (9 min with `--compile`) | 100 designs, ~91 unique loops |
| abguard + novelty + epitope + models | **< 1 min** CPU | ~66 eligible |
| refold top 30 | ~14 s each + 2 min overhead = **9 min** | 30 validated |
| **total** | **~40 min** | 66 eligible, 30 validated |

Reallocating the refold budget into generation gives **+30 designs per batch** for the same
wall-clock — a ~30 % larger eligible pool, ranked by the better discriminator.

---

## 3. Flags we left on the table

`runs/opt-epitope-s21/design/config.json` shows what we actually ran:

```
compile: False          ← docs: "faster from about 8 designs per length"; we ran 100
batch_size: 1           ← never raised
recycle_strategy: high/front_loaded
```

- **`--compile true`** was never enabled on any run, and it is worth **3.3×**:

  | setting | per design | 100 designs | 400 designs |
  | --- | --- | --- | --- |
  | `compile: False` (what we ran) | **18.0 s** | 30 min | 2 h |
  | `compile: True` (measured, `runs/compile-test`) | **5.4 s** | **9 min** | **36 min** |

  Measured over 17 samples, excluding the one-time compilation at startup, which
  amortises away on a 100-design batch. We generated 482 designs today without it:
  ~2.5 h of GPU where ~45 min would have done, and that difference is exactly the
  time refolding was short of. One flag, documented as "faster from about 8 designs
  per length", and every batch ran 100 per length.
- **`--batch-size`** above 1 batches inputs; untested here.
- **`--recycle-strategy`** controls compute per design. `high/front_loaded` is the
  antibody recommendation and we kept it — lowering it trades quality for speed and was
  not explored.

---

## 4. Delayering the model does not help

Tested directly. RFD4-Proteina's `--inverse-folding true` freezes the structure and runs
only the sequence head:

| Mode | Per sequence | Output diversity |
| --- | --- | --- |
| full design (diffusion + decode) | **18.0 s** | 91 unique loops per 100 |
| inverse folding (decode only) | **21.1 s** | 12 unique loops per 20 |

**Decode-only is slower**, because the model runs its full trunk either way — there is no
separable lightweight decoder to isolate. Freezing the structure removes a computation the
architecture performs regardless.

Inverse folding is still useful for **refining a validated backbone**: all 20 sequences it
produced sat 6 mutations from trastuzumab, varying only at loop positions 2–4. That is a
cheap way to get related variants of a design you trust, not a way to generate faster.

**Millisecond sequence design needs a different architecture** — ProteinMPNN is 1.7M
parameters against RFD4's full trunk, roughly 20 ms per sequence. It was cloned with its
soluble weights and torch installed, but the sandbox blocks executing downloaded code.
With it, generation stops being the bottleneck entirely and the ordering in section 1
becomes the only thing that matters.

---

## 5. Scheduling

What we did right, and worth keeping:

- **CPU validation runs during GPU generation.** `pipeline validate` is CPU-only, so
  collecting, filtering and scoring a finished batch costs nothing against the next one.
- **Generation and refolding share the H100.** 3.6 GB of 80 GB was in use during
  generation, so running Boltz-2 alongside cost both jobs a little and saved a lot of
  wall-clock.
- **Off-box validation is free parallelism.** AF2-IG, OpenDDE-abag and BoltzGen ran on
  JapanFold for **$0.20 total** while the H100 kept generating. That is the cheapest
  capacity in the whole pipeline and we barely used it.

What we did wrong:

- **Generated before validating.** Four batches of 100 were queued before the first refold
  returned. When seed 1 came back 0/20, three more batches of a doomed configuration were
  already running or queued.
- **Probe first, always.** A 2-design probe costs 90 s and catches every spec error we hit.
  We adopted this after the first failure, not before.

---

## 6. Concrete changes, in order of value

1. **Rank with the measured-data models, refold to confirm.** Refold 5–10 per batch rather
   than 30 while the pass rate stays high. Frees ~7 min of GPU per 100 designs.
2. **Enable `--compile true`.** Measured 3.3× (18.0 → 5.4 s/design). Biggest single win.
3. **Track the refold pass rate as a control signal.** Above ~80 %, shrink the refold set;
   below ~20 %, something is wrong with the configuration and generating more is waste.
4. **Validate the first batch before queueing the next three.** One refold of 20 designs
   (~5 min) decides whether the configuration is worth 90 more minutes of generation.
5. **Push validation off-box.** JapanFold costs cents and runs in parallel with the GPU.
6. **Use inverse folding only on validated backbones**, for variants — not for throughput.
