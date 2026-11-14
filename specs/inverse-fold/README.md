# Sequence decode only (inverse folding)

`--inverse-folding true` keeps every backbone coordinate and runs only RFD4-Proteina's
sequence head. No diffusion over structure: the fold you give it is the fold you get
back, with new sequence written onto it.

```bash
bash scripts/inverse_fold.sh                              # the HER2 CDRH3 example
SPEC=specs/mine/design.json N=50 bash scripts/inverse_fold.sh
STRUCT=runs/.../x_generated.cif N=20 bash scripts/inverse_fold.sh   # a structure, no spec
```

## You control what is decoded

Through `C_SEQ`, and nothing else:

| In the spec | Effect |
| --- | --- |
| listed in `C_SEQ` | sequence is **fixed**, copied from the input structure |
| **not** listed | sequence is **decoded** for that backbone |

In the shipped spec, chain B residues 99–108 — trastuzumab's CDRH3 — are absent from
`C_SEQ`, so only those ten positions get new sequence:

```json
"C_SEQ": { "true": [
  { "select": "sel('A')" },                      // light chain: fixed
  { "select": "sel('B/*/[1-98, 109-120]')" },    // heavy chain minus CDRH3: fixed
  { "select": "sel('C')" }, { "select": "sel('D')" }   // HER2: fixed
]}
```

To change what is redesigned, change that gap:

- **whole heavy chain** — replace the chain B line with nothing, or `sel('B/*/[1-0]')`
- **several CDRs** — list the framework ranges you want kept, e.g.
  `sel('B/*/[1-25, 36-49, 59-98, 109-120]')` leaves three loops to decode
- **the light chain too** — narrow the `sel('A')` line the same way

## The other conditions

- **`C_CRDGRP`** holds coordinates. With inverse folding nothing moves, so this is
  bookkeeping: group 1 is the antibody, group 2 the antigen backbone.
- **`C_HOT`** is optional here. The decoder is conditioned on the structure it is given,
  so hotspots matter less than in a diffusion run; keep them to bias the decoder toward
  the interface.

## Notes

- **No comments in the JSON.** Every top-level key is read as a design id, so a
  `"_comment"` entry fails validation with `Invalid RFD4-Proteina design input`.
- Residue numbers are the file's own (`label_seq_id`), not what a viewer shows.
- Decoded sequences land in the generated `.cif` like any other run; read them with
  `python -m abhack.designs collect runs/inverse-s<seed>`.
- The designs still have to clear `abguard`: decoding onto trastuzumab's backbone can
  easily produce a CDRH3 too close to a known one.
