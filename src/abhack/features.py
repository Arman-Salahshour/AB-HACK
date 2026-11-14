"""Sequence featurisation for the CDRH3 binder/affinity models.

The FACS and SPR tables vary a fixed 10-residue window, so a positional one-hot
is a strong, cheap baseline. ``aa_property_matrix`` adds physicochemical columns
for models that benefit from them.
"""

from __future__ import annotations

import numpy as np

from abhack.constants import AMINO_ACIDS, TRASTUZUMAB_CDRH3_WINDOW

AA_INDEX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

# Kyte-Doolittle hydropathy, volume (A^3), charge at pH 7, aromaticity.
AA_PROPERTIES: dict[str, tuple[float, float, float, float]] = {
    "A": (1.8, 88.6, 0, 0), "C": (2.5, 108.5, 0, 0), "D": (-3.5, 111.1, -1, 0),
    "E": (-3.5, 138.4, -1, 0), "F": (2.8, 189.9, 0, 1), "G": (-0.4, 60.1, 0, 0),
    "H": (-3.2, 153.2, 0.1, 1), "I": (4.5, 166.7, 0, 0), "K": (-3.9, 168.6, 1, 0),
    "L": (3.8, 166.7, 0, 0), "M": (1.9, 162.9, 0, 0), "N": (-3.5, 114.1, 0, 0),
    "P": (-1.6, 112.7, 0, 0), "Q": (-3.5, 143.8, 0, 0), "R": (-4.5, 173.4, 1, 0),
    "S": (-0.8, 89.0, 0, 0), "T": (-0.7, 116.1, 0, 0), "V": (4.2, 140.0, 0, 0),
    "W": (-0.9, 227.8, 0, 1), "Y": (-1.3, 193.6, 0, 1),
}


def one_hot(sequences: list[str], length: int | None = None) -> np.ndarray:
    """``(n, length * 20)`` positional one-hot. Shorter sequences are zero-padded."""
    length = length or max(len(s) for s in sequences)
    out = np.zeros((len(sequences), length, len(AMINO_ACIDS)), dtype=np.float32)
    for i, seq in enumerate(sequences):
        for j, aa in enumerate(seq[:length]):
            k = AA_INDEX.get(aa)
            if k is not None:
                out[i, j, k] = 1.0
    return out.reshape(len(sequences), -1)


def aa_property_matrix(sequences: list[str], length: int | None = None) -> np.ndarray:
    """``(n, length * 4)`` physicochemical features."""
    length = length or max(len(s) for s in sequences)
    out = np.zeros((len(sequences), length, 4), dtype=np.float32)
    for i, seq in enumerate(sequences):
        for j, aa in enumerate(seq[:length]):
            out[i, j] = AA_PROPERTIES.get(aa, (0.0, 0.0, 0.0, 0.0))
    return out.reshape(len(sequences), -1)


def mutation_features(sequences: list[str], reference: str = TRASTUZUMAB_CDRH3_WINDOW) -> np.ndarray:
    """``(n, 2)``: number of mutations from the reference, and fraction mutated."""
    out = np.zeros((len(sequences), 2), dtype=np.float32)
    for i, seq in enumerate(sequences):
        n = sum(a != b for a, b in zip(seq, reference)) + abs(len(seq) - len(reference))
        out[i] = (n, n / max(len(reference), 1))
    return out


def featurise(sequences: list[str], *, properties: bool = True, length: int | None = None):
    """One-hot (+ properties + mutation count) feature matrix."""
    blocks = [one_hot(sequences, length), mutation_features(sequences)]
    if properties:
        blocks.insert(1, aa_property_matrix(sequences, length))
    return np.hstack(blocks)
