"""Fixed facts about trastuzumab, HER2 and the hackathon filters."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"
DOCS = ROOT / "data" / "docs"
RUNS = ROOT / "runs"
DESIGNS = ROOT / "designs"

FACS_CSV = DATA_RAW / "her2-trastuzumab" / "tz_her2_facs_3classes_524346.csv.zip"
SPR_CSV = DATA_RAW / "her2-trastuzumab" / "trastuzumab_her2_spr_all_768.csv"

# --- trastuzumab (the baseline every design is compared against) -------------
# Heavy and light variable domains, as they appear in the hackathon SPR table.
TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)

#: The 10-residue window the FACS screen mutated. This is what the hackathon
#: tables call ``cdrh3``; it is the IMGT CDR3 without its ``SR``/``Y`` flanks.
TRASTUZUMAB_CDRH3_WINDOW = "WGGDGFYAMD"
#: Index of that window inside :data:`TRASTUZUMAB_VH`.
CDRH3_WINDOW_START = TRASTUZUMAB_VH.index(TRASTUZUMAB_CDRH3_WINDOW)
CDRH3_WINDOW_END = CDRH3_WINDOW_START + len(TRASTUZUMAB_CDRH3_WINDOW)

# HER2 structure used as the design template (Fv + HER2 domain IV, cropped).
HER2_TEMPLATE_CIF = DOCS / "Data" / "12-her2" / "her2_fv_1n8z.cif"
HER2_DESIGN_YAML = DOCS / "Data" / "12-her2" / "her2_cdrh3.yaml"
HER2_PDB_ID = "1N8Z"

# --- hard filters ------------------------------------------------------------
MAX_TOTAL_LENGTH = 250  # aa, single chain or VH+VL combined
MIN_NUMBERING_CONFIDENCE = 0.1
NOVELTY_IDENTITY_THRESHOLD = 0.70  # below this -> de novo track
#: CDRH3 mutations required away from every provided sequence.
MIN_CDRH3_MUTATIONS_OPTIMISATION = 3
MIN_CDRH3_MUTATIONS_DE_NOVO = 2

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# --- submission table --------------------------------------------------------
SUBMISSION_COLUMNS = [
    "sequence_id",
    "kind",  # vhh | scfv | paired
    "heavy_sequence",
    "light_sequence",
    "design_description",
]

# Common scFv linker, if a paired design is submitted as a single chain.
SCFV_LINKER = "GGGGSGGGGSGGGGS"
