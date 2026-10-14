"""Shared settings for the neurophenomenology pipeline. All paths are relative to the folder
above this one (inputs/ and outputs/), so the folder can be moved anywhere.
Run 01-05 in order from this folder."""
from pathlib import Path

ROOT     = Path(__file__).resolve().parents[1]
INPUTS   = ROOT / "inputs"
RATINGS  = INPUTS / "ratings_raw"
EEG_CSV  = INPUTS / "eeg_measures.csv"
OUT      = ROOT / "outputs"
TABLES   = OUT / "tables"
FIGS     = OUT / "figures"
INTERIM  = OUT / "intermediate"
for d in (TABLES, FIGS, INTERIM): d.mkdir(parents=True, exist_ok=True)

# ---------------- phenomenology parsing rules ----------------
# sub-035 completed two ACAM-J sessions; only run 2 has retained EEG, so run 1 is excluded.
EXCLUDE_RATING_FILES = {"Adept-035-Jhana-1_ratings.txt"}
# Scores above the scale maximum are 'not entered / not scored' codes -> missing.
# They are dropped BEFORE averaging across runs (otherwise they contaminate valid run values).
SCALE_MAX = 10

# ---------------- EEG measures (the 16 measures in the PLSC) ----------------
EEG = [  # (label used in figures, column in eeg_measures.csv)
    ("LZsum", "LZsum_c0"), ("LZsumN", "LZsumPN"), ("HFD", "HFD_kmax200"), ("DFA", "DFA_alpha"),
    ("SpecEn", "SpecEn_1to35Hz"), ("PeEn", "PE_1_35_tau4_s2"),
    ("aper delta", "aperiodic_delta"), ("aper theta", "aperiodic_theta"), ("aper alpha", "aperiodic_alpha"),
    ("aper low beta", "aperiodic_low_beta"), ("aper high beta", "aperiodic_high_beta"),
    ("osc delta", "oscillatory_delta"), ("osc theta", "oscillatory_theta"), ("osc alpha", "oscillatory_alpha"),
    ("osc low beta", "oscillatory_low_beta"), ("osc high beta", "oscillatory_high_beta"),
]
XLAB = [a for a, _ in EEG]; XCOL = [b for _, b in EEG]

# ---------------- phenomenology ----------------
# State-level items entering the main PLSC (per-absorption items + block-level content items).
Y = [("quality of ACAM-J", "jhana_intensity"), ("stability of attention", "stability_attention"),
     ("width of attention", "width_attention"), ("sights", "sensations_sights"), ("sounds", "sensations_sounds"),
     ("physical sensations", "sensations_physical"), ("narrative thought", "sensations_narrative")]
YLAB = [a for a, _ in Y]; YCOL = [b for _, b in Y]
# Run-level items -> separate subject-level PLSC.
Q = [("overall quality", "quality_overall"), ("sense of unity", "quality_sense_of_unity"),
     ("sense of self", "quality_sense_of_self"), ("luminosity", "quality_luminosity")]
QLAB = [a for a, _ in Q]; QCOL = [b for _, b in Q]
# ACAM-J-specific factors -> univariate Spearman, FDR within each factor.
FACTORS = [("Bliss/joy\n(ACAM-J 1–2)", "bliss_joy", [1, 2]), ("Contentment\n(ACAM-J 3)", "contentment", [3]),
           ("Equanimity\n(ACAM-J 4)", "jhanic_neutral_qualities", [4]),
           ("Formlessness\n(ACAM-J 5–8)", "formlessness", [5, 6, 7, 8])]

# ---------------- statistics ----------------
N_PERM = 10000
N_BOOT = 10000
SEED = 999
IMPUTE_ITER = 20          # MICE (IterativeImputer) iterations
IMPUTE_SEED = 0
# The original analysis table stored EEG values rounded to 4 decimals; kept for exact reproduction.
ROUND_EEG_DECIMALS = 4
