"""Step 2 — merge EEG measures with parsed phenomenology (ACAM-J rows only).
Block-level content items: the 'early' rating is assigned to ACAM-J1-4, the 'late' rating to ACAM-J5-8.
Output: outputs/intermediate/analysis_table.csv"""
import pandas as pd, numpy as np
from config import EEG_CSV, INTERIM, XCOL

eeg = pd.read_csv(EEG_CSV)
ph = pd.read_csv(INTERIM / "phenomenology_parsed.csv")
acam = eeg[eeg.state.str.startswith("ACAM")].copy()            # original row order preserved
acam["stage"] = acam.state.str.replace("ACAM", "").astype(int)
acam = acam.merge(ph, on=["subject", "state"], how="left", validate="one_to_one")
for t in ["sights", "sounds", "physical", "narrative"]:
    acam[f"sensations_{t}"] = np.where(acam.stage <= 4, acam[f"sensations_early_{t}"], acam[f"sensations_late_{t}"])
acam.insert(1, "subid", acam.subject.map(lambda s: f"sub-{int(s.replace('sub-','')):03d}"))
acam.to_csv(INTERIM / "analysis_table.csv", index=False)
print(f"[02] analysis table: {len(acam)} rows ({acam.subject.nunique()} participants x 8 ACAM-J states)")
print(f"[02] rows with any missing EEG measure: {int(acam[XCOL].isna().any(axis=1).sum())}")
