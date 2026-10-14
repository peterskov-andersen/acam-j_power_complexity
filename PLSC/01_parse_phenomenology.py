"""Step 1 — parse the raw ACAM-J rating files into one row per participant x state.
Input : inputs/ratings_raw/*Jhana*_ratings.txt
Output: outputs/intermediate/phenomenology_parsed.csv, outputs/tables/dropped_out_of_scale_ratings.csv"""
import re, pandas as pd
from config import RATINGS, INTERIM, TABLES, EXCLUDE_RATING_FILES, SCALE_MAX

SCORE = re.compile(r"-\s*([0-9]+(?:\.[0-9]+)?)\s*$")
ORD = [(1, "1st"), (2, "2nd"), (3, "3rd"), (4, "4th"), (5, "5th"), (6, "6th"), (7, "7th"), (8, "8th")]
SENS = {"sights": "sights", "sounds": "sounds", "physical sensations": "physical", "narrative thought stream": "narrative"}
RUNQ = {"overall quality": "quality_overall", "sense of unity": "quality_sense_of_unity",
        "sense of self": "quality_sense_of_self", "luminosity": "quality_luminosity"}

def state_no(label):
    low = label.lower()
    for n, p in ORD:
        if re.search(p + r"\s+jhana", low): return n
    return None

def classify(label):
    """Map an item label to [(variable, [ACAM-J states it applies to])]."""
    L = label.strip(); low = L.lower(); n = state_no(L)
    per_state = {"stabilityofattention": "stability_attention", "widthofattention": "width_attention",
                 "intensity-": "jhana_intensity", "grades of formlessness": "formlessness",
                 "bliss/joy": "bliss_joy", "automatic movements": "automatic_movements",
                 "energetic sensations in the body": "energetic_sensations"}
    for key, var in per_state.items():
        if low.startswith(key): return [(var, [n])] if n else []
    if low.startswith("jhanic neutral qualities"): return [("jhanic_neutral_qualities", [n or 4])]
    if low.startswith("contentment"): return [("contentment", [n or 3])]
    # both spellings ("senstions"/"sensations", "qualilty"/"quality") occur in the rating files
    for pre, block, states in [(("senstionsearly-", "sensationsearly-"), "early", [1, 2, 3, 4]),
                               (("senstionslate-", "sensationslate-"), "late", [5, 6, 7, 8])]:
        if low.startswith(pre):
            t = L.split("-", 1)[1].strip().lower()
            return [(f"sensations_{block}_{SENS[t]}", states)] if t in SENS else []
    if low.startswith("qualilty-") or low.startswith("quality-"):
        r = L.split("-", 1)[1].strip().lower()
        return [(RUNQ[r], list(range(1, 9)))] if r in RUNQ else []
    return []

rows, dropped, excluded = [], [], []
for f in sorted(RATINGS.glob("*Jhana*_ratings.txt")):
    if f.name in EXCLUDE_RATING_FILES: excluded.append(f.name); continue
    for line in open(f, encoding="utf8", errors="ignore"):
        parts = line.rstrip("\n").rstrip("\r").split("\t")
        if len(parts) < 2: continue
        idm = re.search(r"\d+", parts[0]); entry = parts[1].strip()
        if not idm or not entry: continue
        sm = SCORE.search(entry)
        if not sm: continue
        subj = f"sub-{int(idm.group()):02d}"; score = float(sm.group(1)); label = entry[:sm.start()]
        items = classify(label)
        if not items: continue
        if score > SCALE_MAX:      # NA code -> missing, dropped before averaging runs
            dropped.append((f.name, subj, entry)); continue
        for var, states in items:
            for st in states:
                if st: rows.append((subj, f"ACAM{st}", var, score))

long = pd.DataFrame(rows, columns=["subject", "state", "variable", "score"])
wide = (long.groupby(["subject", "state", "variable"], as_index=False)["score"].mean()
            .pivot(index=["subject", "state"], columns="variable", values="score").reset_index())
wide.columns.name = None
wide.to_csv(INTERIM / "phenomenology_parsed.csv", index=False)
pd.DataFrame(dropped, columns=["file", "subject", "entry"]).to_csv(TABLES / "dropped_out_of_scale_ratings.csv", index=False)
print(f"[01] excluded files: {excluded}")
print(f"[01] dropped {len(dropped)} out-of-scale ratings before averaging")
print(f"[01] parsed {wide.shape[0]} participant x state rows, {wide.shape[1]-2} variables")
