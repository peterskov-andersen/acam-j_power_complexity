"""Step 4 — Spearman correlations between each EEG measure and each ACAM-J-specific factor.
Score and EEG are averaged over the states where the factor is rated (one value per participant).
Benjamini-Hochberg FDR within each factor, across all 16 EEG measures; never pooled across factors."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from config import INTERIM, TABLES, EEG, FACTORS

t = pd.read_csv(INTERIM / "analysis_table.csv")
def bh(p):
    p = np.asarray(p, float); q = np.full_like(p, np.nan); ix = np.where(~np.isnan(p))[0]; m = len(ix)
    o = ix[np.argsort(p[ix])]; qv = p[o] * m / np.arange(1, m + 1); qv = np.minimum.accumulate(qv[::-1])[::-1]
    for i, j in enumerate(o): q[j] = min(qv[i], 1.0)
    return q
recs = []
for flab, fcol, states in FACTORS:
    g = t[t.stage.isin(states)].groupby("subject"); y = g[fcol].mean(); block = []
    for lab, col in EEG:
        dd = pd.concat([g[col].mean(), y], axis=1).dropna(); dd.columns = ["x", "y"]
        rho, p = spearmanr(dd.x, dd.y)
        block.append(dict(factor=flab.split("\n")[0], brain=lab, rho=rho, p=p, n=len(dd)))
    q = bh([b["p"] for b in block])
    for b, qq in zip(block, q): b["q_fdr"] = qq
    recs += block
res = pd.DataFrame(recs)[["factor", "brain", "rho", "p", "q_fdr", "n"]]
res.to_csv(TABLES / "univariate_factor_results.csv", index=False)
print(f"[04] {len(res)} correlations; smallest raw p = {res.p.min():.3f}; smallest q = {res.q_fdr.min():.3f}; "
      f"q < .05: {int((res.q_fdr < .05).sum())}")
