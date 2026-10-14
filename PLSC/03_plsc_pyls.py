"""Step 3 — behavioural PLSC with the pyls package (netneurolab/pypyls v0.0.2).
  A) between-participant : observed values, unrestricted permutation, bootstrap over observations
  B) within-participant  : X and Y centred within participant, permutation within participant,
                           bootstrap over participants (custom samples passed to pyls)
  C) run-level quality   : EEG averaged across ACAM-J states, one row per participant (n = 20)
Significance: permutation p < .05 AND Kaiser criterion (variance > 100 / n_components).
NOTE: pyls v0.0.2 needs permindices=True; without it, internally generated permutation indices
are passed in place of Y and the permutation step fails."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyls
from scipy.stats import zscore
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
import sys
from config import (INTERIM, TABLES, XCOL, XLAB, YCOL, YLAB, QCOL, QLAB, N_PERM, N_BOOT, SEED,
                    IMPUTE_ITER, IMPUTE_SEED, ROUND_EEG_DECIMALS)

t = pd.read_csv(INTERIM / "analysis_table.csv")
t[XCOL] = t[XCOL].round(ROUND_EEG_DECIMALS)

def save(tag, r, ylab, extra=None):
    ve, p = r["varexp"], r["permres"]["pvals"]; br = r["bootres"]
    np.savez(INTERIM / f"pyls_{tag}.npz", ve=ve, p=p, bsr=br["x_weights_normed"], yl=br["y_loadings"],
             ci=br["y_loadings_ci"], xs=r["x_scores"], ys=r["y_scores"])
    k = len(ve); kai = 100 / k; lv = [f"LV{i+1}" for i in range(k)]
    tab = pd.DataFrame({"LV": range(1, k + 1), "var_%": (ve * 100).round(2), "perm_p": np.round(p, 4),
                        "kaiser_pass": ve * 100 > kai, "significant": (p < .05) & (ve * 100 > kai)})
    if extra: tab.insert(1, "n_obs", extra)
    tab.to_csv(TABLES / f"{tag}_variance_explained.csv", index=False)
    pd.DataFrame(br["x_weights_normed"], index=XLAB, columns=lv).to_csv(TABLES / f"{tag}_X_bootstrap_ratios.csv")
    pd.DataFrame(br["y_loadings"], index=ylab, columns=lv).to_csv(TABLES / f"{tag}_Y_loadings.csv")
    lo, hi = br["y_loadings_ci"][:, :, 0], br["y_loadings_ci"][:, :, 1]
    pd.DataFrame(lo, index=ylab, columns=lv).to_csv(TABLES / f"{tag}_Y_loadings_CI_lower.csv")
    pd.DataFrame(hi, index=ylab, columns=lv).to_csv(TABLES / f"{tag}_Y_loadings_CI_upper.csv")
    sig = [f"LV{i+1} ({ve[i]*100:.1f}%, p={p[i]:.4f})" for i in range(k) if (p[i] < .05) and (ve[i] * 100 > kai)]
    print(f"[03] {tag:8s}: significant -> {', '.join(sig) if sig else 'none'}")

# ---- state-level data: impute missing phenomenology (MICE), then drop rows with missing EEG ----
d = t.copy()
d[YCOL] = IterativeImputer(max_iter=IMPUTE_ITER, random_state=IMPUTE_SEED).fit_transform(d[YCOL])
d = d.dropna(subset=XCOL).reset_index(drop=True)

# Which analyses to run: default all; e.g. `python 03_plsc_pyls.py within` runs one.
RUN = set(sys.argv[1:]) or {"between", "within", "quality"}

# A) between-participant
if "between" in RUN:
    X = zscore(d[XCOL].to_numpy(), ddof=1); Y = zscore(d[YCOL].to_numpy(), ddof=1)
    rb = pyls.behavioral_pls(X, Y, n_perm=N_PERM, n_boot=N_BOOT, test_split=0, seed=SEED, verbose=False,
                             permindices=True)
    save("between", rb, YLAB, extra=len(d))

# B) within-participant
if "within" in RUN:
    c = d.groupby("subid").size(); dw = d[d.subid.isin(c[c > 1].index)].reset_index(drop=True)
    g = dw["subid"].values; idx = {s: np.where(g == s)[0] for s in np.unique(g)}
    def centre(A):
        A = A.astype(float).copy()
        for ix in idx.values(): A[ix] -= A[ix].mean(0)
        return A
    Xw = zscore(centre(dw[XCOL].to_numpy()), ddof=1); Yw = zscore(centre(dw[YCOL].to_numpy()), ddof=1)
    rng = np.random.default_rng(SEED); n = len(dw); subs = np.unique(g)
    perms = np.zeros((n, N_PERM), dtype=int)
    for j in range(N_PERM):
        pm = np.arange(n)
        for ix in idx.values(): pm[ix] = rng.permutation(ix)
        perms[:, j] = pm
    # bootstrap resamples participants; rows are truncated or repeated to keep n fixed (pyls requirement)
    boots = np.zeros((n, N_BOOT), dtype=int)
    for j in range(N_BOOT):
        rows = np.concatenate([idx[s] for s in rng.choice(subs, size=len(subs), replace=True)])
        boots[:, j] = rows[:n] if len(rows) >= n else np.pad(rows, (0, n - len(rows)), mode="wrap")
    rw = pyls.behavioral_pls(Xw, Yw, n_perm=N_PERM, n_boot=N_BOOT, test_split=0, seed=SEED, verbose=False,
                             permindices=True, permsamples=perms, bootsamples=boots)
    save("within", rw, YLAB, extra=len(dw))

# C) run-level quality: one row per participant, EEG averaged across states
if "quality" in RUN:
    gq = t.groupby("subid"); dq = pd.concat([gq[XCOL].mean(), gq[QCOL].first()], axis=1).dropna()
    Xq = zscore(dq[XCOL].to_numpy(), ddof=1); Yq = zscore(dq[QCOL].to_numpy(), ddof=1)
    rq = pyls.behavioral_pls(Xq, Yq, n_perm=N_PERM, n_boot=N_BOOT, test_split=0, seed=SEED, verbose=False,
                             permindices=True)
    save("quality", rq, QLAB, extra=len(dq))
