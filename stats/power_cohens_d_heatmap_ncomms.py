"""Cohen's d heatmaps of global power change (ACAM-J1-8 vs combined control),
Nature Communications version.

Reads the per-measure LMM results written by entropy_lmm_and_figures_ncomms.R
(per_measure/LMM_<measure>_ACAM_vs_control_FDR_BH_with_d.csv) and draws one
heatmap per power component (total, oscillatory, aperiodic): rows = frequency
bands, columns = ACAM-J states, cell text = Cohen's d, asterisks = BH-FDR-corrected
P across the eight contrasts within each measure.

Layout: 183 mm wide (double column), sans-serif (Arial if installed), 5-7 pt text,
8 pt bold lowercase panel letters, shared symmetric colour scale (blue = lower
power during ACAM-J, red = higher). Saves PDF (editable text), PNG and TIFF.

Place this script in the R output folder (entropy_outputs_ncomms) or set RESULTS_DIR.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(BASE_DIR, "per_measure")   # where the R script wrote the per-measure CSVs
OUT_STEM    = os.path.join(BASE_DIR, "FIG_power_cohens_d_heatmap")

MM = 1 / 25.4
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "mathtext.fontset": "custom", "mathtext.rm": "sans", "mathtext.it": "sans:italic",
    "font.size": 6.5, "axes.titlesize": 7, "axes.titleweight": "normal",
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.facecolor": "white", "savefig.facecolor": "white",
})

bands       = ["delta", "theta", "alpha", "low_beta", "high_beta"]
band_labels = ["Delta", "Theta", "Alpha", "Low beta", "High beta"]
categories  = ["total", "oscillatory", "aperiodic"]
cat_labels  = ["Total power", "Oscillatory power", "Aperiodic power"]
col_labels  = [f"J{i}" for i in range(1, 9)]          # legend: J1-J8, ACAM-J1-8


def load(category, band):
    path = os.path.join(RESULTS_DIR, f"LMM_{category}_{band}_5s_ACAM_vs_control_FDR_BH_with_d.csv")
    df = pd.read_csv(path)
    df["acam_num"] = df["contrast"].str.extract(r"ACAM(\d+)", expand=False).astype(int)
    df = df.set_index("acam_num").reindex(range(1, 9))   # keeps a column even if a state is missing
    return df["d"].to_numpy(), df["p_fdr_bh"].to_numpy()


D = {c: np.full((len(bands), 8), np.nan) for c in categories}
Q = {c: np.full((len(bands), 8), np.nan) for c in categories}
for c in categories:
    for i, b in enumerate(bands):
        D[c][i], Q[c][i] = load(c, b)

# Shared symmetric colour scale, rounded up to the next 0.5
vmax = np.ceil(np.nanmax(np.abs(np.concatenate([m.ravel() for m in D.values()]))) * 2) / 2
star = lambda q: "" if np.isnan(q) else "***" if q < .001 else "**" if q < .01 else "*" if q < .05 else ""

fig = plt.figure(figsize=(183 * MM, 62 * MM))
gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.05], wspace=0.28,
                      left=0.07, right=0.94, top=0.86, bottom=0.12)
axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
cax = fig.add_subplot(gs[0, 3])

for k, (ax, c, title) in enumerate(zip(axes, categories, cat_labels)):
    im = ax.imshow(D[c], aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax, interpolation="nearest")
    ax.set_xticks(range(8)); ax.set_xticklabels(col_labels)
    ax.set_yticks(range(len(bands))); ax.set_yticklabels(band_labels if k == 0 else [])
    ax.set_title(title, pad=4)
    ax.tick_params(length=0)
    # white grid between cells
    ax.set_xticks(np.arange(-.5, 8, 1), minor=True); ax.set_yticks(np.arange(-.5, len(bands), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.8); ax.tick_params(which="minor", length=0)
    for sp in ax.spines.values(): sp.set_visible(False)

    for r in range(len(bands)):
        for col in range(8):
            v, q = D[c][r, col], Q[c][r, col]
            if np.isnan(v):
                continue
            tc = "white" if abs(v) > 0.55 * vmax else "0.1"
            s = star(q)
            txt = f"{v:.2f}".replace("-0.00", "0.00")
            ax.text(col, r - (0.13 if s else 0), txt, ha="center", va="center", fontsize=5, color=tc)
            if s:
                ax.text(col, r + 0.27, s, ha="center", va="center", fontsize=5.5, color=tc)

cb = fig.colorbar(im, cax=cax)
cb.set_label(r"Cohen's $d$ (ACAM-J $-$ control)", fontsize=6.5)
cb.ax.tick_params(labelsize=6, length=2, width=0.5)
cb.outline.set_visible(False)

fig.canvas.draw()
for ax, letter in zip(axes, "abc"):
    bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(max(bb.x0 - 0.01, 0.002), 0.985, letter, fontsize=8, fontweight="bold", va="top", ha="left")

fig.savefig(OUT_STEM + ".pdf")
fig.savefig(OUT_STEM + ".png", dpi=600)
fig.savefig(OUT_STEM + ".tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
print(f"Saved {OUT_STEM}.pdf/.png/.tiff (183 x 62 mm, colour scale +/-{vmax})")
