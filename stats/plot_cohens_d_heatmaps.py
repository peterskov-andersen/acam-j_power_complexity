import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend — no GUI window needed
import matplotlib.pyplot as plt
import os

plt.rcParams["font.family"] = "Times New Roman"

base_dir = os.path.dirname(os.path.abspath(__file__))

bands = ["delta", "theta", "alpha", "low_beta", "high_beta"]
band_labels = ["Delta", "Theta", "Alpha", "Low Beta", "High Beta"]
categories = ["total", "oscillatory", "aperiodic"]
cat_labels = ["Total", "Oscillatory", "Aperiodic"]
col_labels = [f"ACAM-J{i}" for i in range(1, 9)]


def load_d_values(category, band):
    fname = f"LMM_{category}_{band}_5s_ACAM_vs_cm_FDR_BH_with_d.csv"
    path = os.path.join(base_dir, fname)
    df = pd.read_csv(path)
    df["acam_num"] = df["contrast"].str.extract(r"ACAM(\d+)").astype(int)
    df = df.sort_values("acam_num")
    return df["d"].values


# Build matrices: shape (5 bands × 8 states)
matrices = {}
for cat in categories:
    mat = np.zeros((len(bands), 8))
    for i, band in enumerate(bands):
        mat[i, :] = load_d_values(cat, band)
    matrices[cat] = mat

# Shared colour scale across all three matrices
all_vals = np.concatenate([m.ravel() for m in matrices.values()])
vmax = np.nanmax(np.abs(all_vals))
vmin = -vmax

cmap = "RdBu_r"  # blue (negative) → white → red (positive)

fig, axes = plt.subplots(
    1, 3,
    figsize=(20, 7),
    gridspec_kw={"wspace": 0.4}
)

ims = []
for ax, cat, label in zip(axes, categories, cat_labels):
    im = ax.imshow(
        matrices[cat],
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest"
    )
    ims.append(im)

    ax.set_xticks(range(8))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=12, fontweight="bold")
    ax.set_yticks(range(len(bands)))
    ax.set_yticklabels(band_labels, fontsize=12, fontweight="bold")
    ax.set_title(label, fontsize=15, fontweight="bold", pad=10)

    # Annotate each cell with the d value
    for row in range(len(bands)):
        for col in range(8):
            val = matrices[cat][row, col]
            text_color = "white" if abs(val) > 0.6 * vmax else "black"
            ax.text(col, row, f"{val:.2f}", ha="center", va="center",
                    fontsize=10, color=text_color, fontweight="bold")

# Single shared colorbar on the right
cbar = fig.colorbar(ims[0], ax=axes, orientation="vertical",
                    fraction=0.02, pad=0.04, shrink=0.85)
cbar.set_label("Cohen's d", fontsize=13, fontweight="bold")
cbar.ax.tick_params(labelsize=12)
for label in cbar.ax.get_yticklabels():
    label.set_fontweight("bold")

fig.suptitle("Cohen's d: ACAM Jhana States vs Control Meditation",
             fontsize=15, fontweight="bold", y=1.02)

out_path = os.path.join(base_dir, "cohens_d_heatmaps.png")
plt.savefig(out_path, dpi=200, bbox_inches="tight")
print(f"Saved to: {out_path}")
