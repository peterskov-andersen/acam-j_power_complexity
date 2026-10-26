"""Step 6 — Nature Communications versions of the PLSC figures.

Same content as 05_figures.py, re-laid out to the journal's figure specifications:
  * 183 mm wide (double column), sans-serif text at 5-7 pt, panel letters 8 pt bold lowercase (a-f)
  * non-bold titles, thinner lines, italic r and P
  * text kept as editable text in the PDF (TrueType embedded)
Uses Arial when it is installed (e.g. on Windows); otherwise falls back to Liberation Sans
(metrically identical to Arial) or DejaVu Sans.

Outputs (in outputs/figures):
  Figure_PLSC_main_ncomms.pdf / .png / .tiff
  Figure_S_factors_and_quality_ncomms.pdf / .png / .tiff
Run from the scripts folder:  python 06_figures_ncomms.py
"""
import numpy as np, pandas as pd
from scipy.stats import pearsonr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from config import INTERIM, TABLES, FIGS, XLAB, YLAB, QLAB, EEG, FACTORS

MM = 1 / 25.4
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    'mathtext.fontset': 'custom', 'mathtext.rm': 'sans', 'mathtext.it': 'sans:italic',
    'mathtext.bf': 'sans:bold',
    'font.size': 6.5, 'axes.labelsize': 7, 'axes.titlesize': 7, 'axes.titleweight': 'normal',
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'axes.linewidth': 0.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2, 'ytick.major.size': 0,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'figure.facecolor': 'white', 'savefig.facecolor': 'white'})
LETTER = dict(fontsize=8, fontweight='bold', va='top', ha='left')

RED, BLUE, NS = '#B40426', '#3B4CC0', '#D9D9D9'; ERED, EBLUE, ENS = '#7d0219', '#25307d', '#A6A6A6'
CMAP = LinearSegmentedColormap.from_list('c', [BLUE, '#8A9BE0', '#C9D2F0', '#FFFFFF', '#F2C4C0', '#D2604F', RED])
FLIP = {"between": True, "within": True}   # same orientation convention as 05_figures.py

# Readable labels (figure only; the analysis labels in config.py are unchanged)
def pretty_eeg(s):
    s = s.replace('aper ', 'Aperiodic ').replace('osc ', 'Oscillatory ')
    return r'DFA $\alpha$' if s == 'DFA' else s
cap = lambda s: s[0].upper() + s[1:]
XL = [pretty_eeg(x) for x in XLAB]
YL = [cap(y) for y in YLAB]
QL = [cap(q) for q in QLAB]

def p_txt(p):
    return r'$P$ < 0.001' if p < 0.001 else rf'$P$ = {p:.3f}'

def clean(ax):
    for k, sp in ax.spines.items(): sp.set_visible(k in ('left', 'bottom'))

def panelset(axs, z, lv, flip, title, scatter_x='Brain latent score', sig=True):
    s = -1.0 if flip else 1.0
    xs = z['xs'][:, lv] * s; ys = z['ys'][:, lv] * s; bsr = z['bsr'][:, lv] * s; ylo = z['yl'][:, lv] * s
    if flip: lo, hi = -z['ci'][:, lv, 1], -z['ci'][:, lv, 0]
    else:    lo, hi =  z['ci'][:, lv, 0],  z['ci'][:, lv, 1]

    ax = axs[0]; rr, pp = pearsonr(xs, ys)
    ax.scatter(xs, ys, s=7, facecolor=RED, edgecolor='white', lw=.25, alpha=.8)
    m, b = np.polyfit(xs, ys, 1); xx = np.linspace(xs.min(), xs.max(), 50)
    ax.plot(xx, m * xx + b, color='0.15', lw=0.9)
    ax.text(.03, .97, rf'$r$ = {rr:.2f}, ' + p_txt(pp), transform=ax.transAxes, fontsize=6, va='top')
    ax.set_xlabel(scatter_x); ax.set_ylabel('Phenomenology latent score'); ax.set_title(title, loc='center')
    clean(ax)

    ax = axs[1]; o = np.argsort(bsr)
    ax.barh(range(len(bsr)), bsr[o], height=.72,
            color=[RED if v > 1.96 else BLUE if v < -1.96 else NS for v in bsr[o]],
            edgecolor=[ERED if v > 1.96 else EBLUE if v < -1.96 else ENS for v in bsr[o]], lw=.35)
    ax.set_yticks(range(len(bsr))); ax.set_yticklabels(np.array(XL)[o])
    for x in (1.96, -1.96): ax.axvline(x, ls='--', c='0.4', lw=.5)
    ax.axvline(0, c='k', lw=.5)
    ax.set_xlabel('Bootstrap ratio'); ax.set_title('EEG measures'); clean(ax)

    ax = axs[2]; o = np.argsort(ylo); rel = (lo * hi > 0) & sig
    ax.barh(range(len(ylo)), ylo[o], xerr=[ylo[o] - lo[o], hi[o] - ylo[o]], height=.65,
            color=[(RED if ylo[i] > 0 else BLUE) if rel[i] else NS for i in o],
            edgecolor=[(ERED if ylo[i] > 0 else EBLUE) if rel[i] else ENS for i in o], lw=.35,
            error_kw=dict(ecolor='black', lw=.6, capsize=1.5, capthick=.6))
    ax.set_yticks(range(len(ylo))); ax.set_yticklabels(np.array(YL)[o])
    ax.axvline(0, c='k', lw=.5); ax.set_xlabel('Loading ($r$ with brain latent score)')
    ax.set_title('Phenomenology'); clean(ax)

zb = np.load(INTERIM / 'pyls_between.npz'); zw = np.load(INTERIM / 'pyls_within.npz'); zq = np.load(INTERIM / 'pyls_quality.npz')

def save(fig, stem):
    fig.savefig(FIGS / f'{stem}.pdf')
    fig.savefig(FIGS / f'{stem}.png', dpi=600)
    fig.savefig(FIGS / f'{stem}.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})

# ---------------- main figure (183 x 140 mm) ----------------
fig = plt.figure(figsize=(183 * MM, 140 * MM))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1.05, 1.0], wspace=0.62, hspace=0.55,
                      top=0.93, bottom=0.08, left=0.06, right=0.985)
axs_b = [fig.add_subplot(gs[0, i]) for i in range(3)]
axs_w = [fig.add_subplot(gs[1, i]) for i in range(3)]
panelset(axs_b, zb, 0, FLIP['between'],
         f"Between-participant LV1\n{zb['ve'][0]*100:.1f}% of covariance, " + p_txt(zb['p'][0]))
panelset(axs_w, zw, 0, FLIP['within'],
         f"Within-participant LV1\n{zw['ve'][0]*100:.1f}% of covariance, " + p_txt(zw['p'][0]),
         scatter_x='Brain latent score (within-participant)')
fig.canvas.draw()
for axs, letters in ((axs_b, 'abc'), (axs_w, 'def')):
    for ax, l in zip(axs, letters):
        bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
        fig.text(max(bb.x0 - 0.005, 0.002), bb.y1 + 0.012, l, **LETTER)
save(fig, 'Figure_PLSC_main_ncomms'); plt.close(fig)

# ---------------- supplementary figure (183 x 105 mm) ----------------
u = pd.read_csv(TABLES / 'univariate_factor_results.csv')
flabs = [f[0].split('\n')[0] for f in FACTORS]; mlabs = [m[0] for m in EEG]
R = u.pivot(index='brain', columns='factor', values='rho').loc[mlabs, flabs].to_numpy()
Qv = u.pivot(index='brain', columns='factor', values='q_fdr').loc[mlabs, flabs].to_numpy()
star = lambda q: '***' if q < .001 else '**' if q < .01 else '*' if q < .05 else ''

fig, axes = plt.subplots(1, 2, figsize=(183 * MM, 105 * MM),
                         gridspec_kw={'width_ratios': [1.35, 1], 'wspace': 0.75})
fig.subplots_adjust(left=0.13, right=0.98, top=0.9, bottom=0.2)
ax = axes[0]; im = ax.imshow(R, cmap=CMAP, vmin=-.6, vmax=.6, aspect='auto')
ax.set_xticks(range(len(FACTORS))); ax.set_xticklabels([f[0] for f in FACTORS], rotation=90, va='top')
ax.set_yticks(range(len(EEG))); ax.set_yticklabels([pretty_eeg(m) for m in mlabs])
for i in range(len(EEG)):
    for j in range(len(FACTORS)):
        s = star(Qv[i, j])
        ax.text(j, i, f"{s}\n{R[i, j]:.2f}" if s else f"{R[i, j]:.2f}", ha='center', va='center',
                fontsize=5.5, linespacing=1.0, color='white' if abs(R[i, j]) > .42 else '0.1')
ax.axhline(5.5, color='k', lw=.7); ax.axhline(10.5, color='k', lw=.7)
ax.set_xticks(np.arange(-.5, len(FACTORS), 1), minor=True); ax.set_yticks(np.arange(-.5, len(EEG), 1), minor=True)
ax.grid(which='minor', color='white', lw=.8); ax.tick_params(which='both', length=0)
for sp in ax.spines.values(): sp.set_visible(False)
ax.set_title('ACAM-J–specific factors')
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04); cb.set_label(r"Spearman's $\rho$"); cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=6, length=0)

ax = axes[1]; ylq = zq['yl'][:, 0]; lo, hi = zq['ci'][:, 0, 0], zq['ci'][:, 0, 1]; o = np.argsort(ylq)
ax.barh(range(len(ylq)), ylq[o], xerr=[ylq[o] - lo[o], hi[o] - ylq[o]], height=.55, color=NS, edgecolor=ENS, lw=.35,
        error_kw=dict(ecolor='black', lw=.6, capsize=1.5, capthick=.6))
ax.set_yticks(range(len(ylq))); ax.set_yticklabels(np.array(QL)[o])
ax.axvline(0, c='k', lw=.5); ax.set_xlabel('Loading ($r$ with brain latent score)')
ax.set_title(f"Run-level quality items\nLV1: {zq['ve'][0]*100:.1f}% of covariance, " + p_txt(zq['p'][0]))
ax.text(.5, -0.2, 'LV not significant; loadings not interpreted', transform=ax.transAxes, ha='center',
        fontsize=6, style='italic', color='0.35')
clean(ax)
fig.canvas.draw()
for ax, l in zip(axes, 'ab'):
    bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(max(bb.x0 - 0.005, 0.002), min(bb.y1 + 0.03, 0.995), l, **LETTER)
save(fig, 'Figure_S_factors_and_quality_ncomms'); plt.close(fig)
print(f"[06] Nature Communications figures written to {FIGS}")
