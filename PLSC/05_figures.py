"""Step 5 — manuscript figures.
  Figure_PLSC_main            : between-participant LV1 (A-C) and within-participant LV1 (D-F)
  Figure_S_factors_and_quality: univariate factor heatmap (A) and run-level quality PLSC (B)
LV signs are arbitrary; both main LVs are flipped so the positive pole = lower spectral power
(the direction in which ACAM-J differs from the control condition)."""
import numpy as np, pandas as pd
from scipy.stats import pearsonr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from config import INTERIM, TABLES, FIGS, XLAB, YLAB, QLAB, EEG, FACTORS

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans','Arial'],
 'axes.labelsize':11,'axes.titlesize':11,'axes.titleweight':'bold','xtick.labelsize':10,
 'ytick.labelsize':10,'axes.linewidth':0.8,'figure.facecolor':'white','xtick.major.size':0,'ytick.major.size':0})
RED,BLUE,NS='#B40426','#3B4CC0','#D9D9D9'; ERED,EBLUE,ENS='#7d0219','#25307d','#A6A6A6'
CMAP=LinearSegmentedColormap.from_list('c',[BLUE,'#8A9BE0','#C9D2F0','#FFFFFF','#F2C4C0','#D2604F',RED])
FLIP = {"between": True, "within": True}   # orientation convention (see docstring)

def panelset(fig,gs,row,z,lv,flip,xl,yl_,title,scatter_x='Brain latent score',sig=True):
    s=-1.0 if flip else 1.0
    xs=z['xs'][:,lv]*s; ys=z['ys'][:,lv]*s; bsr=z['bsr'][:,lv]*s; ylo=z['yl'][:,lv]*s
    if flip: lo,hi=-z['ci'][:,lv,1],-z['ci'][:,lv,0]
    else:    lo,hi= z['ci'][:,lv,0], z['ci'][:,lv,1]
    ax=fig.add_subplot(gs[row,0]); rr,pp=pearsonr(xs,ys)
    ax.scatter(xs,ys,s=28,facecolor=RED,edgecolor='white',lw=.4,alpha=.75)
    m,b=np.polyfit(xs,ys,1); xx=np.linspace(xs.min(),xs.max(),50); ax.plot(xx,m*xx+b,color='0.15',lw=1.7)
    ax.text(.04,.955,f"r = {rr:.2f}, {'p < 0.001' if pp<.001 else f'p = {pp:.3f}'}",transform=ax.transAxes,
            fontsize=9.5,va='top',bbox=dict(boxstyle='round,pad=0.3',fc='white',ec='0.75',lw=.8))
    ax.set_xlabel(scatter_x); ax.set_ylabel('Phenomenology latent score'); ax.set_title(title,fontsize=10.5)
    for k,sp in ax.spines.items(): sp.set_visible(k in ('left','bottom'))
    ax=fig.add_subplot(gs[row,1]); o=np.argsort(bsr)
    ax.barh(range(len(bsr)),bsr[o],height=.75,
        color=[RED if v>1.96 else BLUE if v<-1.96 else NS for v in bsr[o]],
        edgecolor=[ERED if v>1.96 else EBLUE if v<-1.96 else ENS for v in bsr[o]],lw=.6)
    ax.set_yticks(range(len(bsr))); ax.set_yticklabels(np.array(xl)[o],fontsize=8.5)
    ax.axvline(1.96,ls='--',c='0.4',lw=.9); ax.axvline(-1.96,ls='--',c='0.4',lw=.9); ax.axvline(0,c='k',lw=.8)
    ax.set_xlabel('Bootstrap ratio'); ax.set_title('EEG metrics',fontsize=10.5)
    for k,sp in ax.spines.items(): sp.set_visible(k in ('left','bottom'))
    ax=fig.add_subplot(gs[row,2]); o=np.argsort(ylo); rel=(lo*hi>0)&sig
    ax.barh(range(len(ylo)),ylo[o],xerr=[ylo[o]-lo[o],hi[o]-ylo[o]],height=.7,
        color=[(RED if ylo[i]>0 else BLUE) if rel[i] else NS for i in o],
        edgecolor=[(ERED if ylo[i]>0 else EBLUE) if rel[i] else ENS for i in o],lw=.6,
        error_kw=dict(ecolor='black',lw=1.2,capsize=3,capthick=1.2))
    ax.set_yticks(range(len(ylo))); ax.set_yticklabels(np.array(yl_)[o],fontsize=9)
    ax.axvline(0,c='k',lw=.9); ax.set_xlabel('Loading (r with brain LV)'); ax.set_title('Phenomenology',fontsize=10.5)
    for k,sp in ax.spines.items(): sp.set_visible(k in ('left','bottom'))

zb=np.load(INTERIM/'pyls_between.npz'); zw=np.load(INTERIM/'pyls_within.npz'); zq=np.load(INTERIM/'pyls_quality.npz')

# ---------------- main figure ----------------
fig=plt.figure(figsize=(15,10.6))
gs=fig.add_gridspec(2,3,width_ratios=[1,1.12,1.05],wspace=0.47,hspace=0.42,top=0.925,bottom=0.075,left=0.075,right=0.985)
panelset(fig,gs,0,zb,0,FLIP['between'],XLAB,YLAB,f"Between-participant LV1: {zb['ve'][0]*100:.1f}% of covariance, $p$ = {zb['p'][0]:.4f}")
panelset(fig,gs,1,zw,0,FLIP['within'],XLAB,YLAB,f"Within-participant LV1: {zw['ve'][0]*100:.1f}% of covariance, $p$ = {zw['p'][0]:.4f}",
         scatter_x='Brain latent score (within-participant)')
for x,l in zip([0.012,0.335,0.672],'ABC'): fig.text(x,0.995,l,fontsize=17,fontweight='bold',va='top')
for x,l in zip([0.012,0.335,0.672],'DEF'): fig.text(x,0.508,l,fontsize=17,fontweight='bold',va='top')
fig.savefig(FIGS/'Figure_PLSC_main.png',dpi=400,facecolor='white'); fig.savefig(FIGS/'Figure_PLSC_main.pdf',facecolor='white')
plt.close(fig)

# ---------------- supplementary figure ----------------
u=pd.read_csv(TABLES/'univariate_factor_results.csv')
flabs=[f[0].split('\n')[0] for f in FACTORS]; mlabs=[m[0] for m in EEG]
R=u.pivot(index='brain',columns='factor',values='rho').loc[mlabs,flabs].to_numpy()
Qv=u.pivot(index='brain',columns='factor',values='q_fdr').loc[mlabs,flabs].to_numpy()
star=lambda q:'***' if q<.001 else '**' if q<.01 else '*' if q<.05 else ''
fig,axes=plt.subplots(1,2,figsize=(13.5,7.2),gridspec_kw={'width_ratios':[1.35,1],'wspace':0.6})
ax=axes[0]; im=ax.imshow(R,cmap=CMAP,vmin=-.6,vmax=.6,aspect='auto')
ax.set_xticks(range(len(FACTORS))); ax.set_xticklabels([f[0] for f in FACTORS],fontsize=9,rotation=90,va='top')
ax.set_yticks(range(len(EEG))); ax.set_yticklabels(mlabs,fontsize=9)
for i in range(len(EEG)):
    for j in range(len(FACTORS)):
        s=star(Qv[i,j]); ax.text(j,i,f"{s}\n{R[i,j]:.2f}" if s else f"{R[i,j]:.2f}",ha='center',va='center',
                                 fontsize=8,linespacing=1.0,color='white' if abs(R[i,j])>.42 else '0.1')
ax.axhline(5.5,color='k',lw=1.2); ax.axhline(10.5,color='k',lw=1.2)
ax.set_xticks(np.arange(-.5,len(FACTORS),1),minor=True); ax.set_yticks(np.arange(-.5,len(EEG),1),minor=True)
ax.grid(which='minor',color='white',lw=1.4); ax.tick_params(which='minor',length=0)
for sp in ax.spines.values(): sp.set_visible(False)
ax.set_title('ACAM-J–specific factors',fontsize=11)
cb=fig.colorbar(im,ax=ax,fraction=0.045,pad=0.04); cb.set_label('Spearman rho',fontsize=10); cb.outline.set_visible(False)
ax=axes[1]; ylq=zq['yl'][:,0]; lo,hi=zq['ci'][:,0,0],zq['ci'][:,0,1]; o=np.argsort(ylq)
kai=100/len(zq['ve']); sigq=(zq['p'][0]<.05) and (zq['ve'][0]*100>kai)
ax.barh(range(len(ylq)),ylq[o],xerr=[ylq[o]-lo[o],hi[o]-ylq[o]],height=.6,color=NS,edgecolor=ENS,lw=.6,
        error_kw=dict(ecolor='black',lw=1.2,capsize=3,capthick=1.2))
ax.set_yticks(range(len(ylq))); ax.set_yticklabels(np.array(QLAB)[o],fontsize=10)
ax.axvline(0,c='k',lw=.9); ax.set_xlabel('Loading (r with brain LV)')
ax.set_title(f"Run-level quality (brain averaged)\nLV1 = {zq['ve'][0]*100:.1f}% covariance, $p$ = {zq['p'][0]:.2f} (n.s.)",fontsize=10.5)
ax.text(.5,-0.14,'LV not significant — loadings not interpreted',transform=ax.transAxes,ha='center',fontsize=9,style='italic',color='0.35')
for k,sp in ax.spines.items(): sp.set_visible(k in ('left','bottom'))
for x,l in zip([0.012,0.545],'AB'): fig.text(x,0.985,l,fontsize=17,fontweight='bold',va='top')
fig.savefig(FIGS/'Figure_S_factors_and_quality.png',dpi=400,bbox_inches='tight',facecolor='white')
fig.savefig(FIGS/'Figure_S_factors_and_quality.pdf',bbox_inches='tight',facecolor='white')
plt.close(fig)
print(f"[05] figures written to {FIGS}")
