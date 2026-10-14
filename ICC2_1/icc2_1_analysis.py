#!/usr/bin/env python3
"""ICC(2,1) epoch-stability analysis (two-way random effects, absolute agreement).
Repeats epoch_stability_icc.py with a two-way model. Epoch position (temporal
order in the stored trial dimension) is the occasion/column. For each measure x
condition, a balanced matrix of the first k epochs is used (k = smallest epoch
count across retained subjects, minimum 3; at least 8 subjects). ICC(1,1) and
ICC(2,1) are computed on the same matrix and projected to K epochs with the
Spearman-Brown formula. Values are the global (128-channel mean) value per epoch.
Bootstrap CIs: 300 resamples over subjects."""
import os, numpy as np
# Paths below are from the environment where this was run; edit before use.
os.environ.setdefault('ICC_DATA_DIR','path/to/data')
os.environ['ICC_OUT_DIR']='path/to/icc2_outputs'
import importlib.util
spec=importlib.util.spec_from_file_location("esi","path/to/epoch_stability_icc.py")
esi=importlib.util.module_from_spec(spec); spec.loader.exec_module(esi)
from math import ceil
rng=np.random.default_rng(20260708)
OUT='path/to/icc2_outputs'; TAB=OUT+'/tables'; FIG=OUT+'/figures'
for d in (OUT,TAB,FIG): os.makedirs(d,exist_ok=True)
MEAS=esi.MEASURE_ORDER; CONDS=esi.CONDITIONS; CL=esi.COND_LABEL; M=esi.MEASURES
sb=esi.spearman_brown; mink=esi.min_k_for
KGRID=[1,2,3,4,5,6,8,10,12,15,18,20,25,30]

def balanced(mk,cond):
    d=esi.load_measure_condition(mk,cond)
    vecs=[]
    for s,arr in d.items():
        v=np.nanmean(arr,axis=0); v=v[np.isfinite(v)]
        if len(v)>=1: vecs.append(v)
    if len(vecs)<8: return None
    k=min(len(v) for v in vecs)
    if k<3: 
        vecs=[v for v in vecs if len(v)>=3]; 
        if len(vecs)<8: return None
        k=3
    return np.array([v[:k] for v in vecs])   # n x k (first k = temporal order)

def iccs(Mx):
    n,k=Mx.shape; grand=Mx.mean(); rowm=Mx.mean(1); colm=Mx.mean(0)
    MSR=k*np.sum((rowm-grand)**2)/(n-1)
    MSC=n*np.sum((colm-grand)**2)/(k-1)
    MSE=np.sum((Mx-(rowm[:,None]+colm[None,:]-grand))**2)/((n-1)*(k-1))
    MSW=np.sum((Mx-rowm[:,None])**2)/(n*(k-1))
    icc11=(MSR-MSW)/(MSR+(k-1)*MSW)
    icc21=(MSR-MSE)/(MSR+(k-1)*MSE+(k/n)*(MSC-MSE))
    return icc11,icc21,n,k,MSC/MSE if MSE>0 else np.nan

def boot_ci(Mx,K,nb=300):
    # 95% CI of R at K, bootstrap over subjects
    n=Mx.shape[0]; out=[]
    for _ in range(nb):
        idx=rng.integers(0,n,n); _,i21,_,_,_=iccs(Mx[idx])
        out.append(sb(i21,K))
    return np.nanpercentile(out,[2.5,97.5])

def main():
    import pandas as pd
    grows=[]; crows=[]
    for mk in MEAS:
        cfg=M[mk]; K=cfg['threshold']
        for cond in CONDS:
            Mx=balanced(mk,cond)
            if Mx is None: continue
            i11,i21,n,k,fcr=iccs(Mx)
            R21=float(sb(i21,K)); R11=float(sb(i11,K))
            lo,hi=boot_ci(Mx,K)
            grows.append(dict(measure=mk,kind=cfg['kind'],short=cfg['short'],condition=cond,
                condition_label=CL[cond],threshold_K=K,n_subjects=n,k_balanced=k,occ_F=fcr,
                ICC11=i11,ICC21=i21,R_ICC11_at_thr=R11,R_ICC21_at_thr=R21,
                R_ICC21_lo=float(lo),R_ICC21_hi=float(hi),
                minK_R75=mink(i21,0.75),minK_R90=mink(i21,0.90)))
            for kk in KGRID:
                crows.append(dict(measure=mk,kind=cfg['kind'],short=cfg['short'],condition=cond,
                    condition_label=CL[cond],K=kk,R_K=float(sb(i21,kk)),threshold_K=K))
    gdf=pd.DataFrame(grows); cdf=pd.DataFrame(crows)
    gdf.to_csv(TAB+'/icc21_global_summary.csv',index=False)
    cdf.to_csv(TAB+'/icc21_reliability_curves.csv',index=False)
    print("cells:",len(gdf),"| max |ICC21-ICC11| =",np.nanmax(np.abs(gdf.ICC21-gdf.ICC11)).round(4),
          "| mean |diff| =",np.nanmean(np.abs(gdf.ICC21-gdf.ICC11)).round(4))
    return gdf,cdf

if __name__=='__main__':
    main()
