#!/usr/bin/env python3
"""Build Table S14 source data: per-therapeutic-area bootstrap distributions.

Uses only the frozen Open Targets target vectors produced by
run_open_targets_disease_area_expansion.py. No model fitting, mapping, target
selection, or therapeutic-area definition is changed here.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

KS=(100,500,1000)
MODELS=("SVD","NeuralMF")
NBOOT=2000
BASESEED=20260913


def matrix_for(df, model, k, areas):
    z=df[(df.model==model)&(df.k==k)].copy()
    genes=sorted(z.entrez_gene.astype(str).unique(), key=lambda x:(len(x),x))
    gi={g:i for i,g in enumerate(genes)}; ai={a:i for i,a in enumerate(areas)}
    X=np.zeros((len(genes),len(areas)),float)
    for _,r in z.iterrows():
        X[gi[str(r.entrez_gene)],ai[str(r.therapeutic_area_id)]]=float(r.target_normalized_area_weight)
    return genes,X


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--target-vectors',required=True)
    ap.add_argument('--profiles',required=True)
    ap.add_argument('--out-prefix',required=True)
    a=ap.parse_args()
    tv=pd.read_csv(a.target_vectors,dtype={'entrez_gene':str,'therapeutic_area_id':str})
    prof=pd.read_csv(a.profiles,dtype={'therapeutic_area_id':str})
    areas=sorted(prof.therapeutic_area_id.astype(str).unique())
    names=(prof[['therapeutic_area_id','therapeutic_area_name']]
           .drop_duplicates().set_index('therapeutic_area_id').therapeutic_area_name.to_dict())
    summary=[]; draws=[]
    for k in KS:
        mats={}
        for m in MODELS:
            genes,X=matrix_for(tv,m,k,areas)
            mats[m]=(genes,X)
        observed={m:mats[m][1].mean(axis=0) for m in MODELS}
        rng=np.random.default_rng(BASESEED+k)
        boot={m:np.empty((NBOOT,len(areas)),float) for m in MODELS}
        for b in range(NBOOT):
            for m in MODELS:
                X=mats[m][1]
                idx=rng.integers(0,len(X),size=len(X))
                boot[m][b]=X[idx].mean(axis=0)
        diff=boot['SVD']-boot['NeuralMF']
        for j,area in enumerate(areas):
            s=boot['SVD'][:,j]; n=boot['NeuralMF'][:,j]; d=diff[:,j]
            q=lambda x:np.percentile(x,[2.5,50,97.5])
            qs,qn,qd=q(s),q(n),q(d)
            summary.append({
              'k':k,'therapeutic_area_id':area,'therapeutic_area_name':names.get(area,''),
              'observed_svd_mass':observed['SVD'][j],
              'observed_neuralmf_mass':observed['NeuralMF'][j],
              'observed_difference_svd_minus_neuralmf':observed['SVD'][j]-observed['NeuralMF'][j],
              'svd_bootstrap_median':qs[1],'svd_bootstrap_95ci_low':qs[0],'svd_bootstrap_95ci_high':qs[2],
              'neuralmf_bootstrap_median':qn[1],'neuralmf_bootstrap_95ci_low':qn[0],'neuralmf_bootstrap_95ci_high':qn[2],
              'difference_bootstrap_median':qd[1],'difference_bootstrap_95ci_low':qd[0],'difference_bootstrap_95ci_high':qd[2],
              'bootstrap_replicates':NBOOT
            })
            for b in range(NBOOT):
                draws.append({'k':k,'bootstrap':b,'therapeutic_area_id':area,
                              'svd_mass':s[b],'neuralmf_mass':n[b],
                              'difference_svd_minus_neuralmf':d[b]})
    p=Path(a.out_prefix); p.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(summary).to_csv(str(p)+'_summary.csv',index=False)
    pd.DataFrame(draws).to_csv(str(p)+'_draws.csv',index=False)
    print(f'areas={len(areas)} rows_summary={len(summary)} rows_draws={len(draws)}')

if __name__=='__main__': main()
