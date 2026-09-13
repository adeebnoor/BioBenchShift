#!/usr/bin/env python3
from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

KS=(100,500,1000)
MODELS=("SVD","NeuralMF")


def load_vectors(path):
    df=pd.read_csv(path)
    vectors={}
    names={}
    for (m,k,g),q in df.groupby(["model","k","entrez_gene"], sort=False):
        v={str(r.therapeutic_area_id):float(r.target_normalized_area_weight) for _,r in q.iterrows()}
        vectors[(str(m),int(k),str(g))]=v
        for _,r in q.iterrows(): names[str(r.therapeutic_area_id)]=str(r.therapeutic_area_name)
    genes={(m,k):sorted({str(g) for mm,kk,g in vectors if mm==m and kk==k}) for m in MODELS for k in KS}
    areas=sorted(names)
    return vectors,genes,areas,names


def mean_profile(vs,areas):
    z=defaultdict(float)
    for v in vs:
        for a,w in v.items(): z[a]+=w
    n=len(vs)
    return {a:(z[a]/n if n else float('nan')) for a in areas}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--target-vectors',required=True); ap.add_argument('--out-prefix',required=True); ap.add_argument('--n',type=int,default=2000); a=ap.parse_args()
    vectors,genes,areas,names=load_vectors(a.target_vectors)
    summary=[]; reps=[]
    for k in KS:
        av=[vectors[("SVD",k,g)] for g in genes[("SVD",k)]]
        bv=[vectors[("NeuralMF",k,g)] for g in genes[("NeuralMF",k)]]
        obs_a=mean_profile(av,areas); obs_b=mean_profile(bv,areas)
        rng=np.random.default_rng(20260913+k)
        vals={ar:[] for ar in areas}
        for b in range(a.n):
            aa=[av[int(i)] for i in rng.integers(len(av),size=len(av))]
            bb=[bv[int(i)] for i in rng.integers(len(bv),size=len(bv))]
            pa=mean_profile(aa,areas); pb=mean_profile(bb,areas)
            for ar in areas:
                d=pa[ar]-pb[ar]; vals[ar].append(d)
                reps.append({'k':k,'replicate':b,'therapeutic_area_id':ar,'difference_svd_minus_neuralmf':d})
        for ar in areas:
            x=np.asarray(vals[ar],float)
            q=np.quantile(x,[.025,.25,.5,.75,.975])
            summary.append({'k':k,'therapeutic_area_id':ar,'therapeutic_area_name':names.get(ar,''),
                'svd_observed_mass':obs_a[ar],'neuralmf_observed_mass':obs_b[ar],
                'observed_difference_svd_minus_neuralmf':obs_a[ar]-obs_b[ar],
                'bootstrap_q025':q[0],'bootstrap_q25':q[1],'bootstrap_median':q[2],'bootstrap_q75':q[3],'bootstrap_q975':q[4],
                'bootstrap_p_diff_gt_0':float(np.mean(x>0)),'bootstrap_replicates':a.n})
    p=Path(a.out_prefix); p.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(summary).to_csv(str(p)+'_summary.csv',index=False)
    pd.DataFrame(reps).to_csv(str(p)+'_replicates.csv.gz',index=False,compression='gzip')

if __name__=='__main__': main()
