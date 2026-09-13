#!/usr/bin/env python3
from __future__ import annotations

import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import hypergeom

import run_biogrid_temporal_h4 as h4
import run_biogrid_future_yield as h4b

KS=(100,500,1000,5000,10000,50000)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--old',required=True);ap.add_argument('--new',required=True);ap.add_argument('--out-prefix',required=True);ap.add_argument('--boot',type=int,default=20000);args=ap.parse_args()
    old,_,_,_=h4.load_human_edges(args.old);new,_,_,_=h4.load_human_edges(args.new)
    oldset=set(old);newset=set(new);nodes=sorted({x for e in old for x in e});nodeset=set(nodes);ix={x:i for i,x in enumerate(nodes)};n=len(nodes)
    future=sorted(e for e in newset-oldset if e[0] in nodeset and e[1] in nodeset)
    future_flat=np.array([min(ix[a],ix[b])*n+max(ix[a],ix[b]) for a,b in future],dtype=np.int64)
    N=n*(n-1)//2-len(oldset);M=len(future)
    models=h4.fit_all(old,nodes,oldset,100)
    tops={}
    for name,(model,_) in models.items():
        ids,_=h4b.global_topk(name,model,nodes,oldset,50000);tops[name]=ids
    rng=np.random.default_rng(20260913);rows=[];diffrows=[]
    indicators={}
    for name,ids in tops.items():
        rank={int(q):i+1 for i,q in enumerate(ids)}
        ind={k:np.array([rank.get(int(q),10**18)<=k for q in future_flat],dtype=np.int8) for k in KS}
        indicators[name]=ind
        for k in KS:
            hits=int(ind[k].sum());p=float(hypergeom.sf(hits-1,N,M,k))
            # Bootstrap recall among later-supported edges.
            vals=[]
            for _ in range(args.boot):
                z=rng.integers(0,M,M);vals.append(float(ind[k][z].mean()))
            lo,hi=np.quantile(vals,[.025,.975])
            rows.append({'model':name,'k':k,'hits':hits,'recall':hits/M,'recall_boot_lo':float(lo),'recall_boot_hi':float(hi),'uniform_hypergeom_p':p})
    for k in KS:
        a=indicators['SVD'][k].astype(float);b=indicators['NeuralMF'][k].astype(float);d=a-b
        vals=[]
        for _ in range(args.boot):
            z=rng.integers(0,M,M);vals.append(float(d[z].mean()))
        lo,hi=np.quantile(vals,[.025,.975]);obs=float(d.mean());p_nonpos=(np.sum(np.asarray(vals)<=0)+1)/(args.boot+1)
        diffrows.append({'k':k,'svd_minus_neuralmf_recall':obs,'paired_boot_lo':float(lo),'paired_boot_hi':float(hi),'bootstrap_p_nonpositive':float(p_nonpos),'svd_hits':int(a.sum()),'neuralmf_hits':int(b.sum())})
    p=Path(args.out_prefix);p.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(str(p)+'_models.csv',index=False);pd.DataFrame(diffrows).to_csv(str(p)+'_paired.csv',index=False)
    summary={'bootstrap_replicates':args.boot,'bootstrap_seed':20260913,'later_supported_edges':M,'candidate_universe':N,'paired_contrast':'SVD minus NeuralMF recall among later-added BioGRID edges','rows':diffrows}
    Path(str(p)+'_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# H4c — uncertainty for non-circular BioGRID future yield','',f'Paired bootstrap over **{M:,}** later-added relations; **{args.boot:,}** resamples. All rankings are generated from the frozen historical network before later-edge labels are used.','','| K | SVD hits | NeuralMF hits | Recall difference | Paired bootstrap 95% CI | P(bootstrap difference <= 0) |','|---:|---:|---:|---:|---:|---:|']
    for r in diffrows:lines.append(f"| {r['k']:,} | {r['svd_hits']:,} | {r['neuralmf_hits']:,} | {r['svd_minus_neuralmf_recall']:.4f} | [{r['paired_boot_lo']:.4f}, {r['paired_boot_hi']:.4f}] | {r['bootstrap_p_nonpositive']:.4g} |")
    lines+=['','Hypergeometric enrichment probabilities against a uniformly random top-K set are stored in the model table. The paired bootstrap is the primary uncertainty summary for the prespecified SVD-versus-NeuralMF future-evidence contrast.']
    Path(str(p)+'.md').write_text('\n'.join(lines)+'\n');print('\n'.join(lines))

if __name__=='__main__':main()
