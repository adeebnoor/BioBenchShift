#!/usr/bin/env python3
from __future__ import annotations

import argparse, itertools, json
from pathlib import Path

import numpy as np
import pandas as pd

import run_biosnap_dti_gate2_models as dti
import run_hetionet_gate2_models as het
import run_lightgcn_gate2 as lgcn
import run_h5_hypothesis_turnover as h5

KS=(100,500,1000)


def top_sets(score,mask):
    return {k:set(h5.top_flat(score,mask,k).tolist()) for k in KS}


def turnover(a,b,k):
    return 1-len(a[k]&b[k])/k


def within_stats(top_by_init):
    rows=[]
    for i,j in itertools.combinations(sorted(top_by_init),2):
        for k in KS:
            rows.append({'init_a':i,'init_b':j,'k':k,'turnover':turnover(top_by_init[i],top_by_init[j],k)})
    return pd.DataFrame(rows)


def cross_stats(A,B):
    rows=[]
    for i in sorted(set(A)&set(B)):
        for k in KS:
            rows.append({'init':i,'k':k,'turnover':turnover(A[i],B[i],k)})
    return pd.DataFrame(rows)


def fixed_split(edges,seed=0):
    train,test=h5.split_edges(edges,seed)
    return train,test


def run_dti(path,ninit):
    edges=dti.parse(path); left=sorted({a for a,_ in edges});right=sorted({b for _,b in edges});blocked=set(edges)
    train,_=fixed_split(edges,0);mask=h5.mask_candidates(edges,train,left,right)
    nmf={};svd={}
    for init in range(ninit):
        m=dti.fit_logistic_mf(train,left,right,blocked,np.random.default_rng(10000+init))
        nmf[init]=top_sets(h5.mf_matrix(m),mask)
        s=dti.fit_svd(train,left,right,32,20000+init)
        svd[init]=top_sets(h5.svd_matrix(s),mask)
    return mask.sum(), {'NeuralMF':nmf,'SVD':svd}, 'NeuralMF','SVD'


def run_dag(path,ninit):
    edges=het.parse_edges(path,'DaG','Disease::','Gene::');left=sorted({a for a,_ in edges});right=sorted({b for _,b in edges});blocked=set(edges)
    train,_=fixed_split(edges,0);mask=h5.mask_candidates(edges,train,left,right)
    nmf={};gcn={}
    for init in range(ninit):
        m=het.fit_logistic_mf(train,left,right,blocked,np.random.default_rng(30000+init))
        nmf[init]=top_sets(h5.mf_matrix(m),mask)
        g=lgcn.fit_lightgcn(train,left,right,blocked,40000+init)
        gcn[init]=top_sets(h5.lgcn_matrix(g,left,right),mask)
    return mask.sum(), {'NeuralMF':nmf,'LightGCN':gcn}, 'NeuralMF','LightGCN'


def summarize_family(family,candidate_count,models,random_winner,neutral_winner):
    rows=[]; summary={'candidate_count':int(candidate_count),'random_winner':random_winner,'neutral_winner':neutral_winner,'models':{}}
    for name,tops in models.items():
        w=within_stats(tops);w['family']=family;w['comparison']='within';w['model_a']=name;w['model_b']=name;rows.append(w)
        z={}
        for k in KS:
            q=w[w.k==k].turnover
            z[f'within_ht{k}_mean']=float(q.mean());z[f'within_ht{k}_max']=float(q.max());z[f'within_ht{k}_min']=float(q.min())
        summary['models'][name]=z
    c=cross_stats(models[random_winner],models[neutral_winner]);c['family']=family;c['comparison']='cross_selected';c['model_a']=random_winner;c['model_b']=neutral_winner;rows.append(c)
    cz={}
    for k in KS:
        q=c[c.k==k].turnover;cz[f'cross_ht{k}_mean']=float(q.mean());cz[f'cross_ht{k}_min']=float(q.min());cz[f'cross_ht{k}_max']=float(q.max())
        within_max=max(summary['models'][m][f'within_ht{k}_max'] for m in [random_winner,neutral_winner])
        cz[f'cross_minus_worst_within_ht{k}']=float(q.mean()-within_max)
    summary['cross_selected']=cz
    return pd.concat(rows,ignore_index=True),summary


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dti',required=True);ap.add_argument('--hetionet',required=True);ap.add_argument('--out-prefix',required=True);ap.add_argument('--inits',type=int,default=5);args=ap.parse_args()
    dti_n,dti_models,dti_r,dti_nw=run_dti(args.dti,args.inits)
    dag_n,dag_models,dag_r,dag_nw=run_dag(args.hetionet,args.inits)
    a,sa=summarize_family('DTI',dti_n,dti_models,dti_r,dti_nw);b,sb=summarize_family('Disease-Gene',dag_n,dag_models,dag_r,dag_nw)
    out=pd.concat([a,b],ignore_index=True);p=Path(args.out_prefix);p.parent.mkdir(parents=True,exist_ok=True);out.to_csv(str(p)+'_pairs.csv',index=False)
    S={'DTI':sa,'Disease-Gene':sb,'fixed_split_seed':0,'initializations':args.inits};Path(str(p)+'_summary.json').write_text(json.dumps(S,indent=2)+'\n')
    lines=['# H5 stability control — model-selection turnover versus initialization instability','','All comparisons use one fixed training split and one fixed candidate universe per relation family. Only model initialization/training randomness varies.','','| Family | Comparison | HT@100 | HT@500 | HT@1000 |','|---|---|---:|---:|---:|']
    for fam,s in [('DTI',sa),('Disease-Gene',sb)]:
        for m,z in s['models'].items(): lines.append(f"| {fam} | within {m} (mean pairwise) | {z['within_ht100_mean']:.3f} | {z['within_ht500_mean']:.3f} | {z['within_ht1000_mean']:.3f} |")
        z=s['cross_selected'];lines.append(f"| {fam} | conventional winner vs neutralized winner (mean paired-init) | {z['cross_ht100_mean']:.3f} | {z['cross_ht500_mean']:.3f} | {z['cross_ht1000_mean']:.3f} |")
    lines += ['','A benchmark-selection effect is strongest when cross-selected-model turnover materially exceeds ordinary within-model initialization turnover on the same candidate universe. No biological correctness is inferred from stability alone.']
    md='\n'.join(lines)+'\n';Path(str(p)+'.md').write_text(md);print(md)

if __name__=='__main__': main()
