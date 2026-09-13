#!/usr/bin/env python3
from __future__ import annotations

import argparse, json
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


def ht(A,B,k): return 1-len(A[k]&B[k])/k


def ensemble_and_loo(score_list,mask):
    stack=np.stack(score_list,axis=0)
    full=top_sets(stack.mean(axis=0),mask)
    loo=[]
    for i in range(len(score_list)):
        idx=[j for j in range(len(score_list)) if j!=i]
        loo.append(top_sets(stack[idx].mean(axis=0),mask))
    return full,loo


def loo_instability(loo):
    rows=[]
    for i in range(len(loo)):
        for j in range(i+1,len(loo)):
            for k in KS: rows.append({'i':i,'j':j,'k':k,'turnover':ht(loo[i],loo[j],k)})
    return pd.DataFrame(rows)


def fixed_split(edges): return h5.split_edges(edges,0)


def run_dti(path,ninit):
    edges=dti.parse(path);left=sorted({a for a,_ in edges});right=sorted({b for _,b in edges});blocked=set(edges)
    train,_=fixed_split(edges);mask=h5.mask_candidates(edges,train,left,right)
    mf=[];svd=[]
    for init in range(ninit):
        m=dti.fit_logistic_mf(train,left,right,blocked,np.random.default_rng(51000+init));mf.append(h5.mf_matrix(m))
        s=dti.fit_svd(train,left,right,32,52000+init);svd.append(h5.svd_matrix(s))
    A,Aloo=ensemble_and_loo(mf,mask);B,Bloo=ensemble_and_loo(svd,mask)
    return int(mask.sum()),'NeuralMF','SVD',A,Aloo,B,Bloo


def run_dag(path,ninit):
    edges=het.parse_edges(path,'DaG','Disease::','Gene::');left=sorted({a for a,_ in edges});right=sorted({b for _,b in edges});blocked=set(edges)
    train,_=fixed_split(edges);mask=h5.mask_candidates(edges,train,left,right)
    mf=[];gcn=[]
    for init in range(ninit):
        m=het.fit_logistic_mf(train,left,right,blocked,np.random.default_rng(61000+init));mf.append(h5.mf_matrix(m))
        g=lgcn.fit_lightgcn(train,left,right,blocked,62000+init);gcn.append(h5.lgcn_matrix(g,left,right))
    A,Aloo=ensemble_and_loo(mf,mask);B,Bloo=ensemble_and_loo(gcn,mask)
    return int(mask.sum()),'NeuralMF','LightGCN',A,Aloo,B,Bloo


def summarize(family,result):
    n,a_name,b_name,A,Aloo,B,Bloo=result
    wa=loo_instability(Aloo);wb=loo_instability(Bloo)
    out={'candidate_count':n,'conventional_ensemble':a_name,'neutralized_ensemble':b_name,'initializations':len(Aloo),'metrics':{}}
    rows=[]
    for k in KS:
        cross=ht(A,B,k);qa=wa[wa.k==k].turnover;qb=wb[wb.k==k].turnover
        m={'cross_ht':float(cross),'conventional_loo_ht_mean':float(qa.mean()),'conventional_loo_ht_max':float(qa.max()),'neutralized_loo_ht_mean':float(qb.mean()),'neutralized_loo_ht_max':float(qb.max())}
        m['cross_minus_worst_loo_mean']=float(cross-max(qa.mean(),qb.mean()))
        m['cross_minus_worst_loo_max']=float(cross-max(qa.max(),qb.max()))
        out['metrics'][str(k)]=m
        rows.append({'family':family,'k':k,**m})
    return out,rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dti',required=True);ap.add_argument('--hetionet',required=True);ap.add_argument('--out-prefix',required=True);ap.add_argument('--inits',type=int,default=5);args=ap.parse_args()
    s1,r1=summarize('DTI',run_dti(args.dti,args.inits));s2,r2=summarize('Disease-Gene',run_dag(args.hetionet,args.inits));S={'DTI':s1,'Disease-Gene':s2,'fixed_split_seed':0}
    p=Path(args.out_prefix);p.parent.mkdir(parents=True,exist_ok=True);Path(str(p)+'_summary.json').write_text(json.dumps(S,indent=2)+'\n');pd.DataFrame(r1+r2).to_csv(str(p)+'_table.csv',index=False)
    lines=['# H5b confirmatory ensemble analysis','','Each model family is represented by the mean score of five independently initialized fits on one fixed training split. Ensemble stability is estimated by leave-one-initialization-out (LOO) ensembles on the same candidate universe.','','| Family | K | Cross-selected ensemble HT | Conventional LOO HT mean | Neutralized LOO HT mean | Cross − worst LOO mean |','|---|---:|---:|---:|---:|---:|']
    for fam,s in [('DTI',s1),('Disease-Gene',s2)]:
        for k in KS:
            m=s['metrics'][str(k)];lines.append(f"| {fam} | {k:,} | {m['cross_ht']:.3f} | {m['conventional_loo_ht_mean']:.3f} | {m['neutralized_loo_ht_mean']:.3f} | {m['cross_minus_worst_loo_mean']:.3f} |")
    lines += ['','Interpretation rule: a clean H5 replication requires cross-selected-model turnover to remain materially larger than ensemble instability. Ensembling is confirmatory; the original single-fit H5 and its negative/boundary findings remain reported.']
    md='\n'.join(lines)+'\n';Path(str(p)+'.md').write_text(md);print(md)

if __name__=='__main__': main()
