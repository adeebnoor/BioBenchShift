#!/usr/bin/env python3
from __future__ import annotations

import argparse, json
from pathlib import Path

import numpy as np
from scipy import sparse

import run_biogrid_temporal_h4 as h4

KS=(100,500,1000)
MAXK=max(KS)


def score_block_mf(model, rows):
    ix,U,b=model
    return U[rows] @ U.T + b[rows,None] + b[None,:]


def score_block_svd(model, rows):
    ix,U,V=model
    return 0.5*(U[rows] @ V.T + V[rows] @ U.T)


def keep_top(state_ids,state_scores,new_ids,new_scores,k=MAXK):
    if state_ids is None:
        ids=new_ids; sc=new_scores
    else:
        ids=np.concatenate([state_ids,new_ids]); sc=np.concatenate([state_scores,new_scores])
    if len(sc)<=k:
        o=np.argsort(sc)[::-1]; return ids[o],sc[o]
    sel=np.argpartition(sc,-k)[-k:]; o=sel[np.argsort(sc[sel])[::-1]]
    return ids[o],sc[o]


def ensemble_rank(models, scorer, nodes, known_edges, block=256):
    n=len(nodes);ix={x:i for i,x in enumerate(nodes)}
    rr=[];cc=[]
    for a,b in known_edges:
        i,j=ix[a],ix[b];rr += [i,j];cc += [j,i]
    A=sparse.csr_matrix((np.ones(len(rr),dtype=np.int8),(rr,cc)),shape=(n,n))
    variants=['full']+[f'loo{i}' for i in range(len(models))]
    tops={v:(None,None) for v in variants}
    for start in range(0,n,block):
        stop=min(n,start+block); rows=np.arange(start,stop)
        mats=[np.asarray(scorer(m,rows),dtype=np.float32) for m in models]
        M=np.stack(mats,axis=0)
        vmat={'full':M.mean(axis=0)}
        for i in range(len(models)):
            idx=[j for j in range(len(models)) if j!=i]
            vmat[f'loo{i}']=M[idx].mean(axis=0)
        for name,S in vmat.items():
            S=S.copy()
            for r,i in enumerate(rows):
                S[r,:i+1]=-np.inf
                js=A.indices[A.indptr[i]:A.indptr[i+1]]
                if len(js): S[r,js]=-np.inf
            flat=S.ravel(); valid=np.isfinite(flat); ids=np.flatnonzero(valid)
            if len(ids)==0: continue
            vals=flat[ids]; take=min(MAXK,len(vals)); loc=np.argpartition(vals,-take)[-take:]
            local=ids[loc]; rlocal=local//n; jlocal=local%n
            global_ids=rows[rlocal].astype(np.int64)*n+jlocal.astype(np.int64)
            oldi,olds=tops[name];tops[name]=keep_top(oldi,olds,global_ids,vals[loc])
    return {v:set(tops[v][0][:MAXK].tolist()) for v in variants}


def subset(topset,k):
    # sets lose ordering, so rank output must retain ordered ids; replaced below by ordered arrays.
    raise RuntimeError


def ordered_rank(models,scorer,nodes,known_edges,block=256):
    n=len(nodes);ix={x:i for i,x in enumerate(nodes)}
    rr=[];cc=[]
    for a,b in known_edges:
        i,j=ix[a],ix[b];rr += [i,j];cc += [j,i]
    A=sparse.csr_matrix((np.ones(len(rr),dtype=np.int8),(rr,cc)),shape=(n,n))
    variants=['full']+[f'loo{i}' for i in range(len(models))]
    tops={v:(None,None) for v in variants}
    for start in range(0,n,block):
        stop=min(n,start+block);rows=np.arange(start,stop)
        M=np.stack([np.asarray(scorer(m,rows),dtype=np.float32) for m in models],axis=0)
        V={'full':M.mean(axis=0)}
        for i in range(len(models)):
            idx=[j for j in range(len(models)) if j!=i];V[f'loo{i}']=M[idx].mean(axis=0)
        for name,S0 in V.items():
            S=S0.copy()
            for r,i in enumerate(rows):
                S[r,:i+1]=-np.inf
                js=A.indices[A.indptr[i]:A.indptr[i+1]]
                if len(js):S[r,js]=-np.inf
            flat=S.ravel();ids=np.flatnonzero(np.isfinite(flat));
            if len(ids)==0:continue
            vals=flat[ids];take=min(MAXK,len(vals));loc=np.argpartition(vals,-take)[-take:];local=ids[loc]
            global_ids=rows[local//n].astype(np.int64)*n+(local%n).astype(np.int64)
            oi,os=tops[name];tops[name]=keep_top(oi,os,global_ids,vals[loc])
    return {v:tops[v][0] for v in variants}


def HT(a,b,k): return 1-len(set(a[:k].tolist())&set(b[:k].tolist()))/k


def loo_stats(ranks):
    loo=[ranks[k] for k in sorted(ranks) if k.startswith('loo')]
    out={}
    for k in KS:
        vals=[]
        for i in range(len(loo)):
            for j in range(i+1,len(loo)): vals.append(HT(loo[i],loo[j],k))
        out[str(k)]={'mean':float(np.mean(vals)),'max':float(np.max(vals)),'min':float(np.min(vals))}
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--old',required=True);ap.add_argument('--out-prefix',required=True);ap.add_argument('--inits',type=int,default=3);args=ap.parse_args()
    edges,member,_,_=h4.load_human_edges(args.old);nodes=sorted({x for e in edges for x in e});blocked=set(edges)
    mf=[];svd=[]
    for seed in range(args.inits):
        mf.append(h4.fit_mf(edges,nodes,blocked,71000+seed))
        svd.append(h4.fit_svd(edges,nodes,72000+seed))
    mr=ordered_rank(mf,score_block_mf,nodes,edges);sr=ordered_rank(svd,score_block_svd,nodes,edges)
    ms=loo_stats(mr);ss=loo_stats(sr);cross={str(k):HT(mr['full'],sr['full'],k) for k in KS}
    candidate_count=len(nodes)*(len(nodes)-1)//2-len(edges)
    S={'release':'5.0.250','historical_edges':len(edges),'historical_nodes':len(nodes),'candidate_non_edges':candidate_count,'initializations':args.inits,'conventional_winner':'NeuralMF','neutralized_winner':'SVD','cross_ensemble_ht':cross,'neuralmf_loo':ms,'svd_loo':ss}
    p=Path(args.out_prefix);p.parent.mkdir(parents=True,exist_ok=True);Path(str(p)+'_summary.json').write_text(json.dumps(S,indent=2)+'\n')
    lines=['# PPI H5 — BioGRID hypothesis identity with ensemble stability','','Both selected model families are trained on the complete frozen BioGRID 5.0.250 human MV-Physical graph. Each family is a three-initialization score ensemble. Leave-one-initialization-out ensembles quantify within-family ranking instability on the identical full historical non-edge universe.','','| K | NeuralMF vs SVD ensemble HT | NeuralMF LOO HT mean | SVD LOO HT mean |','|---:|---:|---:|---:|']
    for k in KS: lines.append(f"| {k:,} | {cross[str(k)]:.3f} | {ms[str(k)]['mean']:.3f} | {ss[str(k)]['mean']:.3f} |")
    lines += ['',f"Candidate universe: **{candidate_count:,}** historical non-edges among **{len(nodes):,}** historical proteins.",'','This experiment measures scientific-decision identity only. Later-evidence correctness is evaluated separately in H4/H4b.']
    md='\n'.join(lines)+'\n';Path(str(p)+'.md').write_text(md);print(md)

if __name__=='__main__': main()
