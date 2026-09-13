#!/usr/bin/env python3
from __future__ import annotations

import argparse, json
from pathlib import Path
import numpy as np
from scipy import sparse

import run_biogrid_temporal_h4 as h4


def score_block(model_name, model, rows):
    if model_name == 'SVD':
        ix,U,V=model
        return 0.5*(U[rows] @ V.T + V[rows] @ U.T)
    if model_name == 'NeuralMF':
        ix,U,b=model
        return U[rows] @ U.T + b[rows,None] + b[None,:]
    if model_name == 'LightGCN':
        ix,Z=model
        return Z[rows] @ Z.T
    raise KeyError(model_name)


def global_topk(model_name, model, nodes, known_edges, max_k=50000, block=256):
    ix={x:i for i,x in enumerate(nodes)}
    n=len(nodes)
    rr=[];cc=[]
    for a,b in known_edges:
        i,j=ix[a],ix[b]
        rr += [i,j]; cc += [j,i]
    A=sparse.csr_matrix((np.ones(len(rr),dtype=np.int8),(rr,cc)),shape=(n,n))
    cand_scores=[]; cand_flat=[]
    for start in range(0,n,block):
        stop=min(n,start+block); rows=np.arange(start,stop)
        S=np.asarray(score_block(model_name,model,rows),dtype=np.float64)
        # only strict upper triangle => every unordered pair exactly once
        for r,i in enumerate(rows):
            S[r,:i+1] = -np.inf
            js=A.indices[A.indptr[i]:A.indptr[i+1]]
            if len(js): S[r,js] = -np.inf
        flat=S.ravel(); valid=np.isfinite(flat); nv=int(valid.sum())
        if nv==0: continue
        take=min(max_k,nv)
        ids=np.flatnonzero(valid)
        vals=flat[ids]
        loc=np.argpartition(vals,-take)[-take:]
        local=ids[loc]
        cand_scores.append(flat[local])
        rlocal=local//n; jlocal=local%n
        iglobal=rows[rlocal]
        cand_flat.append(iglobal.astype(np.int64)*n+jlocal.astype(np.int64))
    sc=np.concatenate(cand_scores); ids=np.concatenate(cand_flat)
    take=min(max_k,len(sc)); sel=np.argpartition(sc,-take)[-take:]; order=sel[np.argsort(sc[sel])[::-1]]
    return ids[order],sc[order]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--old',required=True);ap.add_argument('--new',required=True);ap.add_argument('--out-prefix',required=True);ap.add_argument('--max-k',type=int,default=50000);args=ap.parse_args()
    old,old_member,_,_=h4.load_human_edges(args.old); new,new_member,_,_=h4.load_human_edges(args.new)
    oldset=set(old);newset=set(new);nodes=sorted({x for e in old for x in e});nodeset=set(nodes); ix={x:i for i,x in enumerate(nodes)};n=len(nodes)
    future=sorted(e for e in newset-oldset if e[0] in nodeset and e[1] in nodeset); fidx={tuple(sorted((ix[a],ix[b]))) for a,b in future}
    candidate_count=n*(n-1)//2-len(oldset); prevalence=len(future)/candidate_count
    models=h4.fit_all(old,nodes,oldset,100)
    ks=[100,500,1000,5000,10000,args.max_k]
    rows=[]
    for name,(model,fn) in models.items():
        ids,scores=global_topk(name,model,nodes,oldset,args.max_k)
        pairs=[(int(q//n),int(q%n)) for q in ids]
        for k in ks:
            kk=min(k,len(pairs)); hits=sum(p in fidx for p in pairs[:kk]); precision=hits/kk; recall=hits/len(future); expected=kk*prevalence
            rows.append({'model':name,'k':kk,'future_hits':hits,'precision_lower_bound':precision,'recall_of_later_added':recall,'chance_expected_hits':expected,'enrichment_over_uniform':hits/expected if expected else float('nan')})
    import pandas as pd
    out=pd.DataFrame(rows);p=Path(args.out_prefix);p.parent.mkdir(parents=True,exist_ok=True);out.to_csv(str(p)+'_table.csv',index=False)
    summary={'old_release':'5.0.250','new_release':'5.0.261','historical_nodes':n,'historical_edges':len(oldset),'candidate_non_edges':candidate_count,'later_added_closed_world_edges':len(future),'future_prevalence_in_candidate_universe':prevalence,'ranking_seed':100,'note':'Precision is a lower bound because pairs not added by 5.0.261 are persistent-unobserved, not proven negatives.'}
    Path(str(p)+'_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# H4b — non-circular prospective discovery yield','',f"Frozen 2025 candidate universe: **{candidate_count:,}** unknown human protein pairs among historical nodes; **{len(future):,}** became BioGRID MV-Physical edges by 5.0.261.",'','| Model | K | Later-added hits | Precision lower bound | Recall of later-added edges | Enrichment vs uniform |','|---|---:|---:|---:|---:|---:|']
    for _,r in out.iterrows(): lines.append(f"| {r.model} | {int(r.k):,} | {int(r.future_hits):,} | {r.precision_lower_bound:.4f} | {r.recall_of_later_added:.4f} | {r.enrichment_over_uniform:.1f}× |")
    lines += ['','No degree matching or future-control sampling is used in this ranking analysis. Every model ranks the same full historical non-edge universe. Later-added BioGRID edges are treated as future evidence; persistent-unobserved pairs are not called biological negatives.']
    md='\n'.join(lines)+'\n';Path(str(p)+'.md').write_text(md);print(md)

if __name__=='__main__': main()
