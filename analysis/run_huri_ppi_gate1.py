#!/usr/bin/env python3
from __future__ import annotations

import argparse, json, math
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd


def auc_mw(pos, neg):
    p=np.asarray(pos,float); n=np.asarray(neg,float)
    if not len(p) or not len(n): return float('nan')
    ranks=pd.Series(np.concatenate([p,n])).rank(method='average').to_numpy()
    u=ranks[:len(p)].sum()-len(p)*(len(p)+1)/2
    return float(u/(len(p)*len(n)))

def dbin(d):
    return 0 if d<=0 else int(math.floor(math.log2(d)))+1

def norm(a,b):
    a,b=str(a).strip(),str(b).strip()
    return tuple(sorted((a,b)))

def score(e,deg):
    return float(np.log1p(deg.get(e[0],0))*np.log1p(deg.get(e[1],0)))

def parse(path):
    # HuRI.tsv is a simple two-column TSV of Ensembl gene IDs.
    d=pd.read_csv(path,sep='\t',comment='#',header=None,dtype=str)
    if d.shape[1] < 2:
        raise ValueError(f'Expected >=2 TSV columns, got {d.shape}')
    edges=sorted({norm(a,b) for a,b in zip(d.iloc[:,0],d.iloc[:,1]) if pd.notna(a) and pd.notna(b) and str(a).strip()!=str(b).strip()})
    if len(edges)<1000: raise ValueError(f'Parsed only {len(edges)} edges')
    return edges

def sample_random_nonedges(rng,n,nodes,blocked):
    out=set(); attempts=0; limit=max(500000,n*100)
    while len(out)<n and attempts<limit:
        attempts+=1
        a=nodes[int(rng.integers(len(nodes)))]; b=nodes[int(rng.integers(len(nodes)))]
        if a==b: continue
        e=norm(a,b)
        if e in blocked or e in out: continue
        out.add(e)
    if len(out)<n: raise RuntimeError(f'Only sampled {len(out)}/{n} random nonedges')
    return list(out)

def targeted_match(rng,positives,nodes,blocked,deg):
    bybin=defaultdict(list)
    for v in nodes: bybin[dbin(deg.get(v,0))].append(v)
    used=set(); mp=[]; mn=[]
    for i in rng.permutation(len(positives)):
        p=positives[int(i)]
        b1,b2=sorted((dbin(deg.get(p[0],0)),dbin(deg.get(p[1],0))))
        A=bybin.get(b1,[]); B=bybin.get(b2,[])
        if not A or not B: continue
        found=None
        for _ in range(600):
            a=A[int(rng.integers(len(A)))]; b=B[int(rng.integers(len(B)))]
            if a==b: continue
            e=norm(a,b)
            if e in blocked or e in used: continue
            # Ensure unordered signature is exactly preserved.
            if sorted((dbin(deg.get(e[0],0)),dbin(deg.get(e[1],0)))) != [b1,b2]: continue
            found=e; break
        if found is not None:
            used.add(found); mp.append(p); mn.append(found)
    return mp,mn

def one(edges,seed,test_fraction=.2):
    rng=np.random.default_rng(seed); idx=rng.permutation(len(edges)); nt=max(1,int(round(len(edges)*test_fraction)))
    test=[edges[i] for i in idx[:nt]]; train=[edges[i] for i in idx[nt:]]; allset=set(edges)
    deg=Counter();
    for a,b in train: deg[a]+=1; deg[b]+=1
    nodes=sorted({x for e in edges for x in e})
    seen=[e for e in test if deg.get(e[0],0)>0 and deg.get(e[1],0)>0]
    if len(seen)<10: raise ValueError('Too few seen-endpoint test edges')
    rnd=sample_random_nonedges(rng,len(seen),nodes,allset)
    mp,mn=targeted_match(rng,seen,nodes,allset,deg)
    ar=auc_mw([score(e,deg) for e in seen],[score(e,deg) for e in rnd])
    am=auc_mw([score(e,deg) for e in mp],[score(e,deg) for e in mn]) if mp else float('nan')
    return {'seed':seed,'n_train':len(train),'n_test':len(test),'n_seen':len(seen),'seen_fraction':len(seen)/len(test),
            'n_matched':len(mp),'matched_fraction':len(mp)/len(seen),'auc_random':ar,'auc_degree_matched':am,'inflation_auc':ar-am}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input',required=True); ap.add_argument('--out-csv',required=True); ap.add_argument('--out-json',required=True); ap.add_argument('--out-md',required=True); ap.add_argument('--replicates',type=int,default=20); args=ap.parse_args()
    edges=parse(args.input); rows=pd.DataFrame([one(edges,s) for s in range(args.replicates)]); Path(args.out_csv).parent.mkdir(parents=True,exist_ok=True); rows.to_csv(args.out_csv,index=False)
    s={'benchmark':'HuRI human protein-protein interaction network','edge_count':len(edges),'nodes':len({x for e in edges for x in e}),'replicates':args.replicates,
       'auc_random_mean':float(rows.auc_random.mean()),'auc_random_sd':float(rows.auc_random.std(ddof=1)),'auc_degree_matched_mean':float(rows.auc_degree_matched.mean()),'auc_degree_matched_sd':float(rows.auc_degree_matched.std(ddof=1)),
       'inflation_auc_mean':float(rows.inflation_auc.mean()),'inflation_auc_sd':float(rows.inflation_auc.std(ddof=1)),'matched_fraction_mean':float(rows.matched_fraction.mean()),'seen_fraction_mean':float(rows.seen_fraction.mean()),
       'matching_strategy':'targeted same unordered log2 degree-bin nonedges; degrees from training edges only','interpretation_boundary':'Model-free structural diagnostic only.'}
    Path(args.out_json).write_text(json.dumps(s,indent=2)+'\n')
    gate=s['auc_random_mean']>=.60 and s['inflation_auc_mean']>=.08 and s['matched_fraction_mean']>=.50
    md=f"""# Gate 1 — HuRI PPI result\n\n- Edges: **{s['edge_count']:,}**\n- Proteins: **{s['nodes']:,}**\n- Degree-only AUC, conventional random negatives: **{s['auc_random_mean']:.3f} ± {s['auc_random_sd']:.3f}**\n- Degree-only AUC, targeted degree-matched negatives: **{s['auc_degree_matched_mean']:.3f} ± {s['auc_degree_matched_sd']:.3f}**\n- Mean AUC inflation: **{s['inflation_auc_mean']:.3f}**\n- Mean matched fraction: **{s['matched_fraction_mean']:.3f}**\n- Mean evaluable seen-endpoint fraction: **{s['seen_fraction_mean']:.3f}**\n- Gate-1 PPI structural-inflation signal: **{'PASS' if gate else 'FAIL / INCONCLUSIVE'}**\n\nMatching uses nonedges with the same unordered prespecified log2 degree-bin pair as each positive, with degrees computed from training edges only.\n"""
    Path(args.out_md).write_text(md); print(json.dumps(s,indent=2)); print(md)
if __name__=='__main__': main()
