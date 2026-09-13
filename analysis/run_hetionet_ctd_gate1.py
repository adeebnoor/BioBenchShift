#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

def auc_mw(pos,neg):
    p=np.asarray(pos,float); n=np.asarray(neg,float)
    if not len(p) or not len(n): return float('nan')
    r=pd.Series(np.concatenate([p,n])).rank(method='average').to_numpy(); u=r[:len(p)].sum()-len(p)*(len(p)+1)/2
    return float(u/(len(p)*len(n)))
def dbin(d): return 0 if d<=0 else int(math.floor(math.log2(d)))+1
def score(e,dl,dr): return float(np.log1p(dl.get(e[0],0))*np.log1p(dr.get(e[1],0)))
def parse(path):
    d=pd.read_csv(path,sep='\t',compression='infer',dtype=str)
    cols=[c.lower() for c in d.columns]; d.columns=cols
    if not {'source','metaedge','target'}.issubset(d.columns): raise ValueError(f'Unexpected columns {d.columns.tolist()}')
    q=d[d.metaedge.eq('CtD')].copy()
    edges=sorted({(str(a).strip(),str(b).strip()) for a,b in zip(q.source,q.target) if str(a).startswith('Compound::') and str(b).startswith('Disease::')})
    if len(edges)<100: raise ValueError(f'Only {len(edges)} CtD edges')
    return edges
def random_nonedges(rng,n,L,R,blocked):
    out=set(); lim=max(100000,n*500); k=0
    while len(out)<n and k<lim:
        k+=1; e=(L[int(rng.integers(len(L)))],R[int(rng.integers(len(R)))])
        if e in blocked or e in out: continue
        out.add(e)
    if len(out)<n: raise RuntimeError(f'Only {len(out)}/{n} nonedges')
    return list(out)
def matched(rng,pos,L,R,blocked,dl,dr):
    lb=defaultdict(list); rb=defaultdict(list)
    for x in L: lb[dbin(dl.get(x,0))].append(x)
    for x in R: rb[dbin(dr.get(x,0))].append(x)
    used=set(); mp=[]; mn=[]
    for i in rng.permutation(len(pos)):
        p=pos[int(i)]; A=lb.get(dbin(dl.get(p[0],0)),[]); B=rb.get(dbin(dr.get(p[1],0)),[])
        found=None
        for _ in range(1000):
            if not A or not B: break
            e=(A[int(rng.integers(len(A)))],B[int(rng.integers(len(B)))])
            if e not in blocked and e not in used: found=e; break
        if found: used.add(found); mp.append(p); mn.append(found)
    return mp,mn
def one(edges,seed):
    rng=np.random.default_rng(seed); idx=rng.permutation(len(edges)); nt=max(1,int(round(.2*len(edges))))
    test=[edges[i] for i in idx[:nt]]; train=[edges[i] for i in idx[nt:]]; allset=set(edges)
    dl=Counter(a for a,_ in train); dr=Counter(b for _,b in train); L=sorted({a for a,_ in edges}); R=sorted({b for _,b in edges})
    seen=[e for e in test if dl.get(e[0],0)>0 and dr.get(e[1],0)>0]
    rnd=random_nonedges(rng,len(seen),L,R,allset); mp,mn=matched(rng,seen,L,R,allset,dl,dr)
    ar=auc_mw([score(e,dl,dr) for e in seen],[score(e,dl,dr) for e in rnd]); am=auc_mw([score(e,dl,dr) for e in mp],[score(e,dl,dr) for e in mn]) if mp else float('nan')
    return {'seed':seed,'n_seen':len(seen),'seen_fraction':len(seen)/len(test),'n_matched':len(mp),'matched_fraction':len(mp)/len(seen) if seen else 0,'auc_random':ar,'auc_degree_matched':am,'inflation_auc':ar-am}
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input',required=True); ap.add_argument('--out-csv',required=True); ap.add_argument('--out-json',required=True); ap.add_argument('--out-md',required=True); ap.add_argument('--replicates',type=int,default=50); a=ap.parse_args()
    edges=parse(a.input); out=pd.DataFrame([one(edges,s) for s in range(a.replicates)]); Path(a.out_csv).parent.mkdir(parents=True,exist_ok=True); out.to_csv(a.out_csv,index=False)
    s={'benchmark':'Hetionet v1.0 Compound-treats-Disease (CtD)','edge_count':len(edges),'compounds':len({x for x,_ in edges}),'diseases':len({y for _,y in edges}),'replicates':a.replicates,
       'auc_random_mean':float(out.auc_random.mean()),'auc_random_sd':float(out.auc_random.std(ddof=1)),'auc_degree_matched_mean':float(out.auc_degree_matched.mean()),'auc_degree_matched_sd':float(out.auc_degree_matched.std(ddof=1)),'inflation_auc_mean':float(out.inflation_auc.mean()),'inflation_auc_sd':float(out.inflation_auc.std(ddof=1)),'matched_fraction_mean':float(out.matched_fraction.mean()),'seen_fraction_mean':float(out.seen_fraction.mean()),'matching_strategy':'same role-specific log2 training-degree bins'}
    Path(a.out_json).write_text(json.dumps(s,indent=2)+'\n'); gate=s['auc_random_mean']>=.60 and s['inflation_auc_mean']>=.08 and s['matched_fraction_mean']>=.50
    md=f"""# Gate 1 — Hetionet CtD result\n\n- Edges: **{s['edge_count']:,}**\n- Compounds: **{s['compounds']:,}**\n- Diseases: **{s['diseases']:,}**\n- Degree-only AUC, conventional random negatives: **{s['auc_random_mean']:.3f} ± {s['auc_random_sd']:.3f}**\n- Degree-only AUC, targeted degree-matched negatives: **{s['auc_degree_matched_mean']:.3f} ± {s['auc_degree_matched_sd']:.3f}**\n- Mean AUC inflation: **{s['inflation_auc_mean']:.3f}**\n- Mean matched fraction: **{s['matched_fraction_mean']:.3f}**\n- Mean evaluable seen-endpoint fraction: **{s['seen_fraction_mean']:.3f}**\n- Gate-1 drug–disease structural-inflation signal: **{'PASS' if gate else 'FAIL / INCONCLUSIVE'}**\n"""; Path(a.out_md).write_text(md); print(json.dumps(s,indent=2)); print(md)
if __name__=='__main__': main()
