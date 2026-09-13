#!/usr/bin/env python3
"""Reproduce the Anti-DDI ATC5 degree/popularity-bias diagnostic.

This analysis is structural. It does not adjudicate clinical non-interaction.
"""
from __future__ import annotations
import argparse, glob, json
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd

BINS=np.array([0,1,2,3,5,8,12,20,35,60,100,175,300,550,1000,np.inf],dtype=float)

def norm_pair(a,b):
    return tuple(sorted((str(a).strip(),str(b).strip())))

def auc_mw(pos_scores,neg_scores):
    p=np.asarray(pos_scores,float); n=np.asarray(neg_scores,float)
    ranks=pd.Series(np.concatenate([p,n])).rank(method='average').to_numpy()
    u=ranks[:len(p)].sum()-len(p)*(len(p)+1)/2
    return float(u/(len(p)*len(n)))

def load_positives(path=None, pattern=None):
    if path:
        d=pd.read_csv(path)
    else:
        files=sorted(glob.glob(pattern or ''))
        if not files: raise FileNotFoundError(f'no files matched: {pattern}')
        d=pd.concat([pd.read_csv(f) for f in files],ignore_index=True)
    return sorted({norm_pair(a,b) for a,b in zip(d.atc5_a,d.atc5_b) if str(a)!=str(b)})

def build_negative_atc5(df):
    usable=df[df.evidence_tier.isin(['T1_wellpowered','T2_moderate'])]
    out=set()
    for _,r in usable.iterrows():
        aa=[] if pd.isna(r.atc_a) else [x.strip()[:5] for x in str(r.atc_a).split(';') if x.strip()]
        bb=[] if pd.isna(r.atc_b) else [x.strip()[:5] for x in str(r.atc_b).split(';') if x.strip()]
        for a in aa:
            for b in bb:
                if len(a)==5 and len(b)==5 and a!=b: out.add(norm_pair(a,b))
    return sorted(out)

def main():
    ap=argparse.ArgumentParser()
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--positives')
    g.add_argument('--positives-glob')
    ap.add_argument('--antiddi',required=True)
    ap.add_argument('--out',required=True)
    ap.add_argument('--summary',required=True)
    args=ap.parse_args()
    pos=load_positives(args.positives,args.positives_glob); posset=set(pos)
    df=pd.read_csv(args.antiddi)
    neg_all=build_negative_atc5(df); overlap=set(neg_all)&posset
    neg=sorted(set(neg_all)-posset)
    deg=Counter()
    for a,b in pos: deg[a]+=1; deg[b]+=1
    vocab=sorted(deg); blocked=posset|set(neg)
    def score(p): return float(np.log1p(deg.get(p[0],0))*np.log1p(deg.get(p[1],0)))
    def sig(p):
        ix=[]
        for d in p:
            j=int(np.searchsorted(BINS,deg.get(d,0),side='right')-1)
            ix.append(max(0,min(j,len(BINS)-2)))
        return tuple(sorted(ix))
    rows=[]
    for seed in range(20):
        rng=np.random.default_rng(seed); n=min(len(neg),len(pos))
        P=[pos[i] for i in rng.choice(len(pos),size=n,replace=False)]
        rnd=[]; seen=set()
        while len(rnd)<len(neg):
            i,j=rng.integers(len(vocab),size=2)
            if i==j: continue
            q=norm_pair(vocab[int(i)],vocab[int(j)])
            if q in blocked or q in seen: continue
            seen.add(q); rnd.append(q)
        pools={}
        for q in neg: pools.setdefault(sig(q),[]).append(q)
        for k in pools: rng.shuffle(pools[k])
        used={}; mp=[]; mn=[]
        for ii in rng.permutation(len(P)):
            p=P[int(ii)]; s=sig(p); pool=pools.get(s,[]); u=used.get(s,0)
            if u>=len(pool): continue
            used[s]=u+1; mp.append(p); mn.append(pool[u])
        rows.append({'seed':seed,
          'auc_random_unlabelled':auc_mw([score(x) for x in P],[score(x) for x in rnd]),
          'auc_curated_atc5':auc_mw([score(x) for x in P],[score(x) for x in neg]),
          'auc_degree_matched':auc_mw([score(x) for x in mp],[score(x) for x in mn]),
          'n_positive_sample':len(P),'n_curated_atc5':len(neg),'n_matched':len(mp)})
    out=pd.DataFrame(rows); Path(args.out).parent.mkdir(parents=True,exist_ok=True); out.to_csv(args.out,index=False)
    d=out.auc_curated_atc5-out.auc_degree_matched; lo,hi=np.quantile(d,[.025,.975])
    s={'positive_atc5_pairs':len(pos),'positive_atc5_nodes':len(vocab),
       'curated_atc5_pairs_before_overlap_exclusion':len(neg_all),'class_pair_overlap_excluded':len(overlap),
       'curated_atc5_pairs_analyzed':len(neg),'replicates':20,
       'auc_random_mean':float(out.auc_random_unlabelled.mean()),'auc_random_sd':float(out.auc_random_unlabelled.std(ddof=1)),
       'auc_curated_mean':float(out.auc_curated_atc5.mean()),'auc_curated_sd':float(out.auc_curated_atc5.std(ddof=1)),
       'auc_degree_matched_mean':float(out.auc_degree_matched.mean()),'auc_degree_matched_sd':float(out.auc_degree_matched.std(ddof=1)),
       'curated_minus_matched_mean':float(d.mean()),'curated_minus_matched_empirical_95_interval':[float(lo),float(hi)],
       'n_matched_mean':float(out.n_matched.mean()),'n_matched_range':[int(out.n_matched.min()),int(out.n_matched.max())],
       'note':'ATC level-5 structural sensitivity analysis; not pair-level clinical validation.'}
    Path(args.summary).write_text(json.dumps(s,indent=2)+'\n'); print(json.dumps(s,indent=2))
if __name__=='__main__': main()
