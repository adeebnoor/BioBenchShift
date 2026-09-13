#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import roc_auc_score


def _col(df, *needles):
    norm = {c.lower().replace('#','').strip(): c for c in df.columns}
    for n in needles:
        n = n.lower()
        for k, c in norm.items():
            if n in k:
                return c
    raise KeyError(f"missing column matching {needles}; have {list(df.columns)[:20]}")


def load_human_edges(zip_path):
    with zipfile.ZipFile(zip_path) as z:
        names = [n for n in z.namelist() if not n.endswith('/')]
        if not names:
            raise ValueError(f"empty zip: {zip_path}")
        # MV-Physical zip should contain one tab3 text file.
        name = max(names, key=lambda n: z.getinfo(n).file_size)
        with z.open(name) as fh:
            d = pd.read_csv(fh, sep='\t', dtype=str, low_memory=False)
    a = _col(d, 'entrez gene interactor a')
    b = _col(d, 'entrez gene interactor b')
    oa = _col(d, 'organism id interactor a')
    ob = _col(d, 'organism id interactor b')
    q = d[(d[oa].astype(str) == '9606') & (d[ob].astype(str) == '9606')][[a,b]].dropna()
    out = set()
    for x,y in zip(q[a].astype(str), q[b].astype(str)):
        x=x.strip(); y=y.strip()
        if not x or not y or x=='-' or y=='-' or x==y:
            continue
        out.add(tuple(sorted((x,y))))
    return sorted(out), name, len(d), len(q)


def dbin(d):
    return 0 if d <= 0 else int(math.floor(math.log2(d))) + 1


def random_nonedges(rng, n, nodes, blocked):
    out=set(); lim=max(500000,n*100)
    while len(out)<n and lim>0:
        lim-=1
        i,j=rng.integers(len(nodes),size=2)
        if i==j: continue
        e=tuple(sorted((nodes[int(i)],nodes[int(j)])))
        if e not in blocked and e not in out: out.add(e)
    if len(out)<n: raise RuntimeError(f"random nonedges {len(out)}/{n}")
    return list(out)


def matched_nonedges(rng, positives, nodes, blocked, deg):
    bins=defaultdict(list)
    for x in nodes: bins[dbin(deg.get(x,0))].append(x)
    used=set(); P=[]; N=[]
    for ii in rng.permutation(len(positives)):
        a,b=positives[int(ii)]
        A=bins.get(dbin(deg.get(a,0)),[]); B=bins.get(dbin(deg.get(b,0)),[])
        z=None
        for _ in range(2000):
            if not A or not B: break
            x=A[int(rng.integers(len(A)))]; y=B[int(rng.integers(len(B)))]
            if x==y: continue
            e=tuple(sorted((x,y)))
            if e not in blocked and e not in used:
                z=e; break
        if z is not None:
            used.add(z); P.append((a,b)); N.append(z)
    return P,N


def fit_svd(edges,nodes,seed,dim=32):
    ix={x:i for i,x in enumerate(nodes)}
    r=[]; c=[]
    for a,b in edges:
        i,j=ix[a],ix[b]; r += [i,j]; c += [j,i]
    X=sparse.csr_matrix((np.ones(len(r)),(r,c)),shape=(len(nodes),len(nodes)))
    k=max(2,min(dim,len(nodes)-2))
    m=TruncatedSVD(n_components=k,random_state=seed)
    U=m.fit_transform(X); V=m.components_.T
    return ix,U,V


def score_svd(E,m):
    ix,U,V=m
    return np.asarray([0.5*(np.dot(U[ix[a]],V[ix[b]])+np.dot(U[ix[b]],V[ix[a]])) for a,b in E])


def sigmoid(x): return 1/(1+np.exp(-np.clip(x,-30,30)))


def fit_mf(edges,nodes,blocked,seed,dim=32,epochs=18,lr=.035,reg=2e-4):
    rng=np.random.default_rng(seed+91000); ix={x:i for i,x in enumerate(nodes)}
    U=rng.normal(0,.08,(len(nodes),dim)); bias=np.zeros(len(nodes))
    pos=np.asarray([(ix[a],ix[b]) for a,b in edges],dtype=int); bs=2048
    for _ in range(epochs):
        rng.shuffle(pos)
        for s in range(0,len(pos),bs):
            pp=pos[s:s+bs]; n=len(pp); neg=[]
            while len(neg)<n:
                i,j=rng.integers(len(nodes),size=2)
                if i==j: continue
                e=tuple(sorted((nodes[int(i)],nodes[int(j)])))
                if e not in blocked: neg.append((int(i),int(j)))
            nn=np.asarray(neg,dtype=int)
            ii=np.r_[pp[:,0],nn[:,0]]; jj=np.r_[pp[:,1],nn[:,1]]; y=np.r_[np.ones(n),np.zeros(n)]
            pred=sigmoid(np.sum(U[ii]*U[jj],axis=1)+bias[ii]+bias[jj]); g=pred-y
            ui=U[ii].copy(); uj=U[jj].copy()
            np.add.at(U,ii,-lr*(g[:,None]*uj+reg*ui)); np.add.at(U,jj,-lr*(g[:,None]*ui+reg*uj))
            np.add.at(bias,ii,-lr*g); np.add.at(bias,jj,-lr*g)
    return ix,U,bias


def score_mf(E,m):
    ix,U,b=m
    return np.asarray([np.dot(U[ix[a]],U[ix[c]])+b[ix[a]]+b[ix[c]] for a,c in E])


def fit_lightgcn(edges,nodes,blocked,seed,dim=32,layers=2,epochs=18,lr=.025,reg=1e-5):
    import torch
    torch.manual_seed(seed); rng=np.random.default_rng(seed+4567); ix={x:i for i,x in enumerate(nodes)}; n=len(nodes)
    src=[];dst=[];deg=np.zeros(n,dtype=np.float32)
    for a,b in edges:
        i,j=ix[a],ix[b];src += [i,j];dst += [j,i];deg[i]+=1;deg[j]+=1
    vals=[1/math.sqrt(max(deg[i],1)*max(deg[j],1)) for i,j in zip(src,dst)]
    A=torch.sparse_coo_tensor(torch.tensor([src,dst]),torch.tensor(vals,dtype=torch.float32),(n,n)).coalesce()
    emb=torch.nn.Embedding(n,dim);torch.nn.init.normal_(emb.weight,std=.08);opt=torch.optim.Adam(emb.parameters(),lr=lr,weight_decay=reg)
    pi=torch.tensor([ix[a] for a,b in edges]); pj=torch.tensor([ix[b] for a,b in edges])
    def prop():
        z=emb.weight; zs=[z]
        for _ in range(layers): z=torch.sparse.mm(A,z); zs.append(z)
        return torch.stack(zs).mean(0)
    for _ in range(epochs):
        ni=[];nj=[]
        while len(ni)<len(edges):
            i,j=rng.integers(n,size=2)
            if i==j: continue
            e=tuple(sorted((nodes[int(i)],nodes[int(j)])))
            if e not in blocked: ni.append(int(i));nj.append(int(j))
        ni=torch.tensor(ni);nj=torch.tensor(nj);z=prop();ps=(z[pi]*z[pj]).sum(1);ns=(z[ni]*z[nj]).sum(1)
        loss=torch.nn.functional.softplus(-ps).mean()+torch.nn.functional.softplus(ns).mean();opt.zero_grad();loss.backward();opt.step()
    with torch.no_grad(): Z=prop().cpu().numpy()
    return ix,Z


def score_lg(E,m):
    ix,Z=m
    return np.asarray([np.dot(Z[ix[a]],Z[ix[b]]) for a,b in E])


def auc(pos,neg,fn):
    return float(roc_auc_score(np.r_[np.ones(len(pos)),np.zeros(len(neg))],np.r_[fn(pos),fn(neg)]))


def fit_all(train,nodes,blocked,seed):
    return {
      'SVD': (fit_svd(train,nodes,seed),score_svd),
      'NeuralMF': (fit_mf(train,nodes,blocked,seed),score_mf),
      'LightGCN': (fit_lightgcn(train,nodes,blocked,seed),score_lg),
    }


def internal_eval(old,seeds):
    nodes=sorted({x for e in old for x in e}); blocked=set(old); rows=[]
    for seed in seeds:
        rng=np.random.default_rng(seed); perm=rng.permutation(len(old)); nt=max(1,round(.2*len(old)))
        test=[old[i] for i in perm[:nt]]; train=[old[i] for i in perm[nt:]]
        deg=Counter(x for e in train for x in e); seen=[e for e in test if deg[e[0]]>0 and deg[e[1]]>0]
        rnd=random_nonedges(rng,len(seen),nodes,blocked); mp,mn=matched_nonedges(rng,seen,nodes,blocked,deg)
        models=fit_all(train,nodes,blocked,seed)
        for name,(model,fn) in models.items():
            rows.append({'seed':seed,'model':name,'auc_random':auc(seen,rnd,lambda E,m=model,f=fn:f(E,m)),
                         'auc_matched':auc(mp,mn,lambda E,m=model,f=fn:f(E,m)),
                         'n_seen':len(seen),'n_matched':len(mp),'match_fraction':len(mp)/len(seen)})
    return pd.DataFrame(rows)


def future_controls(rng,future,nodes,blocked_current,deg_old):
    rnd=random_nonedges(rng,len(future),nodes,blocked_current)
    mp,mn=matched_nonedges(rng,future,nodes,blocked_current,deg_old)
    return rnd,mp,mn


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--old',required=True);ap.add_argument('--new',required=True);ap.add_argument('--out-prefix',required=True);ap.add_argument('--replicates',type=int,default=3);args=ap.parse_args()
    old,old_member,old_raw,old_h=load_human_edges(args.old); new,new_member,new_raw,new_h=load_human_edges(args.new)
    oldset=set(old);newset=set(new);nodes=sorted({x for e in old for x in e});nodeset=set(nodes)
    future=sorted(e for e in newset-oldset if e[0] in nodeset and e[1] in nodeset)
    if len(future)<100: raise RuntimeError(f'Only {len(future)} closed-world future edges')
    internal=internal_eval(old,list(range(args.replicates)))
    means=internal.groupby('model')[['auc_random','auc_matched']].mean()
    win_random=str(means.auc_random.idxmax());win_matched=str(means.auc_matched.idxmax())

    # Retrain every candidate model on the complete historical graph, without future edges.
    full_models=fit_all(old,nodes,oldset,100)
    deg_old=Counter(x for e in old for x in e);rng=np.random.default_rng(20260913)
    rnd,fp,fm=future_controls(rng,future,nodes,newset,deg_old)
    ext=[]
    for name,(model,fn) in full_models.items():
        ar=auc(future,rnd,lambda E,m=model,f=fn:f(E,m));am=auc(fp,fm,lambda E,m=model,f=fn:f(E,m))
        ext.append({'model':name,'future_auc_random_controls':ar,'future_auc_degree_matched_controls':am,'n_future':len(future),'n_future_matched':len(fp)})
    ext=pd.DataFrame(ext)

    # Compare whether internal benchmark scores identify future performance.
    merged=means.reset_index().merge(ext,on='model')
    rand_future=float(merged.loc[merged.model.eq(win_random),'future_auc_degree_matched_controls'].iloc[0])
    neut_future=float(merged.loc[merged.model.eq(win_matched),'future_auc_degree_matched_controls'].iloc[0])
    verdict='NEUTRALIZED_WINS_FUTURE' if neut_future>rand_future else ('TIE' if abs(neut_future-rand_future)<1e-12 else 'CONVENTIONAL_WINS_FUTURE')

    p=Path(args.out_prefix);p.parent.mkdir(parents=True,exist_ok=True);internal.to_csv(str(p)+'_internal_replicates.csv',index=False);merged.to_csv(str(p)+'_model_table.csv',index=False)
    summary={'old_release':'5.0.250','new_release':'5.0.261','old_member':old_member,'new_member':new_member,
      'old_human_mv_edges':len(old),'new_human_mv_edges':len(new),'historical_nodes':len(nodes),'closed_world_future_edges':len(future),
      'internal_replicates':args.replicates,'internal_matching_fraction_mean':float(internal.match_fraction.mean()),
      'conventional_winner':win_random,'neutralized_winner':win_matched,
      'conventional_winner_future_matched_auc':rand_future,'neutralized_winner_future_matched_auc':neut_future,
      'future_auc_difference_neutralized_minus_conventional':neut_future-rand_future,'verdict':verdict,
      'interpretation':'Future edges are relations absent from BioGRID 5.0.250 and present in 5.0.261, restricted to endpoints already present historically. Persistent unobserved pairs are controls, not proven biological negatives.'}
    Path(str(p)+'_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# H4 pilot — BioGRID temporal external validity','',f"Historical release: **5.0.250**; later release: **5.0.261**.",f"Historical human multi-validated physical network: **{len(old):,} edges / {len(nodes):,} nodes**.",f"Closed-world later-added edges: **{len(future):,}**.",'','## Model selection and future performance','', '| Model | Internal random AUC | Internal degree-matched AUC | Future AUC vs degree-matched persistent controls |','|---|---:|---:|---:|']
    for _,r in merged.sort_values('model').iterrows(): lines.append(f"| {r.model} | {r.auc_random:.3f} | {r.auc_matched:.3f} | {r.future_auc_degree_matched_controls:.3f} |")
    lines += ['',f"Conventional benchmark winner: **{win_random}**.",f"Structure-neutralized winner: **{win_matched}**.",f"Future matched-control AUC difference (neutralized-selected minus conventional-selected): **{neut_future-rand_future:+.3f}**.",f"Pilot verdict: **{verdict}**.",'','Future-added BioGRID edges provide temporal external evidence, not an exhaustive truth set; persistent unobserved pairs are not asserted to be true negatives.']
    md='\n'.join(lines)+'\n';Path(str(p)+'.md').write_text(md);print(md)

if __name__=='__main__': main()
