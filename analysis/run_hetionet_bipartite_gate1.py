#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np,pandas as pd

def auc(pos,neg):
 p=np.asarray(pos,float);n=np.asarray(neg,float)
 if not len(p) or not len(n):return float('nan')
 r=pd.Series(np.r_[p,n]).rank(method='average').to_numpy();u=r[:len(p)].sum()-len(p)*(len(p)+1)/2
 return float(u/(len(p)*len(n)))
def bn(d):return 0 if d<=0 else int(math.floor(math.log2(d)))+1
def sc(e,dl,dr):return float(np.log1p(dl.get(e[0],0))*np.log1p(dr.get(e[1],0)))
def load(path,meta,left_prefix,right_prefix):
 d=pd.read_csv(path,sep='\t',compression='infer',dtype=str);d.columns=[c.lower() for c in d.columns]
 q=d[d.metaedge.eq(meta)]
 E=sorted({(str(a).strip(),str(b).strip()) for a,b in zip(q.source,q.target) if str(a).startswith(left_prefix) and str(b).startswith(right_prefix)})
 if len(E)<100:raise ValueError(f'{meta}: only {len(E)} edges')
 return E
def rnd_non(rng,n,L,R,blocked):
 out=set();k=0;lim=max(200000,n*200)
 while len(out)<n and k<lim:
  k+=1;e=(L[int(rng.integers(len(L)))],R[int(rng.integers(len(R)))])
  if e not in blocked and e not in out:out.add(e)
 if len(out)<n:raise RuntimeError(f'nonedges {len(out)}/{n}')
 return list(out)
def match(rng,pos,L,R,blocked,dl,dr):
 LB=defaultdict(list);RB=defaultdict(list)
 for x in L:LB[bn(dl.get(x,0))].append(x)
 for x in R:RB[bn(dr.get(x,0))].append(x)
 used=set();mp=[];mn=[]
 for i in rng.permutation(len(pos)):
  p=pos[int(i)];A=LB.get(bn(dl.get(p[0],0)),[]);B=RB.get(bn(dr.get(p[1],0)),[]);z=None
  for _ in range(800):
   if not A or not B:break
   e=(A[int(rng.integers(len(A)))],B[int(rng.integers(len(B)))])
   if e not in blocked and e not in used:z=e;break
  if z:used.add(z);mp.append(p);mn.append(z)
 return mp,mn
def rep(E,seed):
 rng=np.random.default_rng(seed);ix=rng.permutation(len(E));nt=max(1,round(.2*len(E)));te=[E[i] for i in ix[:nt]];tr=[E[i] for i in ix[nt:]];blocked=set(E)
 dl=Counter(a for a,_ in tr);dr=Counter(b for _,b in tr);L=sorted({a for a,_ in E});R=sorted({b for _,b in E});seen=[e for e in te if dl.get(e[0],0)>0 and dr.get(e[1],0)>0]
 rn=rnd_non(rng,len(seen),L,R,blocked);mp,mn=match(rng,seen,L,R,blocked,dl,dr);ar=auc([sc(e,dl,dr) for e in seen],[sc(e,dl,dr) for e in rn]);am=auc([sc(e,dl,dr) for e in mp],[sc(e,dl,dr) for e in mn])
 return {'seed':seed,'n_seen':len(seen),'seen_fraction':len(seen)/len(te),'n_matched':len(mp),'matched_fraction':len(mp)/len(seen) if seen else 0,'auc_random':ar,'auc_degree_matched':am,'inflation_auc':ar-am}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--metaedge',required=True);ap.add_argument('--left-prefix',required=True);ap.add_argument('--right-prefix',required=True);ap.add_argument('--label',required=True);ap.add_argument('--out-prefix',required=True);ap.add_argument('--replicates',type=int,default=20);a=ap.parse_args()
 E=load(a.input,a.metaedge,a.left_prefix,a.right_prefix);o=pd.DataFrame([rep(E,s) for s in range(a.replicates)]);p=Path(a.out_prefix);p.parent.mkdir(parents=True,exist_ok=True);o.to_csv(str(p)+'_replicates.csv',index=False)
 s={'benchmark':a.label,'metaedge':a.metaedge,'edge_count':len(E),'left_nodes':len({x for x,_ in E}),'right_nodes':len({y for _,y in E}),'replicates':a.replicates,'auc_random_mean':float(o.auc_random.mean()),'auc_random_sd':float(o.auc_random.std(ddof=1)),'auc_degree_matched_mean':float(o.auc_degree_matched.mean()),'auc_degree_matched_sd':float(o.auc_degree_matched.std(ddof=1)),'inflation_auc_mean':float(o.inflation_auc.mean()),'inflation_auc_sd':float(o.inflation_auc.std(ddof=1)),'matched_fraction_mean':float(o.matched_fraction.mean()),'seen_fraction_mean':float(o.seen_fraction.mean())};Path(str(p)+'_summary.json').write_text(json.dumps(s,indent=2)+'\n')
 gate=s['auc_random_mean']>=.60 and s['inflation_auc_mean']>=.08 and s['matched_fraction_mean']>=.50
 md=f"# Gate 1 — {a.label}\n\n- Edges: **{s['edge_count']:,}**\n- Left nodes: **{s['left_nodes']:,}**\n- Right nodes: **{s['right_nodes']:,}**\n- Degree-only AUC, random negatives: **{s['auc_random_mean']:.3f} ± {s['auc_random_sd']:.3f}**\n- Degree-only AUC, degree-matched negatives: **{s['auc_degree_matched_mean']:.3f} ± {s['auc_degree_matched_sd']:.3f}**\n- Mean AUC inflation: **{s['inflation_auc_mean']:.3f}**\n- Mean matched fraction: **{s['matched_fraction_mean']:.3f}**\n- Mean seen-endpoint fraction: **{s['seen_fraction_mean']:.3f}**\n- Structural-inflation signal: **{'PASS' if gate else 'FAIL / INCONCLUSIVE'}**\n";Path(str(p)+'.md').write_text(md);print(md)
if __name__=='__main__':main()
