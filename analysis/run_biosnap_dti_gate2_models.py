#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import roc_auc_score
from scipy import sparse

def dbin(d):return 0 if d<=0 else int(math.floor(math.log2(d)))+1
def parse(path):
 d=pd.read_csv(path,compression='infer',comment='#',header=None,names=['drug','gene'],sep=',',skip_blank_lines=True,dtype=str).dropna()
 return sorted({(str(a).strip(),str(b).strip()) for a,b in zip(d.drug,d.gene)})
def sample_random(rng,n,L,R,blocked):
 o=set()
 while len(o)<n:
  e=(L[int(rng.integers(len(L)))],R[int(rng.integers(len(R)))])
  if e not in blocked:o.add(e)
 return list(o)
def matched(rng,pos,L,R,blocked,dl,dr):
 lb=defaultdict(list);rb=defaultdict(list)
 for x in L:lb[dbin(dl.get(x,0))].append(x)
 for x in R:rb[dbin(dr.get(x,0))].append(x)
 used=set();P=[];N=[]
 for i in rng.permutation(len(pos)):
  p=pos[int(i)];A=lb[dbin(dl[p[0]])];B=rb[dbin(dr[p[1]])];z=None
  for _ in range(1000):
   e=(A[int(rng.integers(len(A)))],B[int(rng.integers(len(B)))])
   if e not in blocked and e not in used:z=e;break
  if z:used.add(z);P.append(p);N.append(z)
 return P,N
def sigmoid(x):return 1/(1+np.exp(-np.clip(x,-30,30)))
def fit_logistic_mf(train,L,R,blocked,rng,k=32,epochs=25,lr=.04,reg=2e-4):
 li={x:i for i,x in enumerate(L)};ri={x:i for i,x in enumerate(R)};U=rng.normal(0,.08,(len(L),k));V=rng.normal(0,.08,(len(R),k));bu=np.zeros(len(L));bv=np.zeros(len(R))
 pos=np.array([(li[a],ri[b]) for a,b in train],dtype=int);bs=1024
 for ep in range(epochs):
  rng.shuffle(pos)
  for s in range(0,len(pos),bs):
   pp=pos[s:s+bs];n=len(pp);neg=[]
   while len(neg)<n:
    a=L[int(rng.integers(len(L)))];b=R[int(rng.integers(len(R)))];e=(a,b)
    if e not in blocked:neg.append((li[a],ri[b]))
   nn=np.array(neg,dtype=int);ii=np.r_[pp[:,0],nn[:,0]];jj=np.r_[pp[:,1],nn[:,1]];y=np.r_[np.ones(n),np.zeros(n)]
   pred=sigmoid(np.sum(U[ii]*V[jj],axis=1)+bu[ii]+bv[jj]);g=pred-y
   gu=g[:,None]*V[jj]+reg*U[ii];gv=g[:,None]*U[ii]+reg*V[jj]
   np.add.at(U,ii,-lr*gu);np.add.at(V,jj,-lr*gv);np.add.at(bu,ii,-lr*g);np.add.at(bv,jj,-lr*g)
 return li,ri,U,V,bu,bv
def score_mf(E,model):
 li,ri,U,V,bu,bv=model;return np.array([np.dot(U[li[a]],V[ri[b]])+bu[li[a]]+bv[ri[b]] for a,b in E])
def fit_svd(train,L,R,k=32,seed=0):
 li={x:i for i,x in enumerate(L)};ri={x:i for i,x in enumerate(R)};rr=[li[a] for a,b in train];cc=[ri[b] for a,b in train];X=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(L),len(R)))
 svd=TruncatedSVD(n_components=min(k,min(X.shape)-1),random_state=seed);U=svd.fit_transform(X);V=svd.components_.T
 return li,ri,U,V
def score_svd(E,m):
 li,ri,U,V=m;return np.array([np.dot(U[li[a]],V[ri[b]]) for a,b in E])
def A(pos,neg,fn):return roc_auc_score(np.r_[np.ones(len(pos)),np.zeros(len(neg))],np.r_[fn(pos),fn(neg)])
def one(E,seed):
 rng=np.random.default_rng(seed);ix=rng.permutation(len(E));nt=round(.2*len(E));te=[E[i] for i in ix[:nt]];tr=[E[i] for i in ix[nt:]];blocked=set(E);L=sorted({a for a,b in E});R=sorted({b for a,b in E});dl=Counter(a for a,b in tr);dr=Counter(b for a,b in tr);seen=[e for e in te if dl[e[0]]>0 and dr[e[1]]>0];rn=sample_random(rng,len(seen),L,R,blocked);mp,mn=matched(rng,seen,L,R,blocked,dl,dr)
 sm=fit_svd(tr,L,R,32,seed);nm=fit_logistic_mf(tr,L,R,blocked,np.random.default_rng(seed+1000))
 ar_s=A(seen,rn,lambda z:score_svd(z,sm));am_s=A(mp,mn,lambda z:score_svd(z,sm));ar_n=A(seen,rn,lambda z:score_mf(z,nm));am_n=A(mp,mn,lambda z:score_mf(z,nm))
 return {'seed':seed,'n_seen':len(seen),'n_matched':len(mp),'match_fraction':len(mp)/len(seen),'svd_auc_random':ar_s,'svd_auc_matched':am_s,'svd_drop':ar_s-am_s,'neuralmf_auc_random':ar_n,'neuralmf_auc_matched':am_n,'neuralmf_drop':ar_n-am_n}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--out-prefix',required=True);ap.add_argument('--replicates',type=int,default=5);a=ap.parse_args();E=parse(a.input);o=pd.DataFrame([one(E,s) for s in range(a.replicates)]);p=Path(a.out_prefix);p.parent.mkdir(parents=True,exist_ok=True);o.to_csv(str(p)+'_replicates.csv',index=False)
 s={'benchmark':'BioSNAP TargetDecagon DTI','edge_count':len(E),'replicates':a.replicates,'match_fraction_mean':float(o.match_fraction.mean())}
 for m in ['svd','neuralmf']:
  for x in ['auc_random','auc_matched','drop']:s[f'{m}_{x}_mean']=float(o[f'{m}_{x}'].mean());s[f'{m}_{x}_sd']=float(o[f'{m}_{x}'].std(ddof=1))
 Path(str(p)+'_summary.json').write_text(json.dumps(s,indent=2)+'\n');md=f"""# Gate 2 pilot — learned models on BioSNAP DTI\n\n| Model | Random-negative AUC | Degree-matched AUC | Mean drop |\n|---|---:|---:|---:|\n| Truncated-SVD latent factors | {s['svd_auc_random_mean']:.3f} ± {s['svd_auc_random_sd']:.3f} | {s['svd_auc_matched_mean']:.3f} ± {s['svd_auc_matched_sd']:.3f} | {s['svd_drop_mean']:.3f} |\n| Logistic neural matrix factorization | {s['neuralmf_auc_random_mean']:.3f} ± {s['neuralmf_auc_random_sd']:.3f} | {s['neuralmf_auc_matched_mean']:.3f} ± {s['neuralmf_auc_matched_sd']:.3f} | {s['neuralmf_drop_mean']:.3f} |\n\nMean degree-matching coverage: **{s['match_fraction_mean']:.3f}**. These are structural latent models, a Gate-2 pilot rather than the final SOTA model panel.\n""";Path(str(p)+'.md').write_text(md);print(md)
if __name__=='__main__':main()
