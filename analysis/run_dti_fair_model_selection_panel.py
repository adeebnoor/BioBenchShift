#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import roc_auc_score, average_precision_score

# Reuse the exact GraphBAN split/control construction helpers.
from run_graphban_targetdecagon_clean import parse_edges, sample_random, degree_matched


def sigmoid(x):
    return 1/(1+np.exp(-np.clip(x,-30,30)))


def fit_svd(train,L,R,k=32,seed=0):
    li={x:i for i,x in enumerate(L)}; ri={x:i for i,x in enumerate(R)}
    rr=[li[a] for a,b in train]; cc=[ri[b] for a,b in train]
    X=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(L),len(R)))
    svd=TruncatedSVD(n_components=min(k,min(X.shape)-1),random_state=seed)
    U=svd.fit_transform(X); V=svd.components_.T
    return li,ri,U,V

def score_svd(E,m):
    li,ri,U,V=m
    return np.asarray([np.dot(U[li[a]],V[ri[b]]) for a,b in E],dtype=float)


def fit_logistic_mf(train,L,R,blocked,rng,k=32,epochs=25,lr=.04,reg=2e-4):
    li={x:i for i,x in enumerate(L)};ri={x:i for i,x in enumerate(R)}
    U=rng.normal(0,.08,(len(L),k));V=rng.normal(0,.08,(len(R),k));bu=np.zeros(len(L));bv=np.zeros(len(R))
    pos=np.array([(li[a],ri[b]) for a,b in train],dtype=int);bs=1024
    for _ in range(epochs):
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

def score_mf(E,m):
    li,ri,U,V,bu,bv=m
    return np.asarray([np.dot(U[li[a]],V[ri[b]])+bu[li[a]]+bv[ri[b]] for a,b in E],dtype=float)


def fit_lightgcn(train,left,right,blocked,seed,dim=32,layers=2,epochs=35,lr=.03,reg=1e-5):
    import torch
    torch.manual_seed(seed)
    rng=np.random.default_rng(seed+12345)
    li={x:i for i,x in enumerate(left)};ri={x:i+len(left) for i,x in enumerate(right)}
    n=len(left)+len(right);src=[];dst=[];deg=np.zeros(n,dtype=np.float32)
    for a,b in train:
        i,j=li[a],ri[b];src += [i,j];dst += [j,i];deg[i]+=1;deg[j]+=1
    vals=[1.0/math.sqrt(max(deg[i],1.0)*max(deg[j],1.0)) for i,j in zip(src,dst)]
    A=torch.sparse_coo_tensor(torch.tensor([src,dst]),torch.tensor(vals,dtype=torch.float32),(n,n)).coalesce()
    emb=torch.nn.Embedding(n,dim);torch.nn.init.normal_(emb.weight,std=.08)
    opt=torch.optim.Adam(emb.parameters(),lr=lr,weight_decay=reg)
    pi=torch.tensor([li[a] for a,_ in train]);pj=torch.tensor([ri[b] for _,b in train])
    def propagated():
        z0=emb.weight;zs=[z0];z=z0
        for _ in range(layers): z=torch.sparse.mm(A,z);zs.append(z)
        return torch.stack(zs).mean(0)
    for _ in range(epochs):
        na=[];nb=[]
        while len(na)<len(train):
            a=left[int(rng.integers(len(left)))];b=right[int(rng.integers(len(right)))]
            if (a,b) not in blocked: na.append(li[a]);nb.append(ri[b])
        ni=torch.tensor(na);nj=torch.tensor(nb);z=propagated();ps=(z[pi]*z[pj]).sum(1);ns=(z[ni]*z[nj]).sum(1)
        loss=torch.nn.functional.softplus(-ps).mean()+torch.nn.functional.softplus(ns).mean();opt.zero_grad();loss.backward();opt.step()
    with torch.no_grad(): z=propagated().cpu().numpy()
    return li,ri,z

def score_lg(E,m):
    li,ri,z=m
    return np.asarray([np.dot(z[li[a]],z[ri[b]]) for a,b in E],dtype=float)


def metrics(pos,neg,scorefn):
    sp=scorefn(pos);sn=scorefn(neg);y=np.r_[np.ones(len(sp)),np.zeros(len(sn))];s=np.r_[sp,sn]
    return float(roc_auc_score(y,s)),float(average_precision_score(y,s))


def one_seed(E_all,E_mapped,L,R,seed):
    rng=np.random.default_rng(seed);ix=rng.permutation(len(E_all));nt=round(.2*len(E_all))
    test_all={E_all[i] for i in ix[:nt]};train_all={E_all[i] for i in ix[nt:]};mapped=set(E_mapped)
    train=sorted(train_all & mapped);test=sorted(test_all & mapped);blocked=set(E_mapped)
    dl=Counter(a for a,b in train);dr=Counter(b for a,b in train);seen=[e for e in test if dl[e[0]]>0 and dr[e[1]]>0]
    # EXACT same evaluation RNG as the frozen GraphBAN challenge.
    erng=np.random.default_rng(seed+40000);rn=sample_random(erng,len(seen),L,R,blocked);mp,mn=degree_matched(erng,seen,L,R,blocked,dl,dr)
    svd=fit_svd(train,L,R,32,seed)
    nmf=fit_logistic_mf(train,L,R,blocked,np.random.default_rng(seed+1000))
    lg=fit_lightgcn(train,L,R,blocked,seed)
    out={'seed':seed,'n_train':len(train),'n_test_mapped':len(test),'n_seen':len(seen),'n_matched':len(mp),'match_fraction':len(mp)/len(seen)}
    for name,fn in [('svd',lambda z:score_svd(z,svd)),('neuralmf',lambda z:score_mf(z,nmf)),('lightgcn',lambda z:score_lg(z,lg))]:
        ar,pr=metrics(seen,rn,fn);am,pm=metrics(mp,mn,fn)
        out.update({f'{name}_auc_random':ar,f'{name}_auprc_random':pr,f'{name}_auc_matched':am,f'{name}_auprc_matched':pm,f'{name}_auc_drop':ar-am})
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--drug-map',required=True);ap.add_argument('--gene-map',required=True);ap.add_argument('--graphban-replicates',required=True);ap.add_argument('--out-prefix',required=True);ap.add_argument('--seeds',default='0,1,2');a=ap.parse_args()
    E=parse_edges(a.input);dd=pd.read_csv(a.drug_map,dtype=str).fillna('');gd=pd.read_csv(a.gene_map,dtype=str).fillna('')
    dm=set(dd.loc[dd.smiles.astype(bool),'drug']);gm=set(gd.loc[gd.sequence.astype(bool),'gene']);mapped=[e for e in E if e[0] in dm and e[1] in gm]
    L=sorted({x for x,_ in mapped});R=sorted({y for _,y in mapped});seeds=[int(x) for x in a.seeds.split(',') if x.strip()]
    out=pd.DataFrame([one_seed(E,mapped,L,R,s) for s in seeds])
    gb=pd.read_csv(a.graphban_replicates);gb=gb[gb.seed.isin(seeds)].copy()
    if set(gb.seed)!=set(seeds): raise RuntimeError('GraphBAN replicate seeds do not match requested panel seeds')
    for c in ['n_train','n_test_mapped','n_seen','n_matched']:
        z=out[['seed',c]].merge(gb[['seed',c]],on='seed',suffixes=('_base','_gb'))
        if not (z[f'{c}_base']==z[f'{c}_gb']).all(): raise RuntimeError(f'GraphBAN control mismatch in {c}: {z.to_dict("records")}')
    for metric in ['auc_random','auprc_random','auc_matched','auprc_matched','auc_drop']:
        out=out.merge(gb[['seed',metric]].rename(columns={metric:f'graphban_{metric}'}),on='seed',how='left')
    models=['svd','neuralmf','lightgcn','graphban']
    for regime in ['random','matched']:
        cols=[f'{m}_auc_{regime}' for m in models]
        out[f'winner_{regime}']=out[cols].idxmax(axis=1).str.replace(f'_auc_{regime}','',regex=False)
    p=Path(a.out_prefix);p.parent.mkdir(parents=True,exist_ok=True);out.to_csv(str(p)+'_replicates.csv',index=False)
    rows=[]
    for m in models:
        rows.append({'model':m,'auc_random_mean':float(out[f'{m}_auc_random'].mean()),'auc_random_sd':float(out[f'{m}_auc_random'].std(ddof=1)),'auc_matched_mean':float(out[f'{m}_auc_matched'].mean()),'auc_matched_sd':float(out[f'{m}_auc_matched'].std(ddof=1)),'auprc_random_mean':float(out[f'{m}_auprc_random'].mean()),'auprc_matched_mean':float(out[f'{m}_auprc_matched'].mean())})
    tab=pd.DataFrame(rows).sort_values('auc_random_mean',ascending=False);tab.to_csv(str(p)+'_model_table.csv',index=False)
    random_winner=tab.sort_values('auc_random_mean',ascending=False).iloc[0].model
    matched_winner=tab.sort_values('auc_matched_mean',ascending=False).iloc[0].model
    summary={'benchmark':'BioSNAP TargetDecagon DTI mapped GraphBAN universe','edge_count_full':len(E),'edge_count_mapped':len(mapped),'mapping_fraction':len(mapped)/len(E),'seeds':seeds,'match_fraction_mean':float(out.match_fraction.mean()),'winner_random':random_winner,'winner_matched':matched_winner,'winner_reversal':bool(random_winner!=matched_winner),'seedwise_winners_random':out.winner_random.tolist(),'seedwise_winners_matched':out.winner_matched.tolist(),'models':rows}
    Path(str(p)+'_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# Fair DTI model-selection panel on the GraphBAN-mapped universe','','All four models are compared on the same mapped TargetDecagon relation universe, identical frozen splits, identical conventional controls, and identical degree-matched controls used by the leakage-free GraphBAN S6 challenge.','','| Model | Conventional AUROC | Neutralized AUROC | Conventional AUPRC | Neutralized AUPRC |','|---|---:|---:|---:|---:|']
    for r in rows: lines.append(f"| {r['model']} | {r['auc_random_mean']:.3f} ± {r['auc_random_sd']:.3f} | {r['auc_matched_mean']:.3f} ± {r['auc_matched_sd']:.3f} | {r['auprc_random_mean']:.3f} | {r['auprc_matched_mean']:.3f} |")
    lines += ['',f"Mean matching coverage: **{summary['match_fraction_mean']:.3f}**.",f"Conventional winner by mean AUROC: **{random_winner}**.",f"Structure-neutralized winner by mean AUROC: **{matched_winner}**.",f"Full-panel winner reversal: **{'YES' if summary['winner_reversal'] else 'NO'}**.",f"Seed-wise conventional winners: {summary['seedwise_winners_random']}.",f"Seed-wise neutralized winners: {summary['seedwise_winners_matched']}.",'','This panel supersedes comparisons that mixed different relation universes or evaluation-control draws. It does not erase the historical baseline-ladder result; it determines the contemporary full-panel model-selection claim.']
    Path(str(p)+'.md').write_text('\n'.join(lines)+'\n');print('\n'.join(lines))
if __name__=='__main__':main()
