#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from sklearn.preprocessing import StandardScaler
import run_graphban_targetdecagon_clean as gb


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--biosnap',required=True); ap.add_argument('--drug-map',required=True); ap.add_argument('--gene-map',required=True)
    ap.add_argument('--drug-features',required=True); ap.add_argument('--protein-features',required=True)
    ap.add_argument('--fit-index',type=int,required=True); ap.add_argument('--out',required=True)
    a=ap.parse_args()
    edges=gb.parse_edges(a.biosnap)
    dd=pd.read_csv(a.drug_map,dtype=str).fillna(''); gd=pd.read_csv(a.gene_map,dtype=str).fillna('')
    dd=dd[dd.smiles.astype(bool)].drop_duplicates('drug'); gd=gd[gd.sequence.astype(bool)].drop_duplicates('gene')
    dm=set(dd.drug); gm=set(gd.gene); mapped=sorted([e for e in edges if e[0] in dm and e[1] in gm])
    L=sorted({x for x,_ in mapped}); R=sorted({y for _,y in mapped}); blocked=set(mapped)
    Xd0=np.load(a.drug_features).astype('float32'); Xp0=np.load(a.protein_features).astype('float32')
    Xd=StandardScaler().fit_transform(Xd0).astype('float32'); Xp=StandardScaler().fit_transform(Xp0).astype('float32')
    device=torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    data,li,ri=gb.make_graph(L,R,mapped,Xd,Xp,device)
    seed=91000+a.fit_index
    model=gb.train_model(data,li,ri,mapped,blocked,L,R,seed,device,20)
    candidates=[(x,y) for x in L for y in R if (x,y) not in blocked]
    flat=np.asarray([li[x]*len(R)+ri[y] for x,y in candidates],dtype=np.int64)
    scores=gb.model_scores(model,data,candidates,li,ri,device,batch=8192).astype('float32')
    p=Path(a.out); p.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(p,flat=flat,scores=scores,fit_index=np.asarray([a.fit_index]),seed=np.asarray([seed]))
    print('fit',a.fit_index,'seed',seed,'candidates',len(candidates),'score_mean',float(scores.mean()))
if __name__=='__main__': main()
