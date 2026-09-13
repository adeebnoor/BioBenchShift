#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

import run_biogrid_temporal_h4 as h4
import run_biogrid_future_yield as fy


def read_tab3(zip_path):
    with zipfile.ZipFile(zip_path) as z:
        names=[n for n in z.namelist() if not n.endswith('/')]
        name=max(names,key=lambda n:z.getinfo(n).file_size)
        with z.open(name) as fh:
            return pd.read_csv(fh,sep='\t',dtype=str,low_memory=False),name


def exact_ranks_for_pairs(model_name, model, nodes, known_edges, pairs, block=256):
    ix={x:i for i,x in enumerate(nodes)}; n=len(nodes)
    rr=[];cc=[]
    for a,b in known_edges:
        i,j=ix[a],ix[b]; rr += [i,j]; cc += [j,i]
    A=sparse.csr_matrix((np.ones(len(rr),dtype=np.int8),(rr,cc)),shape=(n,n))
    fn={'SVD':h4.score_svd,'NeuralMF':h4.score_mf,'LightGCN':h4.score_lg}[model_name]
    target_scores=np.asarray(fn(pairs,model),dtype=np.float64)
    greater=np.zeros(len(pairs),dtype=np.int64); ties=np.zeros(len(pairs),dtype=np.int64)
    for start in range(0,n,block):
        stop=min(n,start+block); rows=np.arange(start,stop)
        S=np.asarray(fy.score_block(model_name,model,rows),dtype=np.float64)
        for r,i in enumerate(rows):
            S[r,:i+1]=-np.inf
            js=A.indices[A.indptr[i]:A.indptr[i+1]]
            if len(js): S[r,js]=-np.inf
        vals=S[np.isfinite(S)]
        for j,t in enumerate(target_scores):
            greater[j]+=int(np.count_nonzero(vals>t))
            ties[j]+=int(np.count_nonzero(vals==t))
    return target_scores,greater+1,ties


def annotate_new_rows(new_zip, selected):
    d,member=read_tab3(new_zip)
    a=h4._col(d,'entrez gene interactor a'); b=h4._col(d,'entrez gene interactor b')
    oa=h4._col(d,'organism id interactor a'); ob=h4._col(d,'organism id interactor b')
    human=d[(d[oa].astype(str)=='9606') & (d[ob].astype(str)=='9606')].copy()
    wanted={tuple(sorted(x)) for x in selected}
    human['_pair']=[tuple(sorted((str(x).strip(),str(y).strip()))) for x,y in zip(human[a],human[b])]
    q=human[human['_pair'].isin(wanted)].copy()
    def opt(*needles):
        try:return h4._col(q,*needles)
        except Exception:return None
    sa=opt('official symbol interactor a'); sb=opt('official symbol interactor b')
    pub=opt('publication source'); exp=opt('experimental system'); expt=opt('experimental system type'); thr=opt('throughput'); iid=opt('biogrid interaction id')
    out={}
    for pair,g in q.groupby('_pair',sort=False):
        p0,p1=pair; syms=[]
        if sa and sb:
            for _,r in g.iterrows():
                if str(r[a]).strip()==p0: syms.append((str(r[sa]),str(r[sb])))
                else: syms.append((str(r[sb]),str(r[sa])))
        def uniq(c):
            if not c:return ''
            return '; '.join(sorted({str(v).strip() for v in g[c].dropna() if str(v).strip() and str(v).strip()!='-'}))
        out[pair]={
          'symbol_a': next((x[0] for x in syms if x[0] not in ('','-','nan')),''),
          'symbol_b': next((x[1] for x in syms if x[1] not in ('','-','nan')),''),
          'later_record_count':int(len(g)),
          'publication_sources':uniq(pub),'experimental_systems':uniq(exp),
          'experimental_system_types':uniq(expt),'throughput':uniq(thr),
          'biogrid_interaction_ids':uniq(iid),'new_member':member,
        }
    return out


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--old',required=True); ap.add_argument('--new',required=True); ap.add_argument('--out-prefix',required=True); args=ap.parse_args()
    old,old_member,_,_=h4.load_human_edges(args.old); new,new_member,_,_=h4.load_human_edges(args.new)
    oldset=set(old); newset=set(new); nodes=sorted({x for e in old for x in e}); nodeset=set(nodes); ix={x:i for i,x in enumerate(nodes)}; n=len(nodes)
    future=sorted(e for e in newset-oldset if e[0] in nodeset and e[1] in nodeset); fset=set(future)
    models=h4.fit_all(old,nodes,oldset,100)

    svd_model=models['SVD'][0]
    ids,scores=fy.global_topk('SVD',svd_model,nodes,oldset,100)
    svd_top=[tuple(sorted((nodes[int(q//n)],nodes[int(q%n)]))) for q in ids]
    selected=[p for p in svd_top if p in fset]
    neural_ids,_=fy.global_topk('NeuralMF',models['NeuralMF'][0],nodes,oldset,100)
    neural_top=[tuple(sorted((nodes[int(q//n)],nodes[int(q%n)]))) for q in neural_ids]
    neural_future=[p for p in neural_top if p in fset]
    if len(selected)!=7:
        raise RuntimeError(f'Frozen SVD Top-100 later-hit count mismatch: expected 7, got {len(selected)}')
    if len(neural_future)!=0:
        raise RuntimeError(f'Frozen NeuralMF Top-100 later-hit count mismatch: expected 0, got {len(neural_future)}')

    svd_rank={p:i+1 for i,p in enumerate(svd_top)}
    ranks={}; scoremap={}; tiemap={}
    for name in ['NeuralMF','LightGCN']:
        sc,rk,ti=exact_ranks_for_pairs(name,models[name][0],nodes,oldset,selected)
        ranks[name]=dict(zip(selected,map(int,rk))); scoremap[name]=dict(zip(selected,map(float,sc))); tiemap[name]=dict(zip(selected,map(int,ti)))
    ann=annotate_new_rows(args.new,selected)

    rows=[]
    for p in sorted(selected,key=lambda x:svd_rank[x]):
        z=ann.get(p,{})
        rows.append({'gene_a':p[0],'gene_b':p[1],'symbol_a':z.get('symbol_a',''),'symbol_b':z.get('symbol_b',''),
          'svd_rank':svd_rank[p],'neuralmf_rank':ranks['NeuralMF'][p],'lightgcn_rank':ranks['LightGCN'][p],
          'neuralmf_score_tie_count':tiemap['NeuralMF'][p],'lightgcn_score_tie_count':tiemap['LightGCN'][p],
          'later_record_count':z.get('later_record_count',0),'publication_sources':z.get('publication_sources',''),
          'experimental_systems':z.get('experimental_systems',''),'experimental_system_types':z.get('experimental_system_types',''),
          'throughput':z.get('throughput',''),'biogrid_interaction_ids':z.get('biogrid_interaction_ids','')})
    out=pd.DataFrame(rows); p=Path(args.out_prefix); p.parent.mkdir(parents=True,exist_ok=True); out.to_csv(str(p)+'_table.csv',index=False)
    summary={'protocol':'BIOGRID_NAMED_LATER_CASES_AMENDMENT_20260913.md','baseline_commit':'1298e4c83ff4fa488a784d52219d20e706118a9e',
      'old_release':'5.0.250','new_release':'5.0.261','old_member':old_member,'new_member':new_member,'later_added_closed_world_edges':len(future),
      'svd_top100_later_supported_count':len(selected),'neuralmf_top100_later_supported_count':len(neural_future),
      'selection_rule':'all SVD Top-100 pairs in frozen later-added closed-world set; no manual selection',
      'rank_definition':'1 + number of eligible historical candidate pairs with strictly greater score; tie count reported separately'}
    Path(str(p)+'_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# Named later-supported BioGRID cases','',f"Mechanical case set: **{len(selected)}** later-added interactions in the frozen SVD Top-100; NeuralMF reciprocal Top-100 later-added count: **{len(neural_future)}**.",'',
      '| SVD rank | Gene A | Gene B | Symbols | NeuralMF rank | LightGCN rank | Later publication/source | Experimental system |','|---:|---:|---:|---|---:|---:|---|---|']
    for _,r in out.iterrows():
        lines.append(f"| {int(r.svd_rank)} | {r.gene_a} | {r.gene_b} | {r.symbol_a}–{r.symbol_b} | {int(r.neuralmf_rank):,} | {int(r.lightgcn_rank):,} | {r.publication_sources} | {r.experimental_systems} |")
    lines += ['','These are later-supported database interactions, not wet-lab discoveries made by this study. Pair membership was fixed mechanically before annotation.']
    Path(str(p)+'.md').write_text('\n'.join(lines)+'\n'); print('\n'.join(lines))

if __name__=='__main__': main()
