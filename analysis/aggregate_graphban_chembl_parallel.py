#!/usr/bin/env python3
from __future__ import annotations
import argparse, glob, json
from pathlib import Path
import numpy as np
import pandas as pd
import run_chembl_dti_external as ext
import run_graphban_targetdecagon_clean as gb


def load_pairs(path):
    d=pd.read_csv(path,dtype=str); return set(zip(d.drug,d.gene))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--biosnap',required=True); ap.add_argument('--drug-map',required=True); ap.add_argument('--gene-map',required=True); ap.add_argument('--score-glob',required=True); ap.add_argument('--prep-prefix',required=True); ap.add_argument('--out-prefix',required=True); a=ap.parse_args()
    edges=gb.parse_edges(a.biosnap); dd=pd.read_csv(a.drug_map,dtype=str).fillna(''); gd=pd.read_csv(a.gene_map,dtype=str).fillna('')
    dm=set(dd.loc[dd.smiles.astype(bool),'drug']); gm=set(gd.loc[gd.sequence.astype(bool),'gene']); mapped=sorted([e for e in edges if e[0] in dm and e[1] in gm])
    L=sorted({x for x,_ in mapped}); R=sorted({y for _,y in mapped}); li={x:i for i,x in enumerate(L)}; ri={x:i for i,x in enumerate(R)}; blocked=set(mapped)
    files=sorted(glob.glob(a.score_glob));
    if len(files)!=5: raise RuntimeError(f'expected 5 GraphBAN score artifacts, got {files}')
    flats=[]; scores=[]; fits=[]
    for f in files:
        z=np.load(f); flats.append(z['flat']); scores.append(z['scores'].astype('float32')); fits.append(int(z['fit_index'][0]))
    order=np.argsort(fits); flats=[flats[i] for i in order]; scores=[scores[i] for i in order]; fits=[fits[i] for i in order]
    if fits!=[0,1,2,3,4]: raise RuntimeError(f'fit indices {fits}')
    if any(not np.array_equal(flats[0],x) for x in flats[1:]): raise RuntimeError('candidate flat order differs across GraphBAN fits')
    flat=flats[0]; S=np.vstack(scores); mean=S.mean(axis=0); ranked=flat[np.argsort(mean)[::-1]]
    s6=load_pairs(a.prep_prefix+'_support6.csv'); s7=load_pairs(a.prep_prefix+'_support7.csv'); sexact=load_pairs(a.prep_prefix+'_support_exact.csv')
    sf6=ext.flat_support(s6,li,ri,len(R)); sf7=ext.flat_support(s7,li,ri,len(R)); sfexact=ext.flat_support(sexact,li,ri,len(R)); eligible=len(flat)
    gp=pd.DataFrame(ext.metrics_for_rank('GraphBAN',ranked,sf6,eligible)); gs=[]
    for label,S0 in [('pchembl_ge_7',sf7),('exact_relation_pchembl_ge_6',sfexact)]:
        for r in ext.metrics_for_rank('GraphBAN',ranked,S0,eligible): r['sensitivity']=label; gs.append(r)
    baseline=pd.read_csv(a.prep_prefix+'_baseline_primary.csv'); primary=pd.concat([baseline,gp],ignore_index=True)
    bs=pd.read_csv(a.prep_prefix+'_baseline_sensitivity.csv'); sensitivity=pd.concat([bs,pd.DataFrame(gs)],ignore_index=True)
    loo={}
    for k in (100,500,1000):
        top=set(ranked[:k].tolist()); vals=[]
        for i in range(5):
            m=np.delete(S,i,axis=0).mean(axis=0); r=flat[np.argsort(m)[::-1]]; vals.append(1-len(top & set(r[:k].tolist()))/k)
        loo[str(k)]={'mean':float(np.mean(vals)),'max':float(np.max(vals)),'values':[float(v) for v in vals]}
    prep=json.load(open(a.prep_prefix+'_meta.json'))
    p=Path(a.out_prefix); p.parent.mkdir(parents=True,exist_ok=True); primary.to_csv(str(p)+'_primary_table.csv',index=False); sensitivity.to_csv(str(p)+'_sensitivity_table.csv',index=False)
    summary={**prep,'ensemble_fits_per_model':5,'selected_model_conventional':'GraphBAN','selected_model_neutralized':'GraphBAN','benchmark_induced_model_switch':False,'graphban_training_seeds':[91000,91001,91002,91003,91004],'graphban_loo_turnover':loo,'graphban_primary':gp[['k','hits','precision_lower_bound','recall_supported','enrichment_uniform']].to_dict('records')}
    Path(str(p)+'_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    piv=primary.pivot(index='model',columns='k',values='hits'); lines=['# Fair-universe GraphBAN -> ChEMBL 37 follow-through','',f"Mapped TargetDecagon coverage: **{prep['mapping_fraction']:.1%}**. Candidate universe: **{prep['candidate_universe_size']:,}** unknown mapped pairs.",'','GraphBAN is selected under both frozen DTI evaluation regimes; therefore the contemporary DTI panel has **no benchmark-induced model-identity switch**.','','| Model | @100 | @500 | @1,000 | @5,000 | @10,000 | @50,000 |','|---|---:|---:|---:|---:|---:|---:|']
    for model in ['SVD','NeuralMF','LightGCN','GraphBAN']:
        row=piv.loc[model]; vals=[int(row.get(k,0)) for k in (100,500,1000,5000,10000,50000)]; lines.append('| '+model+' | '+' | '.join(map(str,vals))+' |')
    lines += ['',f"GraphBAN leave-one-fit-out turnover: HT@100 **{loo['100']['mean']:.3f}**, HT@500 **{loo['500']['mean']:.3f}**, HT@1000 **{loo['1000']['mean']:.3f}**.",'','**Boundary:** ChEMBL characterizes the independently supported content of the common GraphBAN frontier; it does not create a DTI benchmark-induced switch because GraphBAN is the selected model under both evaluation rules. The historical PPI experiment remains the full benchmark -> winner -> hypotheses -> later-evidence chain.']
    Path(str(p)+'.md').write_text('\n'.join(lines)+'\n'); print('\n'.join(lines))
if __name__=='__main__': main()
