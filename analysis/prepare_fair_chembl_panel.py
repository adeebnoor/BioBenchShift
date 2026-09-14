#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd
import numpy as np
import run_chembl_dti_external as ext
import run_graphban_targetdecagon_clean as gb


def save_pairs(path, pairs):
    pd.DataFrame(sorted(pairs),columns=['drug','gene']).to_csv(path,index=False)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--biosnap',required=True); ap.add_argument('--drug-map',required=True); ap.add_argument('--gene-map',required=True); ap.add_argument('--out-prefix',required=True); ap.add_argument('--inits',type=int,default=5); a=ap.parse_args()
    status=ext.get_json('https://www.ebi.ac.uk/chembl/api/data/status.json'); version=status.get('chembl_db_version') or status.get('version')
    if version!='ChEMBL_37': raise RuntimeError(f'frozen protocol requires ChEMBL_37, got {version!r}')
    edges=gb.parse_edges(a.biosnap)
    dd=pd.read_csv(a.drug_map,dtype=str).fillna(''); gd=pd.read_csv(a.gene_map,dtype=str).fillna('')
    dm=set(dd.loc[dd.smiles.astype(bool),'drug']); gm=set(gd.loc[gd.sequence.astype(bool),'gene'])
    mapped=sorted([e for e in edges if e[0] in dm and e[1] in gm]); L=sorted({x for x,_ in mapped}); R=sorted({y for _,y in mapped})
    drug_map=ext.map_pubchem_to_chembl(L); _,gene_targets=ext.target_maps(R)
    s6,s7,sexact,emeta,raw=ext.evidence_sets(edges,drug_map,gene_targets,6.0)
    L0,R0,scores=ext.fit_ensembles(mapped,a.inits)
    if L0!=L or R0!=R: raise RuntimeError('baseline canonical universe mismatch')
    li={x:i for i,x in enumerate(L)}; ri={x:i for i,x in enumerate(R)}; mask=np.ones((len(L),len(R)),dtype=bool)
    for x,y in mapped: mask[li[x],ri[y]]=False
    eligible=int(mask.sum()); sf6=ext.flat_support(s6,li,ri,len(R)); sf7=ext.flat_support(s7,li,ri,len(R)); sfexact=ext.flat_support(sexact,li,ri,len(R))
    rows=[]; sens=[]
    for name,sc in scores.items():
        ranked=ext.ranked_indices(sc,mask); rows += ext.metrics_for_rank(name,ranked,sf6,eligible)
        for label,S in [('pchembl_ge_7',sf7),('exact_relation_pchembl_ge_6',sfexact)]:
            for r in ext.metrics_for_rank(name,ranked,S,eligible): r['sensitivity']=label; sens.append(r)
    p=Path(a.out_prefix); p.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(str(p)+'_baseline_primary.csv',index=False); pd.DataFrame(sens).to_csv(str(p)+'_baseline_sensitivity.csv',index=False)
    save_pairs(str(p)+'_support6.csv',s6); save_pairs(str(p)+'_support7.csv',s7); save_pairs(str(p)+'_support_exact.csv',sexact)
    raw.to_csv(str(p)+'_eligible_activity_records.csv',index=False)
    meta={'chembl_version':version,'edge_count_full':len(edges),'edge_count_mapped':len(mapped),'mapping_fraction':len(mapped)/len(edges),'candidate_universe_size':eligible,'evidence_meta':emeta,'supported_primary':len(sf6),'supported_strict':len(sf7)}
    Path(str(p)+'_meta.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(meta,indent=2))
if __name__=='__main__': main()
