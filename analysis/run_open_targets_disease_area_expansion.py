#!/usr/bin/env python3
"""Open Targets disease-area context for benchmark-selected DTI queues.

Protocol: OPEN_TARGETS_DISEASE_AREA_AMENDMENT_20260913.md
The aggregation, release, equal-target weighting, and bootstrap rule are frozen
before any Open Targets disease-area outcome is inspected.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.spatial.distance import jensenshannon

import run_biosnap_dti_gate2_models as dti
import run_biological_consequence_expansion as bio

KS=(100,500,1000)
MODELS=("SVD","NeuralMF")
GCONVERT="https://biit.cs.ut.ee/gprofiler/api/convert/convert/"
OT="https://api.platform.opentargets.org/api/v4/graphql"
PAGE_SIZE=3000


def post_json(url,payload,tries=6,timeout=180):
    err=None
    for i in range(tries):
        try:
            r=requests.post(url,json=payload,timeout=timeout,headers={"User-Agent":"BioBenchShift-Science-expansion/1.0"})
            r.raise_for_status(); j=r.json()
            if isinstance(j,dict) and j.get("errors"):
                raise RuntimeError(str(j["errors"])[:2000])
            return j
        except Exception as e:
            err=e; time.sleep(min(2**i,20))
    raise RuntimeError(f"POST failed {url}: {err}")


def ot_meta():
    q="query { meta { apiVersion { x y z } dataVersion { year month } } }"
    j=post_json(OT,{"query":q,"variables":{}})
    m=j.get("data",{}).get("meta") or {}
    dv=m.get("dataVersion") or {}; rel=f"{dv.get('year','')}.{dv.get('month','')}"
    if rel!="26.06":
        raise RuntimeError(f"Protocol requires Open Targets data 26.06, got {rel!r}: {m}")
    return m


def convert_entrez(ids):
    payload={"organism":"hsapiens","target":"ENSG","numeric_ns":"ENTREZGENE_ACC","query":[str(x) for x in sorted(set(ids),key=lambda x:(len(str(x)),str(x)))],"output":"json"}
    j=post_json(GCONVERT,payload)
    by=defaultdict(set)
    for r in j.get("result",[]):
        incoming=str(r.get("incoming") or '').strip(); conv=str(r.get("converted") or '').strip()
        if incoming and conv and conv not in {"N/A","None","nan"} and conv.startswith("ENSG"):
            by[incoming].add(conv.split('.')[0])
    rows=[]; chosen={}
    for g in payload["query"]:
        vals=sorted(by.get(g,set())); pick=vals[0] if vals else None
        if pick: chosen[g]=pick
        rows.append({"entrez_gene":g,"ensembl_gene":pick or "","n_ensembl_mappings":len(vals),"all_ensembl_mappings":";".join(vals)})
    return chosen,pd.DataFrame(rows),j.get("meta",{})


ASSOC_QUERY="""
query TargetAreas($id:String!,$index:Int!,$size:Int!){
  target(ensemblId:$id){
    id
    approvedSymbol
    associatedDiseases(enableIndirect:false,page:{index:$index,size:$size}){
      count
      rows{
        score
        disease{
          id
          name
          therapeuticAreas{ id name }
        }
      }
    }
  }
}
"""


def fetch_target(eid):
    page=0; rows=[]; symbol=None; total=None
    while True:
        j=post_json(OT,{"query":ASSOC_QUERY,"variables":{"id":eid,"index":page,"size":PAGE_SIZE}})
        t=(j.get("data") or {}).get("target")
        if t is None:
            return {"ensembl_gene":eid,"symbol":"","count":0,"rows":[],"unresolved":True}
        symbol=t.get("approvedSymbol") or ''
        a=t.get("associatedDiseases") or {}; total=int(a.get("count") or 0); z=a.get("rows") or []; rows.extend(z)
        if len(rows)>=total or not z: break
        page+=1
        if page>100: raise RuntimeError(f"pagination runaway for {eid}: {len(rows)}/{total}")
    return {"ensembl_gene":eid,"symbol":symbol,"count":total,"rows":rows,"unresolved":False}


def target_area_vector(rec):
    best={}; area_names={}
    for r in rec.get("rows",[]):
        try:s=float(r.get("score") or 0.0)
        except Exception:s=0.0
        if not np.isfinite(s) or s<=0: continue
        disease=r.get("disease") or {}
        for a in disease.get("therapeuticAreas") or []:
            aid=str(a.get("id") or '').strip(); name=str(a.get("name") or '').strip()
            if not aid: continue
            if s>best.get(aid,0.0): best[aid]=s
            if name: area_names[aid]=name
    tot=sum(best.values())
    if tot<=0: return {},area_names
    return {k:v/tot for k,v in best.items()},area_names


def js(a,b,areas):
    x=np.asarray([a.get(k,0.0) for k in areas],float); y=np.asarray([b.get(k,0.0) for k in areas],float)
    if x.sum()<=0 or y.sum()<=0:return float('nan')
    x/=x.sum(); y/=y.sum(); return float(jensenshannon(x,y,base=2.0)**2)


def aggregate(vectors, genes):
    vv=[vectors[g] for g in genes if g in vectors and vectors[g]]
    if not vv:return {},0
    out=defaultdict(float)
    for v in vv:
        for a,w in v.items(): out[a]+=float(w)
    n=len(vv)
    return {a:w/n for a,w in out.items()},n


def boot_js(vectors,a_genes,b_genes,areas,n=2000,seed=0):
    av=[vectors[g] for g in a_genes if g in vectors and vectors[g]]; bv=[vectors[g] for g in b_genes if g in vectors and vectors[g]]
    if not av or not bv:return [float('nan')]*3
    arr=np.asarray(areas,object); rng=np.random.default_rng(seed); vals=[]
    for _ in range(n):
        aa=[av[int(i)] for i in rng.integers(len(av),size=len(av))]; bb=[bv[int(i)] for i in rng.integers(len(bv),size=len(bv))]
        def mean(vs):
            z=defaultdict(float)
            for v in vs:
                for k,w in v.items():z[k]+=w
            return {k:w/len(vs) for k,w in z.items()}
        vals.append(js(mean(aa),mean(bb),areas))
    q=np.nanpercentile(np.asarray(vals,float),[2.5,50,97.5])
    return [float(q[1]),float(q[0]),float(q[2])]


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--biosnap',required=True); ap.add_argument('--out-prefix',required=True); ap.add_argument('--inits',type=int,default=5); ap.add_argument('--workers',type=int,default=6); a=ap.parse_args()
    meta=ot_meta()
    edges=dti.parse(a.biosnap); background,sets=bio.target_sets(edges,a.inits)
    all_entrez=sorted({g for z in sets.values() for g in z},key=lambda x:(len(str(x)),str(x)))
    mapping,mapdf,gmeta=convert_entrez(all_entrez)
    eids=sorted(set(mapping.values()))

    records={}
    with ThreadPoolExecutor(max_workers=max(1,a.workers)) as ex:
        fut={ex.submit(fetch_target,eid):eid for eid in eids}
        for f in as_completed(fut):
            rec=f.result(); records[rec['ensembl_gene']]=rec

    vectors={}; area_names={}; target_assoc_rows=[]
    for entrez,eid in mapping.items():
        rec=records.get(eid,{})
        v,names=target_area_vector(rec); vectors[entrez]=v; area_names.update(names)
        target_assoc_rows.append({"entrez_gene":entrez,"ensembl_gene":eid,"approved_symbol":rec.get('symbol',''),"direct_association_count":rec.get('count',0),"n_therapeutic_areas":len(v),"has_disease_area_annotation":bool(v)})
    all_areas=sorted(area_names)

    profile_rows=[]; summary=[]
    for k in KS:
        genes_by={m:list(sets[(m,k)]) for m in MODELS}
        profiles={}; annotated={}
        for m in MODELS:
            prof,n=aggregate(vectors,genes_by[m]); profiles[m]=prof; annotated[m]=n
            for ar in all_areas:
                profile_rows.append({"model":m,"k":k,"therapeutic_area_id":ar,"therapeutic_area_name":area_names.get(ar,''),"mass":prof.get(ar,0.0)})
        j=js(profiles['SVD'],profiles['NeuralMF'],all_areas)
        top={m:[x for x,_ in sorted(profiles[m].items(),key=lambda z:(-z[1],z[0]))[:5]] for m in MODELS}
        denom=min(len(top['SVD']),len(top['NeuralMF'])); overlap=(len(set(top['SVD'])&set(top['NeuralMF']))/denom) if denom else float('nan')
        med,lo,hi=boot_js(vectors,genes_by['SVD'],genes_by['NeuralMF'],all_areas,2000,20260913+k)
        summary.append({"k":k,"svd_unique_targets":len(genes_by['SVD']),"neuralmf_unique_targets":len(genes_by['NeuralMF']),
          "svd_mapped_targets":sum(g in mapping for g in genes_by['SVD']),"neuralmf_mapped_targets":sum(g in mapping for g in genes_by['NeuralMF']),
          "svd_annotated_targets":annotated['SVD'],"neuralmf_annotated_targets":annotated['NeuralMF'],
          "svd_annotation_fraction":annotated['SVD']/len(genes_by['SVD']) if genes_by['SVD'] else float('nan'),
          "neuralmf_annotation_fraction":annotated['NeuralMF']/len(genes_by['NeuralMF']) if genes_by['NeuralMF'] else float('nan'),
          "js_divergence_disease_area_profile":j,"bootstrap_js_median":med,"bootstrap_js_95ci_low":lo,"bootstrap_js_95ci_high":hi,
          "top5_area_overlap_coefficient":overlap,"svd_top5_areas":";".join(top['SVD']),"neuralmf_top5_areas":";".join(top['NeuralMF'])})

    p=Path(a.out_prefix); p.parent.mkdir(parents=True,exist_ok=True)
    mapdf.merge(pd.DataFrame(target_assoc_rows),on=['entrez_gene','ensembl_gene'],how='left').to_csv(str(p)+'_mapping.csv',index=False)
    pd.DataFrame(profile_rows).to_csv(str(p)+'_profiles.csv',index=False)
    sdf=pd.DataFrame(summary); sdf.to_csv(str(p)+'_summary.csv',index=False)
    target_rows=[]
    for m in MODELS:
        for k in KS:
            for g in sets[(m,k)]:
                for ar,w in vectors.get(g,{}).items(): target_rows.append({"model":m,"k":k,"entrez_gene":g,"ensembl_gene":mapping.get(g,''),"therapeutic_area_id":ar,"therapeutic_area_name":area_names.get(ar,''),"target_normalized_area_weight":w})
    pd.DataFrame(target_rows).to_csv(str(p)+'_target_vectors.csv',index=False)
    md={"protocol":"OPEN_TARGETS_DISEASE_AREA_AMENDMENT_20260913.md","baseline_commit":"1298e4c83ff4fa488a784d52219d20e706118a9e","open_targets_endpoint":OT,"open_targets_meta":meta,"gprofiler_convert_endpoint":GCONVERT,"gprofiler_meta":gmeta,"requested_open_targets_release":"26.06","total_unique_queue_targets":len(all_entrez),"mapped_ensembl_targets":len(mapping),"unique_ensembl_targets_queried":len(eids),"n_therapeutic_areas":len(all_areas),"bootstrap_replicates":2000,"aggregation":"per target: max direct association score within therapeutic area, normalize target vector to unit mass, equal-weight average across annotated targets","interpretation_boundary":"Disease-area profiles measure different allocation of attention across Open Targets therapeutic areas; they do not validate individual DTI predictions or establish clinical efficacy."}
    Path(str(p)+'_metadata.json').write_text(json.dumps(md,indent=2,default=str)+'\n')
    print(sdf.to_string(index=False))

if __name__=='__main__':main()
