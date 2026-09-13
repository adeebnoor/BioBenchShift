#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import time
from collections import defaultdict
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests

import run_biosnap_dti_gate2_models as dti
import run_lightgcn_gate2 as lgcn

KS=(100,500,1000,5000,10000,50000)
SESSION=requests.Session()
SESSION.headers.update({'User-Agent':'ANTI-DDI-Science-external-validation/1.0'})


def get_json(url, params=None, tries=6, timeout=90):
    err=None
    for i in range(tries):
        try:
            r=SESSION.get(url,params=params,timeout=timeout)
            if r.status_code in (429,500,502,503,504):
                time.sleep(min(2**i,20)); continue
            r.raise_for_status(); return r.json()
        except Exception as e:
            err=e; time.sleep(min(2**i,20))
    raise RuntimeError(f'GET failed {url}: {err}')


def post_json(url,data,tries=6,timeout=90):
    err=None
    for i in range(tries):
        try:
            r=SESSION.post(url,data=data,timeout=timeout)
            if r.status_code in (429,500,502,503,504):
                time.sleep(min(2**i,20)); continue
            r.raise_for_status(); return r.json()
        except Exception as e:
            err=e; time.sleep(min(2**i,20))
    raise RuntimeError(f'POST failed {url}: {err}')


def chunks(xs,n):
    xs=list(xs)
    for i in range(0,len(xs),n): yield xs[i:i+n]


def parse_cid(x):
    s=str(x).strip()
    if s.upper().startswith('CID'): s=s[3:]
    return str(int(s))


def pubchem_inchikeys(cids):
    out={}
    for ch in chunks(cids,80):
        url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+','.join(ch)+'/property/InChIKey/JSON'
        j=get_json(url)
        for p in j.get('PropertyTable',{}).get('Properties',[]):
            out[str(p['CID'])]=p.get('ConnectivitySMILES') and p.get('InChIKey') or p.get('InChIKey')
    return out


def map_pubchem_to_chembl(drugs):
    cid_to_drug={parse_cid(d):d for d in drugs}
    ik=pubchem_inchikeys(cid_to_drug)
    drug_to=set()
    out=defaultdict(set)
    for cid,key in ik.items():
        if not key: continue
        j=get_json(f'https://www.ebi.ac.uk/unichem/rest/verbose_inchikey/{key}')
        for rec in j if isinstance(j,list) else []:
            if rec.get('name')=='chembl' or int(rec.get('src_id',-1))==1:
                vals=rec.get('src_compound_id') or []
                if isinstance(vals,str): vals=[vals]
                for v in vals:
                    if str(v).startswith('CHEMBL'): out[cid_to_drug[cid]].add(str(v))
    return {k:sorted(v) for k,v in out.items()}


def uniprot_map(from_db,to_db,ids):
    ids=sorted(set(map(str,ids)))
    if not ids: return []
    j=post_json('https://rest.uniprot.org/idmapping/run',{'from':from_db,'to':to_db,'ids':','.join(ids)})
    job=j['jobId']
    for _ in range(120):
        s=get_json(f'https://rest.uniprot.org/idmapping/status/{job}')
        if s.get('jobStatus') in ('NEW','RUNNING'):
            time.sleep(2); continue
        if s.get('jobStatus') and s.get('jobStatus')!='FINISHED': raise RuntimeError(str(s))
        break
    details=get_json(f'https://rest.uniprot.org/idmapping/details/{job}')
    url=details['redirectURL']
    # Use the paged endpoint to avoid stressing the stream service.
    sep='&' if '?' in url else '?'
    url=url+sep+'size=500'
    allres=[]
    while url:
        r=SESSION.get(url,timeout=120);r.raise_for_status();j=r.json();allres.extend(j.get('results',[]))
        nxt=None
        link=r.headers.get('Link','')
        if 'rel="next"' in link:
            nxt=link.split('<',1)[1].split('>',1)[0]
        url=nxt
    return allres


def target_maps(genes):
    # GeneID -> human UniProt accessions.
    r1=uniprot_map('GeneID','UniProtKB',genes)
    gene_to_acc=defaultdict(set)
    for rec in r1:
        g=str(rec.get('from'))
        to=rec.get('to')
        if isinstance(to,dict):
            org=to.get('organism') or {}
            tax=org.get('taxonId')
            acc=to.get('primaryAccession') or to.get('uniProtkbId') or to.get('id')
            if tax not in (None,9606): continue
        else: acc=str(to) if to else None
        if acc: gene_to_acc[g].add(str(acc))
    accs=sorted({a for v in gene_to_acc.values() for a in v})
    r2=uniprot_map('UniProtKB_AC-ID','ChEMBL',accs)
    acc_to_target=defaultdict(set)
    for rec in r2:
        a=str(rec.get('from'));to=rec.get('to')
        if isinstance(to,dict): t=to.get('id') or to.get('primaryAccession')
        else: t=to
        if t and str(t).startswith('CHEMBL'): acc_to_target[a].add(str(t))
    gene_to_target=defaultdict(set)
    for g,aa in gene_to_acc.items():
        for a in aa: gene_to_target[g].update(acc_to_target.get(a,set()))
    return {k:sorted(v) for k,v in gene_to_acc.items()},{k:sorted(v) for k,v in gene_to_target.items()}


def chembl_pages(endpoint,params=None):
    base='https://www.ebi.ac.uk/chembl/api/data/'
    url=urljoin(base,endpoint)
    first=True
    while url:
        j=get_json(url,params=params if first else None);first=False
        key=endpoint.split('.',1)[0]+'s'
        if key not in j:
            # ChEMBL plural keys are not perfectly regular.
            for k,v in j.items():
                if k!='page_meta' and isinstance(v,list): key=k;break
        for x in j.get(key,[]): yield x
        nxt=(j.get('page_meta') or {}).get('next')
        url=urljoin('https://www.ebi.ac.uk',nxt) if nxt else None


def fetch_target_meta(target_ids):
    out={}
    for ch in chunks(sorted(target_ids),40):
        vals=','.join(ch)
        for x in chembl_pages('target.json',{'target_chembl_id__in':vals,'limit':1000}):
            tid=x.get('target_chembl_id')
            if tid: out[str(tid)]=x
    return out


def fetch_activities(mol_ids,min_p=6.0):
    out=[]
    for ch in chunks(sorted(mol_ids),12):
        params={'molecule_chembl_id__in':','.join(ch),'pchembl_value__gte':str(min_p),'limit':1000}
        out.extend(chembl_pages('activity.json',params))
    return out


def fetch_assays(assay_ids):
    out={}
    for ch in chunks(sorted(assay_ids),40):
        for x in chembl_pages('assay.json',{'assay_chembl_id__in':','.join(ch),'limit':1000}):
            aid=x.get('assay_chembl_id')
            if aid: out[str(aid)]=x
    return out


def valid_data_comment(x):
    v=x.get('data_validity_comment')
    return v in (None,'','None','null')


def evidence_sets(edges,drug_map,gene_targets,primary=6.0):
    all_mols={m for vs in drug_map.values() for m in vs}
    all_tgts={t for vs in gene_targets.values() for t in vs}
    tmeta=fetch_target_meta(all_tgts)
    eligible_t={t for t,x in tmeta.items() if str(x.get('target_type','')).upper()=='SINGLE PROTEIN' and str(x.get('organism','')).lower()=='homo sapiens'}
    acts=fetch_activities(all_mols,primary)
    acts=[a for a in acts if str(a.get('target_chembl_id')) in eligible_t and valid_data_comment(a) and a.get('pchembl_value') not in (None,'')]
    assays=fetch_assays({str(a.get('assay_chembl_id')) for a in acts if a.get('assay_chembl_id')})
    good_assays={aid for aid,x in assays.items() if str(x.get('assay_type','')).upper()=='B' and int(x.get('confidence_score') or -1)==9}
    revd=defaultdict(set);revg=defaultdict(set)
    for d,ms in drug_map.items():
        for m in ms: revd[m].add(d)
    for g,ts in gene_targets.items():
        for t in ts:
            if t in eligible_t: revg[t].add(g)
    blocked=set(edges);support6=set();support7=set();exact6=set();records=[]
    for a in acts:
        if str(a.get('assay_chembl_id')) not in good_assays: continue
        m=str(a.get('molecule_chembl_id'));t=str(a.get('target_chembl_id'))
        try: pv=float(a.get('pchembl_value'))
        except Exception: continue
        relation=str(a.get('standard_relation') or '').strip()
        for d in revd.get(m,[]):
            for g in revg.get(t,[]):
                e=(d,g)
                if e in blocked: continue
                support6.add(e)
                if pv>=7.0: support7.add(e)
                if relation in ('=','<','<=','≤') and pv>=6.0: exact6.add(e)
        records.append({'molecule_chembl_id':m,'target_chembl_id':t,'assay_chembl_id':a.get('assay_chembl_id'),'pchembl_value':pv,'standard_relation':relation,'standard_type':a.get('standard_type')})
    meta={'n_mapped_molecules':len(all_mols),'n_mapped_targets_raw':len(all_tgts),'n_eligible_single_protein_human_targets':len(eligible_t),'n_activities_pchembl_ge_6_on_mapped_targets':len(acts),'n_good_binding_conf9_assays':len(good_assays),'n_supported_candidate_pairs_pchembl_ge_6':len(support6),'n_supported_candidate_pairs_pchembl_ge_7':len(support7),'n_supported_candidate_pairs_exact_relation_ge_6':len(exact6)}
    return support6,support7,exact6,meta,pd.DataFrame(records)


def score_mf_matrix(m,L,R):
    li,ri,U,V,bu,bv=m
    return U@V.T+bu[:,None]+bv[None,:]

def score_svd_matrix(m,L,R):
    li,ri,U,V=m
    return U@V.T

def score_lg_matrix(m,L,R):
    li,ri,z=m
    A=np.vstack([z[li[x]] for x in L]);B=np.vstack([z[ri[y]] for y in R]);return A@B.T


def fit_ensembles(edges,ninit=5):
    L=sorted({a for a,_ in edges});R=sorted({b for _,b in edges});blocked=set(edges)
    mf=[];svd=[];gcn=[]
    for i in range(ninit):
        mf.append(score_mf_matrix(dti.fit_logistic_mf(edges,L,R,blocked,np.random.default_rng(81000+i)),L,R))
        svd.append(score_svd_matrix(dti.fit_svd(edges,L,R,32,82000+i),L,R))
        gcn.append(score_lg_matrix(lgcn.fit_lightgcn(edges,L,R,blocked,83000+i),L,R))
    return L,R,{'NeuralMF':np.mean(mf,axis=0),'SVD':np.mean(svd,axis=0),'LightGCN':np.mean(gcn,axis=0)}


def ranked_indices(score,mask):
    flat=score.ravel();idx=np.flatnonzero(mask.ravel());return idx[np.argsort(flat[idx])[::-1]]


def metrics_for_rank(name,ranked,support_flat,eligible_n):
    support_flat=set(support_flat);tot=len(support_flat);rows=[]
    base=tot/eligible_n if eligible_n else float('nan')
    for k in KS:
        kk=min(k,len(ranked));hits=sum(int(x in support_flat) for x in ranked[:kk]);prec=hits/kk if kk else float('nan');rec=hits/tot if tot else float('nan');enr=prec/base if base>0 else float('nan')
        rows.append({'model':name,'k':kk,'hits':hits,'precision_lower_bound':prec,'recall_supported':rec,'enrichment_uniform':enr})
    return rows


def flat_support(S,li,ri,nr):
    return {li[a]*nr+ri[b] for a,b in S if a in li and b in ri}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--biosnap',required=True);ap.add_argument('--out-prefix',required=True);ap.add_argument('--inits',type=int,default=5);args=ap.parse_args()
    edges=dti.parse(args.biosnap);drugs=sorted({a for a,_ in edges});genes=sorted({b for _,b in edges})
    status=get_json('https://www.ebi.ac.uk/chembl/api/data/status.json')
    version=status.get('chembl_db_version') or status.get('version')
    if version!='ChEMBL_37': raise RuntimeError(f'Protocol requires ChEMBL_37, got {version!r}')
    drug_map=map_pubchem_to_chembl(drugs);gene_acc,gene_targets=target_maps(genes)
    S6,S7,Sexact,emeta,raw=evidence_sets(edges,drug_map,gene_targets,6.0)
    L,R,scores=fit_ensembles(edges,args.inits);li={x:i for i,x in enumerate(L)};ri={x:i for i,x in enumerate(R)};blocked=set(edges)
    mask=np.ones((len(L),len(R)),dtype=bool)
    for a,b in blocked:mask[li[a],ri[b]]=False
    mapped_drugs={d for d in L if drug_map.get(d)};mapped_genes={g for g in R if gene_targets.get(g)}
    mapped_mask=np.zeros_like(mask)
    for d in mapped_drugs:
        mapped_mask[li[d],:]=True
    gene_cols=[ri[g] for g in mapped_genes]
    colmask=np.zeros(len(R),dtype=bool);colmask[gene_cols]=True;mapped_mask &= colmask[None,:];mapped_mask &= mask
    sf6=flat_support(S6,li,ri,len(R));sf7=flat_support(S7,li,ri,len(R));sfexact=flat_support(Sexact,li,ri,len(R))
    all_rows=[];rankings={}
    for name,sc in scores.items():
        r=ranked_indices(sc,mask);rankings[name]=r
        all_rows+=metrics_for_rank(name,r,sf6,int(mask.sum()))
    mapped_rows=[]
    for name,sc in scores.items():
        r=ranked_indices(sc,mapped_mask);mapped_rows+=metrics_for_rank(name,r,sf6,int(mapped_mask.sum()))
    sens=[]
    for label,S in [('pchembl_ge_7',sf7),('exact_relation_pchembl_ge_6',sfexact)]:
        for name,r in rankings.items():
            for x in metrics_for_rank(name,r,S,int(mask.sum())):x['sensitivity']=label;sens.append(x)
    p=Path(args.out_prefix);p.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(all_rows).to_csv(str(p)+'_primary_table.csv',index=False);pd.DataFrame(mapped_rows).to_csv(str(p)+'_mapped_table.csv',index=False);pd.DataFrame(sens).to_csv(str(p)+'_sensitivity_table.csv',index=False);raw.to_csv(str(p)+'_eligible_activity_records.csv',index=False)
    maps={'drug_mapping_coverage':len(mapped_drugs)/len(drugs),'gene_mapping_coverage':len(mapped_genes)/len(genes),'mapped_drugs':len(mapped_drugs),'total_drugs':len(drugs),'mapped_genes':len(mapped_genes),'total_genes':len(genes),'candidate_universe':int(mask.sum()),'mapped_candidate_universe':int(mapped_mask.sum()),**emeta}
    Path(str(p)+'_mapping_summary.json').write_text(json.dumps(maps,indent=2)+'\n')
    pri=pd.DataFrame(all_rows);comp={}
    for k in KS:
        z=pri[pri.k==min(k,int(mask.sum()))].set_index('model')
        if 'SVD' in z.index and 'NeuralMF' in z.index:comp[str(k)]={'svd_hits':int(z.loc['SVD','hits']),'neuralmf_hits':int(z.loc['NeuralMF','hits']),'difference':int(z.loc['SVD','hits']-z.loc['NeuralMF','hits'])}
    summary={'chembl_version':version,'protocol':'EXTERNAL_DTI_CHEMBL_PROTOCOL.md','initializations_per_model':args.inits,'mapping':maps,'confirmatory_SVD_minus_NeuralMF_hits':comp}
    Path(str(p)+'_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# Independent DTI external validation — ChEMBL 37','',f"ChEMBL release: **{version}**. Models were fitted on BioSNAP only; ChEMBL was opened only after rankings were defined according to the frozen protocol.",'',f"Mapping coverage: **{len(mapped_drugs)}/{len(drugs)} drugs ({len(mapped_drugs)/len(drugs):.1%})** and **{len(mapped_genes)}/{len(genes)} genes ({len(mapped_genes)/len(genes):.1%})**.",f"Primary externally supported unknown BioSNAP candidate pairs (binding, confidence 9, human single protein, pChEMBL >= 6): **{len(S6):,}**.",'','| Model | K | ChEMBL-supported hits | Precision lower bound | Recall of supported pairs | Enrichment vs uniform |','|---|---:|---:|---:|---:|---:|']
    for _,x in pri.iterrows():lines.append(f"| {x.model} | {int(x.k):,} | {int(x.hits):,} | {x.precision_lower_bound:.4f} | {x.recall_supported:.4f} | {x.enrichment_uniform:.1f}× |")
    lines+=['','Primary confirmatory contrast was frozen before outcome inspection: SVD (structure-neutralized-selected) minus NeuralMF (conventional-selected) supported-hit yield at prespecified K. Unlabeled ChEMBL pairs are not called negatives.']
    Path(str(p)+'.md').write_text('\n'.join(lines)+'\n');print('\n'.join(lines))

if __name__=='__main__':main()
