#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, time
from collections import defaultdict
from pathlib import Path
import pandas as pd
import requests

S=requests.Session(); S.headers.update({'User-Agent':'ANTI-DDI-Science-GraphBAN-mapping/1.0'})

def get(url, params=None, tries=6, timeout=120):
    err=None
    for i in range(tries):
        try:
            r=S.get(url,params=params,timeout=timeout)
            if r.status_code in (429,500,502,503,504): time.sleep(min(2**i,20)); continue
            r.raise_for_status(); return r
        except Exception as e: err=e; time.sleep(min(2**i,20))
    raise RuntimeError(f'GET failed {url}: {err}')

def post(url,data,tries=6,timeout=120):
    err=None
    for i in range(tries):
        try:
            r=S.post(url,data=data,timeout=timeout)
            if r.status_code in (429,500,502,503,504): time.sleep(min(2**i,20)); continue
            r.raise_for_status(); return r
        except Exception as e: err=e; time.sleep(min(2**i,20))
    raise RuntimeError(f'POST failed {url}: {err}')

def chunks(xs,n):
    xs=list(xs)
    for i in range(0,len(xs),n): yield xs[i:i+n]

def parse(path):
    d=pd.read_csv(path,compression='infer',comment='#',header=None,names=['drug','gene'],sep=',',dtype=str).dropna()
    return sorted({(str(a).strip(),str(b).strip()) for a,b in zip(d.drug,d.gene)})

def cid_norm(x):
    s=str(x).strip()
    if s.upper().startswith('CID'): s=s[3:]
    return str(int(s))

def pubchem_smiles(drugs):
    c2d={cid_norm(d):d for d in drugs}; out={}
    for ch in chunks(c2d,80):
        url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+','.join(ch)+'/property/CanonicalSMILES,ConnectivitySMILES,IsomericSMILES/JSON'
        j=get(url).json()
        for p in j.get('PropertyTable',{}).get('Properties',[]):
            cid=str(p.get('CID'))
            sm=p.get('CanonicalSMILES') or p.get('ConnectivitySMILES') or p.get('IsomericSMILES')
            if cid in c2d and sm: out[c2d[cid]]=str(sm)
    return out

def uniprot_gene_map(genes):
    ids=sorted(set(map(str,genes)))
    job=post('https://rest.uniprot.org/idmapping/run',{'from':'GeneID','to':'UniProtKB','ids':','.join(ids)}).json()['jobId']
    for _ in range(180):
        j=get(f'https://rest.uniprot.org/idmapping/status/{job}').json()
        st=j.get('jobStatus')
        if st in ('NEW','RUNNING'): time.sleep(2); continue
        if st and st!='FINISHED': raise RuntimeError(str(j))
        break
    detail=get(f'https://rest.uniprot.org/idmapping/details/{job}').json(); url=detail['redirectURL']
    url += ('&' if '?' in url else '?')+'size=500'
    cand=defaultdict(list)
    while url:
        r=get(url); j=r.json()
        for rec in j.get('results',[]):
            gid=str(rec.get('from')); to=rec.get('to')
            if not isinstance(to,dict): continue
            org=to.get('organism') or {}
            if org.get('taxonId') not in (None,9606): continue
            acc=to.get('primaryAccession') or to.get('uniProtkbId') or to.get('id')
            if not acc: continue
            entry=str(to.get('entryType') or '')
            reviewed=('reviewed' in entry.lower() and 'unreviewed' not in entry.lower())
            seq=(to.get('sequence') or {}).get('value') if isinstance(to.get('sequence'),dict) else None
            cand[gid].append({'accession':str(acc),'reviewed':reviewed,'sequence':str(seq) if seq else ''})
        link=r.headers.get('Link',''); nxt=None
        if 'rel="next"' in link: nxt=link.split('<',1)[1].split('>',1)[0]
        url=nxt
    chosen={}; multiplicity={}
    for g,rows in cand.items():
        rows=sorted(rows,key=lambda x:(not x['reviewed'],not bool(x['sequence']),x['accession']))
        chosen[g]=rows[0]
        multiplicity[g]=len(rows)
    return chosen,multiplicity

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input',required=True); ap.add_argument('--out-prefix',required=True); a=ap.parse_args()
    E=parse(a.input); drugs=sorted({x for x,_ in E}); genes=sorted({y for _,y in E})
    sm=pubchem_smiles(drugs); gm,mult=uniprot_gene_map(genes)
    gm_seq={g:x for g,x in gm.items() if x.get('sequence')}
    mapped_edges=[e for e in E if e[0] in sm and e[1] in gm_seq]
    summary={
      'benchmark':'BioSNAP TargetDecagon DTI','edges':len(E),'drugs':len(drugs),'genes':len(genes),
      'mapped_drugs':len(sm),'mapped_drug_fraction':len(sm)/len(drugs),
      'mapped_genes':len(gm),'mapped_gene_fraction':len(gm)/len(genes),
      'mapped_genes_with_sequence':len(gm_seq),'mapped_gene_with_sequence_fraction':len(gm_seq)/len(genes),
      'mapped_positive_edges':len(mapped_edges),'mapped_positive_edge_fraction':len(mapped_edges)/len(E),
      'genes_with_multiple_human_uniprot_mappings':sum(v>1 for v in mult.values()),
      'selected_reviewed_genes':sum(bool(x['reviewed']) for x in gm.values()),
      'mapping_gate_threshold':0.90,
      'mapping_gate_pass':len(mapped_edges)/len(E)>=0.90,
      'note':'Mapping-only preflight; positive-edge gate requires both SMILES and an actual human UniProt sequence. No model outcomes inspected.'
    }
    p=Path(a.out_prefix); p.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame([{'drug':d,'smiles':sm[d]} for d in sorted(sm)]).to_csv(str(p)+'_drug_map.csv',index=False)
    pd.DataFrame([{'gene':g,**gm[g],'n_human_mappings':mult.get(g,0)} for g in sorted(gm)]).to_csv(str(p)+'_gene_map.csv',index=False)
    Path(str(p)+'_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    md=f"""# GraphBAN TargetDecagon mapping gate\n\n- Positive edges: **{len(E):,}**\n- Drug mapping: **{len(sm)}/{len(drugs)} ({len(sm)/len(drugs):.1%})**\n- Human UniProt mapping: **{len(gm)}/{len(genes)} ({len(gm)/len(genes):.1%})**\n- Human UniProt mapping with sequence: **{len(gm_seq)}/{len(genes)} ({len(gm_seq)/len(genes):.1%})**\n- Positive-edge feature mapping: **{len(mapped_edges)}/{len(E)} ({len(mapped_edges)/len(E):.1%})**\n- Predeclared >=90% edge-coverage gate: **{'PASS' if summary['mapping_gate_pass'] else 'FAIL'}**\n\nThis is a mapping-only preflight. It contains no GraphBAN performance result.\n"""
    Path(str(p)+'.md').write_text(md); print(md)
if __name__=='__main__': main()
