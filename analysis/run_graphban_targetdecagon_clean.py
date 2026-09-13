#!/usr/bin/env python3
"""Leakage-free GraphBAN-style contemporary challenge on frozen TargetDecagon DTI.

Scientific rules are frozen in CONTEMPORARY_MODEL_CHALLENGE.md and
GRAPHBAN_PROTOCOL_AMENDMENT_20260913.md. Held-out positive edges never enter
message passing. The same trained model is evaluated against conventional
random-unlabelled and training-degree-matched controls.
"""
from __future__ import annotations
import argparse, json, math, random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import StandardScaler
from torch import nn
from torch.nn import Linear, Dropout
from torch_geometric.data import HeteroData
from torch_geometric.nn import SAGEConv, to_hetero
import torch_geometric.transforms as T


def parse_edges(path):
    d=pd.read_csv(path,compression='infer',comment='#',header=None,names=['drug','gene'],sep=',',dtype=str).dropna()
    return sorted({(str(a).strip(),str(b).strip()) for a,b in zip(d.drug,d.gene)})

def dbin(d): return 0 if d<=0 else int(math.floor(math.log2(d)))+1

def sample_random(rng,n,L,R,blocked):
    out=set()
    while len(out)<n:
        e=(L[int(rng.integers(len(L)))],R[int(rng.integers(len(R)))])
        if e not in blocked: out.add(e)
    return list(out)

def degree_matched(rng,pos,L,R,blocked,dl,dr):
    lb=defaultdict(list); rb=defaultdict(list)
    for x in L: lb[dbin(dl.get(x,0))].append(x)
    for x in R: rb[dbin(dr.get(x,0))].append(x)
    used=set(); P=[]; N=[]
    for ii in rng.permutation(len(pos)):
        p=pos[int(ii)]; A=lb.get(dbin(dl[p[0]]),[]); B=rb.get(dbin(dr[p[1]]),[])
        z=None
        for _ in range(1000):
            if not A or not B: break
            e=(A[int(rng.integers(len(A)))],B[int(rng.integers(len(B)))])
            if e not in blocked and e not in used: z=e; break
        if z is not None: used.add(z); P.append(p); N.append(z)
    return P,N

class GNNEncoder(torch.nn.Module):
    def __init__(self, hidden_channels, out_channels):
        super().__init__()
        self.conv1=SAGEConv((-1,-1),out_channels)
        self.conv2=SAGEConv((-1,-1),out_channels)
        self.conv3=SAGEConv((-1,-1),out_channels)
    def forward(self,x,edge_index):
        x=self.conv1(x,edge_index).relu()
        x=self.conv2(x,edge_index).relu()
        x=self.conv3(x,edge_index).relu()
        return x

class EdgeDecoder(torch.nn.Module):
    def __init__(self,hidden_channels):
        super().__init__()
        self.lin1=Linear(2*hidden_channels,hidden_channels)
        self.lin2=Linear(hidden_channels,hidden_channels)
        self.lin3=Linear(hidden_channels,hidden_channels)
        self.lin4=Linear(hidden_channels,1)
        self.dropout=Dropout(p=0.2)
    def forward(self,z_dict,edge_label_index):
        row,col=edge_label_index
        x=torch.cat([z_dict['chem'][row],z_dict['protein'][col]],dim=-1)
        z=self.dropout(self.lin1(x).relu())
        z=self.dropout(self.lin2(z).relu())
        z=self.dropout(self.lin3(z))
        return self.lin4(z).view(-1)

class GraphBANTransductive(torch.nn.Module):
    def __init__(self,data,hidden_channels=256):
        super().__init__()
        self.encoder_base=GNNEncoder(hidden_channels,hidden_channels)
        self.encoder=to_hetero(self.encoder_base,data.metadata(),aggr='sum')
        self.decoder=EdgeDecoder(hidden_channels)
    def forward(self,x_dict,edge_index_dict,edge_label_index):
        z=self.encoder(x_dict,edge_index_dict)
        return self.decoder(z,edge_label_index)

def extract_features(drug_df,gene_df,device):
    from transformers import AutoTokenizer, RobertaModel
    import esm

    torch.set_num_threads(max(1,min(4,torch.get_num_threads())))
    # ChemBERTa feature definition follows the public GraphBAN transductive code.
    tok=AutoTokenizer.from_pretrained('DeepChem/ChemBERTa-77M-MTR')
    chem=RobertaModel.from_pretrained('DeepChem/ChemBERTa-77M-MTR',add_pooling_layer=True).eval().to(device)
    smiles=drug_df.smiles.astype(str).tolist(); cfeat=[]
    for s0 in range(0,len(smiles),16):
        batch=smiles[s0:s0+16]
        enc=tok(batch,return_tensors='pt',padding='max_length',max_length=290,truncation=True).to(device)
        with torch.no_grad(): out=chem(**enc).last_hidden_state[:,0,:]
        cfeat.append(out.cpu().numpy().astype('float32'))
    cfeat=np.vstack(cfeat)
    del chem

    esm_model,alphabet=esm.pretrained.esm1b_t33_650M_UR50S(); esm_model=esm_model.eval().to(device)
    converter=alphabet.get_batch_converter(); seqs=gene_df.sequence.fillna('').astype(str).tolist(); pfeat=[]
    # Preserve row order; use small batches as in the public implementation.
    for s0 in range(0,len(seqs),4):
        b=seqs[s0:s0+4]; dat=[(f'p{s0+i}',seq[:1022]) for i,seq in enumerate(b)]
        _,_,tokens=converter(dat); tokens=tokens.to(device)
        with torch.no_grad(): rep=esm_model(tokens,repr_layers=[33],return_contacts=False)['representations'][33]
        for i,(_,seq) in enumerate(dat):
            pfeat.append(rep[i,1:len(seq)+1].mean(0).cpu().numpy().astype('float32'))
    pfeat=np.vstack(pfeat)
    return cfeat,pfeat

def make_graph(L,R,train,Xd,Xp,device):
    li={x:i for i,x in enumerate(L)}; ri={x:i for i,x in enumerate(R)}
    e=torch.tensor([[li[a] for a,b in train],[ri[b] for a,b in train]],dtype=torch.long)
    data=HeteroData(); data['chem'].x=torch.tensor(Xd,dtype=torch.float32); data['protein'].x=torch.tensor(Xp,dtype=torch.float32)
    data['chem','CPI','protein'].edge_index=e
    data=T.ToUndirected()(data)
    return data.to(device),li,ri

def edge_tensor(E,li,ri,device):
    return torch.tensor([[li[a] for a,b in E],[ri[b] for a,b in E]],dtype=torch.long,device=device)

def model_scores(model,data,E,li,ri,device,batch=8192):
    model.eval(); out=[]
    with torch.no_grad():
        for s in range(0,len(E),batch):
            ix=edge_tensor(E[s:s+batch],li,ri,device)
            out.append(torch.sigmoid(model(data.x_dict,data.edge_index_dict,ix)).cpu().numpy())
    return np.concatenate(out) if out else np.array([],dtype=float)

def metrics(pos,neg,model,data,li,ri,device):
    sp=model_scores(model,data,pos,li,ri,device); sn=model_scores(model,data,neg,li,ri,device)
    y=np.r_[np.ones(len(sp)),np.zeros(len(sn))]; s=np.r_[sp,sn]
    return float(roc_auc_score(y,s)),float(average_precision_score(y,s))

def train_model(data,li,ri,train,blocked,L,R,seed,device,epochs=20):
    torch.manual_seed(seed); np.random.seed(seed); random.seed(seed)
    model=GraphBANTransductive(data,256).to(device)
    # Initialize lazy heterogeneous modules before constructing optimizer.
    with torch.no_grad(): model.encoder(data.x_dict,data.edge_index_dict)
    opt=torch.optim.Adam(model.parameters(),lr=0.001); lossfn=nn.BCEWithLogitsLoss(); rng=np.random.default_rng(seed+90000)
    pos=list(train)
    for _ep in range(epochs):
        model.train(); neg=sample_random(rng,len(pos),L,R,blocked)
        order=rng.permutation(len(pos)); bs=2048
        for s0 in range(0,len(pos),bs):
            pp=[pos[int(i)] for i in order[s0:s0+bs]]; nn0=neg[s0:s0+len(pp)]
            E=pp+nn0; y=torch.cat([torch.ones(len(pp)),torch.zeros(len(nn0))]).to(device)
            opt.zero_grad(); logits=model(data.x_dict,data.edge_index_dict,edge_tensor(E,li,ri,device)); loss=lossfn(logits,y); loss.backward(); opt.step()
    return model

def one_seed(E_mapped,E_all,L,R,baseXd,baseXp,drug_index,gene_index,seed,device):
    # Preserve project split convention: split full frozen edge list first, then apply mapping.
    rng=np.random.default_rng(seed); ix=rng.permutation(len(E_all)); nt=round(.2*len(E_all)); test_all={E_all[i] for i in ix[:nt]}; train_all={E_all[i] for i in ix[nt:]}
    mapped=set(E_mapped); train=sorted(train_all & mapped); test=sorted(test_all & mapped); blocked=set(E_mapped)
    dl=Counter(a for a,b in train); dr=Counter(b for a,b in train); seen=[e for e in test if dl[e[0]]>0 and dr[e[1]]>0]
    # Fit feature scaling on unique training endpoints only, then apply to all mapped nodes.
    td={a for a,b in train}; tg={b for a,b in train}
    dsc=StandardScaler().fit(baseXd[[drug_index[x] for x in sorted(td)]])
    psc=StandardScaler().fit(baseXp[[gene_index[x] for x in sorted(tg)]])
    Xd=dsc.transform(baseXd).astype('float32'); Xp=psc.transform(baseXp).astype('float32')
    data,li,ri=make_graph(L,R,train,Xd,Xp,device)
    model=train_model(data,li,ri,train,blocked,L,R,seed,device,20)
    erng=np.random.default_rng(seed+40000); rn=sample_random(erng,len(seen),L,R,blocked); mp,mn=degree_matched(erng,seen,L,R,blocked,dl,dr)
    ar,pr=metrics(seen,rn,model,data,li,ri,device); am,pm=metrics(mp,mn,model,data,li,ri,device)
    return {'seed':seed,'n_train':len(train),'n_test_mapped':len(test),'n_seen':len(seen),'n_matched':len(mp),'match_fraction':len(mp)/len(seen) if seen else 0,
            'auc_random':ar,'auprc_random':pr,'auc_matched':am,'auprc_matched':pm,'auc_drop':ar-am,'auprc_drop':pr-pm}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input',required=True); ap.add_argument('--drug-map',required=True); ap.add_argument('--gene-map',required=True); ap.add_argument('--out-prefix',required=True); ap.add_argument('--seeds',default='0,1,2'); a=ap.parse_args()
    E=parse_edges(a.input); dd=pd.read_csv(a.drug_map,dtype=str); gd=pd.read_csv(a.gene_map,dtype=str).fillna('')
    dd=dd[dd.smiles.astype(bool)].drop_duplicates('drug'); gd=gd[gd.sequence.astype(bool)].drop_duplicates('gene')
    dm=dict(zip(dd.drug,dd.smiles)); gm=dict(zip(gd.gene,gd.sequence)); mapped=[e for e in E if e[0] in dm and e[1] in gm]
    coverage=len(mapped)/len(E); 
    if coverage<0.90: raise RuntimeError(f'predeclared GraphBAN mapping gate failed: {coverage:.4f}')
    L=sorted({a0 for a0,b0 in mapped}); R=sorted({b0 for a0,b0 in mapped})
    dd=dd.set_index('drug').loc[L].reset_index(); gd=gd.set_index('gene').loc[R].reset_index()
    device=torch.device('cuda' if torch.cuda.is_available() else 'cpu'); print('device',device)
    Xd,Xp=extract_features(dd,gd,device); di={x:i for i,x in enumerate(L)}; gi={x:i for i,x in enumerate(R)}
    rows=[one_seed(mapped,E,L,R,Xd,Xp,di,gi,int(s),device) for s in a.seeds.split(',') if s.strip()]
    out=pd.DataFrame(rows); p=Path(a.out_prefix); p.parent.mkdir(parents=True,exist_ok=True); out.to_csv(str(p)+'_replicates.csv',index=False)
    summary={'benchmark':'BioSNAP TargetDecagon DTI','architecture':'leakage-free GraphBAN-style transductive GraphSAGE + ChemBERTa + ESM-1b','mapping_edge_fraction':coverage,'seeds':[int(x) for x in out.seed],
             'auc_random_mean':float(out.auc_random.mean()),'auc_random_sd':float(out.auc_random.std(ddof=1)),'auc_matched_mean':float(out.auc_matched.mean()),'auc_matched_sd':float(out.auc_matched.std(ddof=1)),
             'auc_drop_mean':float(out.auc_drop.mean()),'auprc_random_mean':float(out.auprc_random.mean()),'auprc_matched_mean':float(out.auprc_matched.mean()),'match_fraction_mean':float(out.match_fraction.mean()),
             'heldout_edges_in_message_passing':False,'epochs':20,'hidden_dim':256}
    Path(str(p)+'_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    md=f"""# Leakage-free contemporary GraphBAN challenge — TargetDecagon DTI\n\nMapping edge coverage: **{coverage:.1%}**. Held-out positives in message passing: **no**.\n\n| Metric | Conventional | Structure-neutralized | Difference |\n|---|---:|---:|---:|\n| AUROC | {summary['auc_random_mean']:.3f} ± {summary['auc_random_sd']:.3f} | {summary['auc_matched_mean']:.3f} ± {summary['auc_matched_sd']:.3f} | {summary['auc_drop_mean']:+.3f} conventional-minus-neutralized |\n| AUPRC | {summary['auprc_random_mean']:.3f} | {summary['auprc_matched_mean']:.3f} | {float(out.auprc_drop.mean()):+.3f} |\n\nMean matching coverage: **{summary['match_fraction_mean']:.3f}** across seeds {summary['seeds']}.\n\nThis is the primary S6 contemporary-model challenge defined before outcome inspection. The public GraphBAN transductive test-graph construction is not used here because held-out positive edges must not enter message passing.\n"""
    Path(str(p)+'.md').write_text(md); print(md)
if __name__=='__main__': main()
