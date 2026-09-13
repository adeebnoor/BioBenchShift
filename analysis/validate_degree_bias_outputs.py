#!/usr/bin/env python3
"""Validate frozen ATC5 structural-bias outputs and manifest invariants."""
import hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parents[1]
r=pd.read_csv(root/'data/degree_bias_atc5_replicates.csv')
s=json.loads((root/'data/degree_bias_atc5_summary.json').read_text())
assert len(r)==20
assert r.seed.tolist()==list(range(20))
checks={
 'auc_random_mean':r.auc_random_unlabelled.mean(),
 'auc_random_sd':r.auc_random_unlabelled.std(ddof=1),
 'auc_curated_mean':r.auc_curated_atc5.mean(),
 'auc_curated_sd':r.auc_curated_atc5.std(ddof=1),
 'auc_degree_matched_mean':r.auc_degree_matched.mean(),
 'auc_degree_matched_sd':r.auc_degree_matched.std(ddof=1),
 'n_matched_mean':r.n_matched.mean(),
}
for k,v in checks.items():
    assert np.isclose(v,s[k],rtol=0,atol=1e-12),(k,v,s[k])
d=r.auc_curated_atc5-r.auc_degree_matched
lo,hi=np.quantile(d,[.025,.975])
assert np.isclose(d.mean(),s['curated_minus_matched_mean'],atol=1e-12)
assert np.allclose([lo,hi],s['curated_minus_matched_empirical_95_interval'],atol=1e-12)
assert [int(r.n_matched.min()),int(r.n_matched.max())]==s['n_matched_range']
assert s['positive_atc5_pairs']==8094 and s['curated_atc5_pairs_analyzed']==992
manifest=(root/'data/FROZEN_REFERENCE_MANIFEST.md').read_text()
assert '1f96a2c3197db78ea7fe31444efec2e4c4cfffb19a5f1fed7bd2d0ca56166018' in manifest
print('ATC5 structural-bias frozen outputs: PASS')
