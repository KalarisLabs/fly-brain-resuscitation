"""Compile keystone tables, write per-dataset 'extra types' lists (union of other datasets' top-50 names present
in this dataset but not yet evaluated), and (mode=final) build the replication table.
python keystone_replication.py prep|final
Ranking: primary = rel_reachN_L6 (relative drop in strict sensory->sink neuron-pair reach, L<=6);
tie-break = rel_influence_L6.  Name matching: exact, case-insensitive cell-type name.
"""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src",):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import os, sys, re
import numpy as np, pandas as pd
from scipy.stats import hypergeom
ROOT = os.path.join(_P.OUT, 'structure'); T = os.path.join(ROOT, 'tables')
DSS = ['FAFB', 'BANC', 'MCNS', 'MANC']
TOP = 50
mode = sys.argv[1]
K = {}
for ds in DSS:
    p = os.path.join(T, f'{ds}_keystone.csv')
    if os.path.exists(p):
        d = pd.read_csv(p).drop_duplicates('type')
        d['key'] = d.type.astype(str).str.lower()
        d = d.sort_values(['rel_reachN_L6', 'rel_influence_L6'], ascending=False).reset_index(drop=True)
        d['rank'] = np.arange(1, len(d) + 1)
        K[ds] = d
top = {ds: d[d['rank'] <= TOP] for ds, d in K.items()}
types_all = {}
for ds in DSS:
    fc = os.path.join(T, f'{ds}_type_flow_centrality.csv')
    if os.path.exists(fc):
        types_all[ds] = pd.read_csv(fc)
if mode == 'prep':
    union = set().union(*[set(t.type.astype(str)) for t in top.values()])
    for ds, tf in types_all.items():
        lower = {str(t).lower(): str(t) for t in tf.type}
        done = set(K[ds].key) if ds in K else set()
        extra = sorted({lower[u.lower()] for u in union if u.lower() in lower and u.lower() not in done})
        open(os.path.join(T, f'{ds}_extra_types.txt'), 'w', encoding='utf8').write('\n'.join(extra))
        print(ds, 'extra types to evaluate:', len(extra))
else:
    rows = []
    union = sorted(set().union(*[set(t.key) for t in top.values()]))
    for u in union:
        r = {'type_key': u}
        for ds, d in K.items():
            hit = d[d.key == u]
            if len(hit):
                h = hit.iloc[0]
                r[f'{ds}_rank'] = int(h['rank']); r[f'{ds}_relReachL6'] = h.rel_reachN_L6
                r[f'{ds}_relInfl'] = h.rel_influence_L6; r[f'{ds}_n'] = int(h.n_neurons); r[f'{ds}_top50'] = bool(h['rank'] <= TOP)
            else:
                present = ds in types_all and u in set(types_all[ds].type.astype(str).str.lower())
                r[f'{ds}_rank'] = np.nan; r[f'{ds}_top50'] = False; r[f'{ds}_present'] = present
        r['n_top50_brain'] = sum(r.get(f'{ds}_top50', False) for ds in ['FAFB', 'BANC', 'MCNS'])
        r['n_top50_vnc'] = sum(r.get(f'{ds}_top50', False) for ds in ['MANC', 'BANC', 'MCNS'])
        r['n_top50_all'] = sum(r.get(f'{ds}_top50', False) for ds in DSS)
        rows.append(r)
    R = pd.DataFrame(rows).sort_values(['n_top50_all', 'FAFB_relReachL6'], ascending=[False, False])
    R.to_csv(os.path.join(T, 'keystone_replication.csv'), index=False)
    # pairwise overlap of top-50 with hypergeometric p (universe = types present in both datasets' typed removable pool)
    ov = []
    for i, a in enumerate(DSS):
        for b in DSS[i + 1:]:
            if a not in K or b not in K: continue
            ua = set(types_all[a].type.astype(str).str.lower()); ub = set(types_all[b].type.astype(str).str.lower())
            U = ua & ub
            ta = set(top[a].key) & U; tb = set(top[b].key) & U
            k = len(ta & tb)
            p = hypergeom.sf(k - 1, len(U), len(ta), len(tb)) if len(U) else np.nan
            ov.append(dict(a=a, b=b, shared_type_universe=len(U), top50_a_in_universe=len(ta), top50_b_in_universe=len(tb),
                           overlap=k, expected=len(ta) * len(tb) / max(len(U), 1), p_hypergeom=p,
                           shared=';'.join(sorted(ta & tb))))
    pd.DataFrame(ov).to_csv(os.path.join(T, 'keystone_overlap.csv'), index=False)
    for ds, d in K.items():
        d.head(TOP).to_csv(os.path.join(T, f'{ds}_keystone_top50.csv'), index=False)
    print(pd.DataFrame(ov)[['a', 'b', 'shared_type_universe', 'overlap', 'expected', 'p_hypergeom']].to_string())
    print(R.head(30)[['type_key', 'n_top50_all', 'n_top50_brain', 'n_top50_vnc'] + [c for c in R.columns if c.endswith('_rank')]].to_string())

