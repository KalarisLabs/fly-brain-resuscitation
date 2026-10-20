"""Load all mechanism-track runs (main + follow-up + sign arms) into one DataFrame with FII_norm.
exec(open(r'$RESEARCH_ROOT\\followup\\mechanism\\code\\analyze_load.py').read())  -> df, Rs, fii_matrix, boot_mean
"""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src", _os.path.join("experiments", "dynamics")):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import os, glob, json
import numpy as np, pandas as pd
OUTD = os.path.join(_P.OUT, 'followup', 'mechanism'); RAW = os.path.join(OUTD, 'raw'); DYN = os.path.join(_P.OUT, 'dynamics')
CEIL = 0.9940
rng = np.random.default_rng(0)


def fii_matrix(Rb, Rc):
    U = (Rb > 0) | (Rc > 0)
    if U.sum() < 3 or not (Rc[U] != 0).any():
        return 0.0
    x, y = Rb[U].ravel(), Rc[U].ravel()
    if x.std() == 0 or y.std() == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def boot_mean(x, n=2000):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) == 0:
        return (np.nan, np.nan, np.nan)
    b = rng.choice(x, (n, len(x))).mean(1)
    return (x.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5))


Rbase = np.mean([np.load(f)['R'].astype(np.float32) for f in sorted(glob.glob(os.path.join(DYN, 'raw', 'base', '*.npz')))], 0)
rows = []; Rs = {}
for f in glob.glob(os.path.join(RAW, 'results*', '*', '*.json')):
    o = json.load(open(f))
    if o['id'] in Rs or any(r['id'] == o['id'] for r in rows[-1:]):
        pass
    d = {k: v for k, v in o.items() if not isinstance(v, (dict, list))}
    if 'lhs' in o:
        d.update({'p_' + k: v for k, v in o['lhs'].items()})
    if o['arm'] not in ('N', 'Y', 'SN', 'SY'):
        R = np.load(f[:-5] + '.npz')['R'].astype(np.float32); Rs[o['id']] = R
        d['FII_def'] = fii_matrix(Rbase, R)
    d['_src'] = os.path.basename(os.path.dirname(os.path.dirname(f)))
    rows.append(d)
df = pd.DataFrame(rows).drop_duplicates('id', keep='last').reset_index(drop=True)
# per-set baselines (L) and per-variant baselines (S)
df['ceil'] = np.nan; df['FII'] = np.nan
Lb = df[(df.arm == 'L') & (df.insult == 'none')]
for key, sel_base, sel_runs in [
        ('L', lambda s: (df.arm == 'L') & (df.insult == 'none') & (df.p_set == s), lambda s: df.arm.isin(['L', 'F']) & (df.p_set == s) & (df.insult != 'none')),
]:
    for s in Lb.p_set.unique():
        RR = [Rs[i] for i in df[sel_base(s)].id]
        if len(RR) < 2:
            continue
        Rb = np.mean(RR, 0)
        loo = [fii_matrix(np.mean([RR[j] for j in range(len(RR)) if j != i], 0), RR[i]) for i in range(len(RR))]
        m = sel_runs(s) | sel_base(s)
        df.loc[m, 'ceil'] = float(np.median(loo))
        df.loc[m, 'FII'] = [fii_matrix(Rb, Rs[i]) for i in df.loc[m, 'id']]
if 'signvar' in df:
    for var in ('A', 'B', 'C'):
        mb = (df.arm == 'S') & (df.signvar == var) & (df.insult == 'none')
        RR = [Rs[i] for i in df[mb].id]
        if len(RR) < 2:
            continue
        Rb = np.mean(RR, 0)
        loo = [fii_matrix(np.mean([RR[j] for j in range(len(RR)) if j != i], 0), RR[i]) for i in range(len(RR))]
        m = df.arm.isin(['S', 'ST']) & (df.signvar == var)
        df.loc[m, 'ceil'] = float(np.median(loo))
        df.loc[m, 'FII'] = [fii_matrix(Rb, Rs[i]) for i in df.loc[m, 'id']]
        df.loc[m, 'FII_vs_default'] = df.loc[m, 'FII_def']
m2 = df.arm.isin(['X', 'T', 'M'])
df.loc[m2, 'ceil'] = CEIL; df.loc[m2, 'FII'] = df.loc[m2, 'FII_def']
df['FII_norm'] = df.FII / df.ceil
df['runaway_cls'] = (df.obs_late_Hz > 0.5) | (df.FSA > 0.05)
