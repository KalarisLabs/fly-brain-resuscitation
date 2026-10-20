"""Sign-robustness analysis (local reduced arm).  exec() in a kernel; writes tables/S_*.csv and returns SG, SYT, CORE."""
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
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
OUTD = os.path.join(_P.OUT, 'followup', 'mechanism'); RES = os.path.join(OUTD, 'raw', 'sign_res'); CA = _P.CACHE
DYN = os.path.join(_P.OUT, 'dynamics')
rng = np.random.default_rng(0)


def fii_matrix(Rb, Rc):
    U = (Rb > 0) | (Rc > 0)
    if U.sum() < 3 or not (Rc[U] != 0).any():
        return 0.0
    x, y = Rb[U].ravel(), Rc[U].ravel()
    if x.std() == 0 or y.std() == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def boot(x, n=2000):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) == 0:
        return (np.nan, np.nan, np.nan)
    b = rng.choice(x, (n, len(x))).mean(1); return (x.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5))


rows = []; Rs = {}
for f in glob.glob(os.path.join(RES, '*', '*.json')):
    o = json.load(open(f)); d = {k: v for k, v in o.items() if not isinstance(v, (dict, list))}
    if o['arm'] in ('S', 'ST'):
        z = np.load(f[:-5] + '.npz'); Rs[o['id']] = z['R'].astype(np.float32)
    rows.append(d)
SD = pd.DataFrame(rows)
Rdef = np.mean([np.load(f)['R'].astype(np.float32) for f in sorted(glob.glob(os.path.join(DYN, 'raw', 'base', '*.npz')))], 0)
SD['FII'] = np.nan; SD['ceil'] = np.nan; SD['FII_vs_default_base'] = np.nan
for v in ('A', 'B'):
    mb = (SD.arm == 'S') & (SD.signvar == v) & (SD.insult == 'none')
    RR = [Rs[i] for i in SD[mb].id]
    if len(RR) < 2:
        continue
    Rb = np.mean(RR, 0)
    loo = [fii_matrix(np.mean([RR[j] for j in range(len(RR)) if j != i], 0), RR[i]) for i in range(len(RR))]
    m = SD.arm.isin(['S', 'ST']) & (SD.signvar == v)
    SD.loc[m, 'ceil'] = float(np.median(loo))
    SD.loc[m, 'FII'] = [fii_matrix(Rb, Rs[i]) for i in SD.loc[m, 'id']]
    SD.loc[m, 'FII_vs_default_base'] = [fii_matrix(Rdef, Rs[i]) for i in SD.loc[m, 'id']]
    SD.loc[mb, 'base_vs_default_FII'] = [fii_matrix(Rdef, Rs[i]) for i in SD.loc[mb, 'id']]
SD['FII_norm'] = SD.FII / SD.ceil
SD['runaway_cls'] = (SD.obs_late_Hz > 0.5) | (SD.FSA > 0.05)
SD.to_csv(os.path.join(OUTD, 'tables', 'S_all_runs.csv'), index=False)
