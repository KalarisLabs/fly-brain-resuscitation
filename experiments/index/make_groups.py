"""Neuron group masks (FAFB princeton5 order) used to attribute runaway spiking (Modal job spont_attrib.py)."""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src",):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import os, sys, math, numpy as np, scipy.sparse as sp
D = os.environ.get('RESEARCH_ROOT', _P.OUT); OUT = os.path.join(D, 'followup', 'index')
sys.path.insert(0, OUT); from light import load_light, LIF_DEFAULTS
M = load_light(); N = M.N
L = LIF_DEFAULTS; dV = L['v_th'] - L['v_rest']; GMF = 0.185 * L['tau_syn'] / (L['tau_m'] * dV)
rows = np.repeat(np.arange(N, dtype=np.int32), np.diff(M.indptr)); cols = M.indices; w = M.weight
exc = w > 0
WeT = sp.csr_matrix((w[exc], (cols[exc], rows[exc])), shape=(N, N))
WiT = sp.csr_matrix((-w[~exc], (cols[~exc], rows[~exc])), shape=(N, N))
# Perron vector
x = np.ones(N, np.float32) / math.sqrt(N)
for k in range(400):
    y = WeT @ x + x; x = y / np.linalg.norm(y)
g = {}
g['sign_uncertain'] = ~M.fast
order = np.argsort(-x); g['perron_top100'] = np.zeros(N, bool); g['perron_top100'][order[:100]] = True
g['perron_top1000'] = np.zeros(N, bool); g['perron_top1000'][order[:1000]] = True
# supercritical 2-loops (reciprocal exc pairs with loop gain > 1)
key = rows[exc].astype(np.int64) * N + cols[exc]; ks = np.sort(key); wk = w[exc][np.argsort(key)]
rk = cols[exc].astype(np.int64) * N + rows[exc]; p = np.searchsorted(ks, rk); p[p >= ks.size] = 0
ok = ks[p] == rk; gain2 = np.where(ok, w[exc] * wk[p] * GMF ** 2, 0)
sel = gain2 > 1; g['loop2_super'] = np.zeros(N, bool); g['loop2_super'][rows[exc][sel]] = True; g['loop2_super'][cols[exc][sel]] = True
for r in (100, 200):
    th = dV / (0.185 * L['tau_syn'] * r / 1000.0); S = np.ones(N, bool)
    for _ in range(500):
        S2 = S & (WeT @ S.astype(np.float32) >= th)
        if S2.sum() == S.sum(): break
        S = S2
    g[f'ecore{r}'] = S
deg = np.diff(M.indptr) + np.bincount(M.indices, minlength=N); o = np.argsort(-deg, kind='stable')
g['hub1'] = np.zeros(N, bool); g['hub1'][o[:int(round(0.01 * N))]] = True
g['inhibitory'] = M.sign < 0
sc = M.neurons.super_class.fillna('NA').values; names = np.unique(sc)
for s in names: g['sc_' + s] = sc == s
al = M.neurons.primary_type.fillna('').str.contains(r'^l?v?\d?LN|LN\d', regex=True).values
g['AL_LN_like'] = al
np.savez_compressed(os.path.join(OUT, 'data', 'groups.npz'), names=np.array(list(g)), masks=np.packbits(np.stack([g[k] for k in g]), axis=1), N=N)
print({k: int(v.sum()) for k, v in g.items()})
