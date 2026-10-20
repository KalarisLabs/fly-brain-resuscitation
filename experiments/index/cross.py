"""Cross-connectome recoverability map: candidate indices on intact and damaged FAFB / BANC / MCNS / MANC graphs.
Usage: python cross.py <DS> <regime> [n_orders]
Graph: Structure-track caches (structure/cache/<DS>_edges.npz, pairs >= 5 synapses; perclib.Graph roles / removal orders).
Sign: presynaptic predicted NT (GABA, GLUT -> -1; ACH, modulatory, unknown -> +1), the simulator's princeton5 convention.
Damage regimes and the f grid follow structure/perclib.removal_order / neuropil_order (removable interneurons only).
Writes followup/index/data/cross_<DS>_<regime>.csv
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
import os, sys, math, time
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.sparse.linalg as sla
from scipy.sparse.csgraph import connected_components
D = os.environ.get('RESEARCH_ROOT', _P.OUT); RAW = os.environ.get('RAW_ROOT', _P.RAW)
OUT = os.path.join(D, 'followup', 'index')
sys.path.insert(0, os.path.join(D, 'structure')); import perclib as PL

W_SYN, TS, TM, DV = 0.185, 5.0, 20.0, 7.0
GMF = W_SYN * TS / (TM * DV)
FGRID = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
NEUR = {'FAFB': 'neurons (2).csv.gz', 'BANC': 'neurons (3).csv.gz', 'MANC': 'neurons (4).csv.gz', 'MCNS': 'neurons.csv.gz'}


def nt_sign(ds, root_ids):
    if ds == 'FAFB':
        t = pd.read_csv(os.path.join(RAW, NEUR[ds]), usecols=['root_id', 'nt_type']).rename(columns={'nt_type': 'nt'})
    else:
        t = pd.read_csv(os.path.join(RAW, NEUR[ds]), usecols=['Root ID', 'Predicted NT type'], low_memory=False)
        t = t.rename(columns={'Root ID': 'root_id', 'Predicted NT type': 'nt'})
    t = t.drop_duplicates('root_id').set_index('root_id').nt
    nt = t.reindex(root_ids).astype(str).str.upper().values
    known = np.array([any(k in x for k in ('ACH', 'ACETYL', 'GABA', 'GLU')) for x in nt])
    inh = np.array([('GABA' in x) or ('GLU' in x) for x in nt])
    return np.where(inh, -1, 1).astype(np.int8), known


class XIndex:
    def __init__(self, ds):
        self.ds = ds
        g = PL.Graph(ds); g.setup_bits(); self.g = g; n = self.n = g.n
        e = np.load(os.path.join(D, 'structure', 'cache', f'{ds}_edges.npz'))
        m = e['syn'] >= 5; pre, post, syn = e['pre'][m], e['post'][m], e['syn'][m].astype(np.float32)
        k = pre != post; pre, post, syn = pre[k], post[k], syn[k]
        sign, known = nt_sign(ds, g.nodes.root_id.values); self.known = known
        s = sign[pre]
        self.WeT = sp.csr_matrix((syn[s > 0], (post[s > 0], pre[s > 0])), shape=(n, n))
        self.WiT = sp.csr_matrix((syn[s < 0], (post[s < 0], pre[s < 0])), shape=(n, n))
        deg = np.bincount(pre, minlength=n) + np.bincount(post, minlength=n); rk = np.argsort(-deg, kind='stable')
        self.HUB = {p: np.zeros(n, bool) for p in (1, 5, 10)}
        for p in self.HUB: self.HUB[p][rk[:int(round(n * p / 100))]] = True
        self.info = dict(n=n, n_edges=int(pre.size), syn_total=float(syn.sum()), mean_syn_per_edge=float(syn.mean()),
                         frac_inh_neurons=float((sign < 0).mean()), frac_unknown_nt=float(1 - known.mean()),
                         exc_syn=float(syn[s > 0].sum()), inh_syn=float(syn[s < 0].sum()))
        self._v = None; self._vk = None

    def perron(self, a, x0, it=1500, tol=1e-5):
        x = (np.ones(self.n, np.float32) if x0 is None else x0.copy()) * a + 1e-6 * a
        nx = np.linalg.norm(x)
        if nx == 0: return 0.0, x
        x /= nx; lo = 0
        for k in range(it):
            y = a * (self.WeT @ x) + x; lam = float(x @ y); x = y / np.linalg.norm(y)
            if k > 20 and abs(lam - lo) < tol * max(1, abs(lam)): break
            lo = lam
        return lam - 1, x

    def compute(self, alive):
        a = alive.astype(np.float32); o = {'alive_frac': float(alive.mean())}
        r, v = self.perron(a, self._v)
        if self._v is None: self._v = v
        o['rho_exc'] = r; o['Reff_mf'] = r * GMF
        o['ipr_exc'] = float((v ** 4).sum() / max((v ** 2).sum() ** 2, 1e-30))
        rk_, vk = self.perron(a * self.known, self._vk)
        if self._vk is None: self._vk = vk
        o['Reff_mf_known'] = rk_ * GMF
        Wd = sp.diags(a.astype(np.float64)) @ (self.WeT - self.WiT).astype(np.float64) @ sp.diags(a.astype(np.float64))
        try:
            ev = sla.eigs(Wd, k=1, which='LR', tol=1e-5, maxiter=20000, return_eigenvectors=False)
            o['lam_signed'] = float(ev.real.max())
        except Exception:
            o['lam_signed'] = np.nan
        o['Reff_signed'] = o['lam_signed'] * GMF
        Ein = self.WeT @ a * a; Iin = self.WiT @ a * a
        o['EI_syn'] = float(Ein.sum() / max(Iin.sum(), 1e-9))
        for p, h in self.HUB.items(): o[f'hub{p}'] = float(alive[h].mean())
        for rr in (100, 200):
            th = DV / (W_SYN * TS * rr / 1000.0); S = alive.copy()
            for _ in range(500):
                S2 = S & (self.WeT @ S.astype(np.float32) >= th)
                if S2.sum() == S.sum(): break
                S = S2
            o[f'ecore{rr}'] = float(S.mean())
        Ab = self.WeT.copy(); Ab.data[:] = 1; Ab = sp.diags(a) @ Ab @ sp.diags(a); Ab.eliminate_zeros()
        nc, lab = connected_components(Ab, directed=True, connection='strong')
        o['gscc_exc'] = float(np.bincount(lab[alive]).max() / self.n) if alive.any() else 0.0
        # Reach / influence are taken from the Structure-track curves (structure/tables/curves_normalised_all.csv.gz)
        pass
        return o


def main(ds, regime, n_orders=2):
    t0 = time.time(); X = XIndex(ds); g = X.g
    print(ds, 'init', round(time.time() - t0), X.info, flush=True)
    rows = []
    n_orders = 2 if regime in ('focal', 'neuropil') else 1   # local-compute budget: 1 order for random/metabolic/degree
    for oi in range(n_orders):
        rng = np.random.default_rng(1000 * oi + 7)
        if regime in ('neuropil', 'focal'):
            order, core = PL.neuropil_order(ds, g, rng)
        else:
            dem = g.nodes.syn_total.values if regime == 'metabolic_syn' else None
            order = PL.removal_order(g, regime, rng, demand=dem); core = ''
        N = g.removable.size
        for f in FGRID:
            pass
            alive = np.ones(g.n, bool); alive[order[:int(round(f * N))]] = False
            t = time.time(); o = X.compute(alive)
            o.update(ds=ds, regime=regime, order=oi, core=core, f=f, t_s=round(time.time() - t, 1), **{'info_' + k: v for k, v in X.info.items()})
            rows.append(o); print(ds, regime, oi, f, round(o['Reff_mf'], 3), round(o['ecore100'], 4), o['t_s'], flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(OUT, 'data', f'cross_{ds}_{regime}.csv'), index=False)
    print('done', round(time.time() - t0), flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 3)
