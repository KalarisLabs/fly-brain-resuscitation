"""Recoverability-index track: connectome-only indices of the SURVIVING FAFB v783 (princeton5) network
for every dynamics run (arms A,B,C,D,E,H,base). Alive masks come from Modal Volume flysim-data:/results/<arm>/<id>.npz
(packed 'alive'), extracted by job 8f28b96c-e7de-4a84-937b-3c33524da666 into data/alive_masks.npz.
Ageing arm E: synapse(edge)-loss masks are regenerated exactly as dynamics/sweep.py make_aged_model (rng 10000+seed).

Indices (all on the post-insult surviving graph; D = diag(alive)):
  rho_exc       Perron root of D*We*D (We = excitatory synapse counts, sign as in the simulator: ACh/unknown/modulatory = +1)
  Reff_mf       rho_exc * w_syn*tau_syn/(tau_m*dV)   mean-field LIF slope (suprathreshold regime), = linear 'reproduction number'
  Reff_psp      rho_exc * peak-PSP/dV                single-spike peak PSP version
  rho_exc_known Perron root restricted to fast-transmitter-typed neurons (DA/SER/OCT/unknown, signed +1 by convention, excluded)
  ipr_exc       inverse participation ratio of the Perron vector (localisation)
  lam_signed    leading real eigenvalue of D*(We - Wi)*D (ARPACK, which='LR')
  EI_syn        surviving excitatory / inhibitory synapse count;  exc_in_frac = mean over alive neurons with input of E/(E+I)
  hub1/5/10     surviving fraction of top-1/5/10 % degree neurons (same sets as arm H)
  rich1_norm    directed edge density among alive top-1 % hubs / intact value
  recip_exc     surviving reciprocal excitatory pairs / intact;  loop2_super = surviving 2-loops with loop gain >1
  cyc3_exc      surviving directed excitatory 3-cycles / intact; loop3_super = surviving 3-loops with loop gain >1
  gscc_exc      giant SCC of the binary excitatory surviving graph (fraction of N)
  ecore_<r>     size (fraction of N) of the excitatory bootstrap core: maximal set where every member receives
                >= theta = dV/(w*tau_syn*r) excitatory synapses from other members (r = assumed saturated rate, Hz)
  score_<r>     signed version (synchronous iteration from all alive, inhibition subtracts; converged flag)
  reachN_L6, reachG_L6, influence_L6, redund_L6   perclib (Structure track) strict sensory->sink metrics
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
import os, sys, json, time, math
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.sparse.linalg as sla
from scipy.sparse.csgraph import connected_components

D = os.environ.get('RESEARCH_ROOT', _P.OUT)
OUT = os.path.join(D, 'followup', 'index')
sys.path.insert(0, os.path.join(D, 'structure')); sys.path.insert(0, os.path.join(D, 'followup', 'index'))
from light import load_light, LIF_DEFAULTS
import perclib as PL

W_SYN = 0.185
RSAT = (100, 200)


class FAFBIndex:
    def __init__(self, use_perclib=True):
        self.use_perclib = use_perclib
        M = load_light(); self.M = M
        N = self.N = M.N
        L = LIF_DEFAULTS; dV = L['v_th'] - L['v_rest']; ts, tm = L['tau_syn'], L['tau_m']
        self.GMF = W_SYN * ts / (tm * dV)
        tp = math.log(tm / ts) * tm * ts / (tm - ts)
        self.GPSP = W_SYN * ts / (tm - ts) * (math.exp(-tp / tm) - math.exp(-tp / ts)) / dV
        self.theta = {r: dV / (W_SYN * ts * r / 1000.0) for r in RSAT}
        rows = np.repeat(np.arange(N, dtype=np.int32), np.diff(M.indptr)); cols = M.indices.astype(np.int32)
        self.rows, self.cols, self.wt = rows, cols, M.weight.astype(np.float32)
        ne = self.ne = len(cols)
        exc = self.wt > 0; inh = ~exc
        # transposed CSRs (rows = post) whose data = edge id, so weights can be re-derived for synapse-loss masks
        self.ETid = sp.csr_matrix((np.flatnonzero(exc) + 1, (cols[exc], rows[exc])), shape=(N, N)).tocsr()
        self.ITid = sp.csr_matrix((np.flatnonzero(inh) + 1, (cols[inh], rows[inh])), shape=(N, N)).tocsr()
        self.ETid.data -= 1; self.ITid.data -= 1
        self.known = M.fast   # fast-transmitter-typed (ACh/GABA/Glu); DA/SER/OCT/unknown = sign-uncertain
        deg = np.diff(M.indptr) + np.bincount(M.indices, minlength=N); self.deg = deg
        rk = np.argsort(-deg, kind='stable')
        self.HUB = {p: np.zeros(N, bool) for p in (1, 5, 10)}
        for p in self.HUB: self.HUB[p][rk[:int(round(N * p / 100))]] = True
        self.base_demand = M.base_demand
        self._loops()
        if use_perclib:
            self.G = PL.Graph('FAFB'); self.G.setup_bits()
            assert (self.G.nodes.root_id.values == M.neurons.root_id.values).all()
        self.ref = None
        self.ref = self.compute(np.ones(N, bool))

    # ------------------------------------------------------------------ loops
    def _loops(self):
        N, rows, cols, wt = self.N, self.rows, self.cols, self.wt
        cp = os.path.join(OUT, 'data', 'loops_cache.npz')
        if os.path.exists(cp):
            z = np.load(cp); self.R2, self.C3, self.E3 = z['R2'], z['C3'], z['E3']
            self.g2 = np.prod(wt[self.R2] * self.GMF, 1); self.g3 = np.prod(wt[self.E3] * self.GMF, 1)
            return
        exc = wt > 0; eid_exc = np.flatnonzero(exc)
        key = rows[exc].astype(np.int64) * N + cols[exc]; srt = np.argsort(key); ks = key[srt]
        def eid(a, b):
            q = a.astype(np.int64) * N + b; p = np.searchsorted(ks, q); p[p >= ks.size] = 0
            ok = ks[p] == q; return np.where(ok, eid_exc[srt[p]], -1)
        rev = eid(cols[exc], rows[exc]); m = (rev >= 0) & (rows[exc] < cols[exc])
        self.R2 = np.c_[eid_exc[m], rev[m]]
        Ab = sp.csr_matrix((np.ones(exc.sum(), np.float32), (rows[exc], cols[exc])), shape=(N, N))
        AbT = Ab.T.tocsr(); cyc = []
        for s in range(0, N, 8000):
            P = (Ab[s:s + 8000] @ Ab).multiply(AbT[s:s + 8000]).tocoo()
            for i, k in zip(P.row + s, P.col):
                if not (i < k): continue
                js = np.intersect1d(Ab.indices[Ab.indptr[i]:Ab.indptr[i + 1]], AbT.indices[AbT.indptr[k]:AbT.indptr[k + 1]], assume_unique=True)
                for j in js:
                    if i < j: cyc.append((i, j, k))
        C3 = np.array(cyc, np.int32).reshape(-1, 3)
        self.C3 = C3; self.E3 = np.c_[eid(C3[:, 0], C3[:, 1]), eid(C3[:, 1], C3[:, 2]), eid(C3[:, 2], C3[:, 0])]
        assert (self.E3 >= 0).all()
        self.g2 = np.prod(wt[self.R2] * self.GMF, 1); self.g3 = np.prod(wt[self.E3] * self.GMF, 1)
        os.makedirs(os.path.dirname(cp), exist_ok=True)
        np.savez_compressed(cp, R2=self.R2, C3=self.C3, E3=self.E3)

    # ------------------------------------------------------------------ helpers
    def mats(self, keep=None):
        if keep is None and getattr(self, '_m0', None) is not None:
            return self._m0
        w = self.wt if keep is None else self.wt * keep
        WeT = sp.csr_matrix((w[self.ETid.data], self.ETid.indices, self.ETid.indptr), shape=(self.N, self.N))
        WiT = sp.csr_matrix((-w[self.ITid.data], self.ITid.indices, self.ITid.indptr), shape=(self.N, self.N))
        if keep is None: self._m0 = (WeT, WiT)
        return WeT, WiT

    def perron(self, AT, a, x0=None, it=1500, tol=1e-5):
        x = (np.ones(self.N, np.float32) if x0 is None else x0.astype(np.float32).copy()) * a + 1e-6 * a
        nx = np.linalg.norm(x)
        if nx == 0: return 0.0, x
        x /= nx; lo = 0.0
        for k in range(it):
            y = a * (AT @ x) + x
            lam = float(x @ y); ny = np.linalg.norm(y); x = y / ny
            if k > 20 and abs(lam - lo) < tol * max(1.0, abs(lam)): break
            lo = lam
        return lam - 1.0, x

    def ecore(self, WeT, alive, theta):
        S = alive.copy()
        for _ in range(500):
            inp = WeT @ S.astype(np.float32)
            S2 = S & (inp >= theta)
            if S2.sum() == S.sum(): break
            S = S2
        return S

    def score(self, WeT, WiT, alive, theta):
        S = alive.copy(); conv = False
        for _ in range(40):
            s = S.astype(np.float32); inp = WeT @ s - WiT @ s
            S2 = alive & (inp >= theta)
            if (S2 == S).all(): conv = True; break
            S = S2
        return S, conv

    # ------------------------------------------------------------------ main
    def compute(self, alive, keep=None, perclib=True):
        N = self.N; a = alive.astype(np.float32); o = {}
        WeT, WiT = self.mats(keep)
        ref = self.ref
        o['dead'] = float(1 - alive.mean())
        x0 = None if ref is None else self._v
        r, v = self.perron(WeT, a, x0)
        if ref is None: self._v = v
        o['rho_exc'] = r; o['Reff_mf'] = r * self.GMF; o['Reff_psp'] = r * self.GPSP
        o['ipr_exc'] = float((v ** 4).sum() / max((v ** 2).sum() ** 2, 1e-30))
        ak = a * self.known
        rk_, vk = self.perron(WeT, ak, None if ref is None else self._vk)
        if ref is None: self._vk = vk
        o['rho_exc_known'] = rk_; o['Reff_mf_known'] = rk_ * self.GMF
        # signed leading eigenvalue
        Wd = sla.LinearOperator((N, N), matvec=lambda x: (a * (WeT @ (a * x.astype(np.float32))) - a * (WiT @ (a * x.astype(np.float32)))).astype(np.float64), dtype=np.float64)
        try:
            ev, evec = sla.eigs(Wd, k=1, which='LR', tol=1e-5, maxiter=20000, v0=(None if ref is None else (self._vs * a + 1e-6).astype(np.float64)))
            o['lam_signed'] = float(ev.real.max())
            if ref is None: self._vs = np.abs(evec[:, 0].real)
        except Exception as e:
            o['lam_signed'] = np.nan
        o['Reff_signed'] = o['lam_signed'] * self.GMF
        # E/I among survivors
        sa = a
        Ein = WeT @ sa * sa; Iin = WiT @ sa * sa
        o['EI_syn'] = float(Ein.sum() / max(Iin.sum(), 1e-9))
        tot = Ein + Iin; m = tot > 0
        o['exc_in_frac'] = float((Ein[m] / tot[m]).mean()) if m.any() else np.nan
        for p, h in self.HUB.items(): o[f'hub{p}'] = float(alive[h].mean())
        h = self.HUB[1] & alive; nh = h.sum()
        e_in = h[self.rows] & h[self.cols]
        if keep is not None: e_in &= keep.astype(bool)
        o['rich1'] = float(e_in.sum() / max(nh * (nh - 1), 1))
        kp = np.ones(self.ne, bool) if keep is None else keep.astype(bool)
        l2 = alive[self.rows[self.R2[:, 0]]] & alive[self.cols[self.R2[:, 0]]] & kp[self.R2[:, 0]] & kp[self.R2[:, 1]]
        l3 = alive[self.C3].all(1) & kp[self.E3].all(1)
        o['recip_exc'] = float(l2.sum()); o['loop2_super'] = float((l2 & (self.g2 > 1)).sum())
        o['cyc3_exc'] = float(l3.sum()); o['loop3_super'] = float((l3 & (self.g3 > 1)).sum())
        Ab = WeT.copy(); Ab.data = (Ab.data > 0).astype(np.float32)
        Ab = sp.diags(a) @ Ab @ sp.diags(a); Ab.eliminate_zeros()
        nc, lab = connected_components(Ab, directed=True, connection='strong')
        o['gscc_exc'] = float(np.bincount(lab[alive]).max() / N) if alive.any() else 0.0
        for rr in RSAT:
            o[f'ecore{rr}'] = float(self.ecore(WeT, alive, self.theta[rr]).mean())
        o['score200'] = np.nan; o['score200_conv'] = False   # dropped after pass 1b (non-convergent, costly)
        if perclib and self.use_perclib:
            mm = self.G.evaluate(alive, do_inf=False, do_scc=False, maxL=6)
            for k in ['reachN_L6', 'reachN_L4', 'reachG_L6', 'influence_L6', 'redund_L6']: o[k] = float(mm[k])
        return o

    def aged_keep(self, age, mode, seed):
        p = 0.3 * age; rng = np.random.default_rng(10_000 + seed)
        if mode == 'random':
            pe = np.full(self.ne, p)
        else:
            b = np.repeat(self.base_demand, np.diff(self.M.indptr)); pe = np.clip(p * b / b.mean(), 0, 0.95)
        return (rng.random(self.ne) >= pe)


def main(shard=0, nshard=1, use_perclib=True):
    z = np.load(os.path.join(OUT, 'data', 'alive_masks.npz')); ids = z['ids']; packed = z['packed']
    runs = pd.read_csv(os.path.join(OUT, 'stage', 'runs_slim.csv')).set_index('id')
    t0 = time.time(); X = FAFBIndex(use_perclib); print('init', round(time.time() - t0, 1), 's', X.ref, flush=True)
    json.dump({k: (float(v) if not isinstance(v, (bool, np.bool_)) else bool(v)) for k, v in X.ref.items()},
              open(os.path.join(OUT, 'data', 'fafb_intact_index.json'), 'w'), indent=1)
    outp = os.path.join(OUT, 'data', f'run_indices_{shard}.csv' if nshard > 1 else 'run_indices.csv')
    import glob as _g
    done = set()
    for fp in _g.glob(os.path.join(OUT, 'data', 'run_indices*.csv')):
        done |= set(pd.read_csv(fp).id)
    # unique (mask, keep-spec) jobs
    jobs = {}
    for i, rid in enumerate(ids):
        r = runs.loc[rid]
        kspec = None
        if r.arm == 'E' and r.get('age', 0) and float(r.age) > 0:
            kspec = (float(r.age), str(r.age_mode), int(r.seed))
        jobs.setdefault((packed[i].tobytes(), kspec), []).append(rid)
    print('unique jobs', len(jobs), flush=True)
    cache = {}
    buf = []
    prio = {'D': 0, 'E': 1, 'H': 2, 'B': 3, 'C': 4}
    items = sorted(jobs.items(), key=lambda kv: prio.get(runs.loc[kv[1][0]].arm, 5))
    for n, ((pb, kspec), rids) in enumerate(items):
        if n % nshard != shard: continue
        if all(r in done for r in rids): continue
        alive = np.unpackbits(np.frombuffer(pb, np.uint8))[:X.N].astype(bool)
        keep = X.aged_keep(*kspec) if kspec else None
        t = time.time()
        for attempt in range(8):
            try:
                o = X.compute(alive, keep); break
            except MemoryError:
                import gc; gc.collect(); print('MemoryError, retry', attempt, flush=True); time.sleep(30)
        else:
            raise MemoryError('persistent')
        o['syn_lost_edges'] = float(1 - keep.mean()) if keep is not None else np.nan
        o['t_s'] = round(time.time() - t, 2)
        for rid in rids: buf.append(dict(id=rid, **o))
        if len(buf) >= 10:
            df = pd.DataFrame(buf)
            if os.path.exists(outp): df = df.reindex(columns=list(pd.read_csv(outp, nrows=0).columns))
            df.to_csv(outp, mode='a', header=not os.path.exists(outp), index=False); buf = []
            print(n, len(jobs), round(time.time() - t0), flush=True)
    if buf:
        df = pd.DataFrame(buf)
        if os.path.exists(outp): df = df.reindex(columns=list(pd.read_csv(outp, nrows=0).columns))
        df.to_csv(outp, mode='a', header=not os.path.exists(outp), index=False)
    print('done', round(time.time() - t0), flush=True)


if __name__ == '__main__':
    a = sys.argv[1:]
    if a and a[0] == 'precompute':
        X = FAFBIndex(use_perclib=False); print('precomputed loops', X.C3.shape, X.R2.shape, flush=True); sys.exit(0)
    main(int(a[0]), int(a[1])) if len(a) == 2 else main()
