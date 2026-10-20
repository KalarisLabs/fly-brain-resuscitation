"""Percolation / reachability engine (pure numpy + scipy; no networkx on whole graphs).

Reach is computed with bitset propagation: each sensory modality group is one bit of a
uint64; one hop = gather over in-edges + np.bitwise_or.reduceat. Reach_L = fraction of
(modality group, sink) pairs with a directed path of <= L hops through surviving neurons.
Path redundancy = sum over sinks of log1p(# walks of length 1..6 from the source set) on the
binarized graph; weighted influence = summed product of input-normalised weights over walks
<= 6 hops (prereg 'weighted variant'). Giant SCC via scipy.sparse.csgraph.
"""
import os, gc, json
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

ROOT = os.path.join(os.environ.get("RESEARCH_ROOT", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs")), 'structure')
LS = (2, 4, 6)


_T8 = np.array([bin(i).count('1') for i in range(256)], np.int64)


def popcount64(x):
    """Population count of uint64 words (np.bitwise_count needs numpy>=2.0; byte-table fallback otherwise)."""
    if hasattr(np, "bitwise_count"):
        return np.bitwise_count(x).astype(np.int64)
    x = np.ascontiguousarray(x)
    return _T8[x.view(np.uint8)].reshape(x.shape + (8,)).sum(-1)


class Graph:
    def __init__(self, ds, thr=5, rewire_seed=None):
        self.ds = ds
        self.nodes = pd.read_csv(os.path.join(ROOT, 'cache', f'{ds}_nodes.csv.gz'), low_memory=False)
        e = np.load(os.path.join(ROOT, 'cache', f'{ds}_edges.npz'))
        n = int(e['n']); self.n = n
        m = e['syn'] >= thr
        pre, post, syn = e['pre'][m], e['post'][m], e['syn'][m].astype(np.float64)
        del e
        if rewire_seed is not None:  # directed configuration null: shuffle targets, keep in/out stubs
            rng = np.random.default_rng(rewire_seed)
            post = post[rng.permutation(post.size)]
            k = pre != post
            pre, post, syn = pre[k], post[k], syn[k]
            key = pre.astype(np.int64) * n + post
            _, u = np.unique(key, return_index=True)
            pre, post, syn = pre[u], post[u], syn[u]
        self.n_edges = int(pre.size)
        role = self.nodes.role.values
        excl = role == 'excluded'
        k = ~(excl[pre] | excl[post])
        pre, post, syn = pre[k], post[k], syn[k]
        self.pre, self.post = pre.astype(np.int32), post.astype(np.int32)
        # AT: rows=post, cols=pre (in-edges)
        AT = sp.csr_matrix((np.ones(pre.size), (post, pre)), shape=(n, n))
        AT.sum_duplicates(); AT.sort_indices()
        self.AT = AT
        insum = np.bincount(post, weights=syn, minlength=n)
        ATw = sp.csr_matrix((syn / np.maximum(insum[post], 1), (post, pre)), shape=(n, n))
        ATw.sum_duplicates(); ATw.sort_indices()
        self.ATw = ATw
        self.src = np.where(role == 'source')[0]
        self.sink = np.where(role == 'sink')[0]
        self.removable = np.where(role == 'removable')[0]
        self.excluded = np.where(excl)[0]
        grp = self.nodes.modality_group.values
        g_src = grp[self.src]
        self.groups = sorted(pd.unique(g_src))
        gi = {g: i for i, g in enumerate(self.groups)}
        self.G = len(self.groups)
        bits = np.zeros(n, np.uint64)
        bits[self.src] = np.left_shift(np.uint64(1), np.array([gi[g] for g in g_src], dtype=np.uint64))
        self.bits0 = bits
        sk = self.nodes.sink_kind.values[self.sink]
        self.sink_sub = {'all': np.ones(self.sink.size, bool), 'DN': sk == 'DN', 'MN': sk == 'MN'}
        self.outdeg = np.bincount(pre, minlength=n); self.indeg = np.bincount(post, minlength=n)
        self.deg = self.outdeg + self.indeg

    # ---------------------------------------------------------------
    def _masked(self, M, alive, rows=None):
        """Return CSR with all edges touching dead nodes removed (rows = precomputed row index of each nnz)."""
        if alive.all():
            return M
        if rows is None:
            rows = np.repeat(np.arange(self.n, dtype=np.int32), np.diff(M.indptr))
        keep = alive[M.indices] & alive[rows]
        cnt = np.bincount(rows[keep], minlength=self.n)
        indptr = np.r_[0, np.cumsum(cnt)].astype(np.int64)
        out = sp.csr_matrix((M.data[keep], M.indices[keep], indptr), shape=M.shape)
        out.has_sorted_indices = True
        return out

    def setup_bits(self, n_sample=128, seed=12345):
        """Words: 0 = modality-group bits (strict relay), 1 = modality-group bits (loose relay),
        2.. = individual source-neuron sample bits (strict relay), stratified by modality group."""
        rng = np.random.default_rng(seed)
        grp = self.nodes.modality_group.values[self.src]
        # proportional allocation with >=1 per group
        vc = pd.Series(grp).value_counts()
        alloc = np.maximum(1, np.floor(vc / vc.sum() * n_sample)).astype(int)
        while alloc.sum() > n_sample: alloc[alloc.idxmax()] -= 1
        while alloc.sum() < n_sample: alloc[(vc / alloc).idxmax()] += 1
        pick = np.concatenate([rng.choice(self.src[grp == k], size=min(a, (grp == k).sum()), replace=False) for k, a in alloc.items()])
        self.sample_src = pick
        nw = 2 + int(np.ceil(pick.size / 64))
        B = np.zeros((self.n, nw), np.uint64)
        B[:, 0] = self.bits0; B[:, 1] = self.bits0
        for i, s in enumerate(pick):
            B[s, 2 + i // 64] |= np.left_shift(np.uint64(1), np.uint64(i % 64))
        self.B0 = B; self.nw = nw; self.n_sample = pick.size
        self.strict_cols = np.array([0] + list(range(2, nw)))
        self.is_src = np.zeros(self.n, bool); self.is_src[self.src] = True
        self.is_sink = np.zeros(self.n, bool); self.is_sink[self.sink] = True

    def evaluate(self, alive, do_inf=True, do_scc=True, do_walk=True, maxL=60):
        """alive: bool array (n). Returns dict of raw metrics.
        Strict relay: paths may pass only through non-source, non-sink neurons (sources start paths, sinks end them).
        Frontier-delta propagation: OR is monotone, so each hop only gathers edges out of neurons whose relay state changed."""
        if not hasattr(self, 'B0'): self.setup_bits()
        if not hasattr(self, 'AT_rows'):
            self.AT_rows = np.repeat(np.arange(self.n, dtype=np.int32), np.diff(self.AT.indptr))
            self.B0T = np.ascontiguousarray(self.B0.T)
        AT = self._masked(self.AT, alive, self.AT_rows)
        rows = np.repeat(np.arange(self.n, dtype=np.int32), np.diff(AT.indptr))
        indices = AT.indices
        out = {}
        G = self.G; S = self.sink.size; sc = self.strict_cols; nw = self.nw
        R = self.B0T.copy(); R[:, ~alive] = 0
        srcB = self.B0T[np.ix_(sc, self.src)]
        def mk_relay(R):
            rl = R.copy(); rl[np.ix_(sc, self.src)] = srcB; rl[np.ix_(sc, self.sink)] = 0
            return rl
        def record(tag, R):
            rs = R[:, self.sink]
            pg = popcount64(rs[0]); pl = popcount64(rs[1]); pn = popcount64(rs[2:]).sum(0)
            for k, msk in self.sink_sub.items():
                if not msk.sum(): continue
                suf = '' if k == 'all' else f'_{k}'
                out[f'reachN_{tag}{suf}'] = pn[msk].sum() / (self.n_sample * msk.sum())   # primary: neuron pairs, strict
                out[f'reachG_{tag}{suf}'] = pg[msk].sum() / (G * msk.sum())              # modality-group pairs, strict
                out[f'reachGloose_{tag}{suf}'] = pl[msk].sum() / (G * msk.sum())         # modality-group pairs, loose relay
        prev = np.zeros_like(R)
        relay = mk_relay(R)
        L = 0
        while True:
            ch = (relay != prev).any(0)
            if not ch.any() or L >= maxL:
                break
            L += 1
            sel = np.flatnonzero(ch[indices])
            if sel.size:
                rs_ = rows[sel]; cs_ = indices[sel]
                bnd = np.flatnonzero(np.r_[True, rs_[1:] != rs_[:-1]])
                ur = rs_[bnd]
                for w in range(nw):
                    R[w, ur] |= np.bitwise_or.reduceat(relay[w, cs_], bnd)
            prev = relay; relay = mk_relay(R)
            if L in LS: record(f'L{L}', R)
        for l in LS:
            if l > L: record(f'L{l}', R)
        if do_inf:
            record('Linf', R); out['L_converge'] = L
        if do_walk:
            ATw = self._masked(self.ATw, alive, self.AT_rows)
            x = np.zeros(self.n); x[self.src] = 1.0; x[~alive] = 0
            xw = x.copy(); W = np.zeros(S); Iw = np.zeros(S)
            relaymask = self.is_src | self.is_sink
            for k in range(6):
                x = AT @ x; xw = ATw @ xw
                W += x[self.sink]; Iw += xw[self.sink]
                x[relaymask] = 0; xw[relaymask] = 0   # strict: interneuron relays only
            out['redund_L6'] = float(np.log1p(W).sum())
            out['influence_L6'] = float(Iw.sum())
        if do_scc:
            nc, lab = connected_components(AT, directed=True, connection='strong')
            lab = lab[alive]
            out['gscc'] = float(np.bincount(lab).max() / self.n) if lab.size else 0.0
        return out


# ---------------------------------------------------------------
def neuropil_order(ds, g, rng, core=None):
    """Growing focal lesion: grow from a core neuropil along the neuropil co-innervation graph."""
    z = np.load(os.path.join(ROOT, 'cache', f'{ds}_neuron_neuropil.npz'), allow_pickle=True)
    names = z['names']; P = names.size
    B = sp.csr_matrix((z['s'].astype(float), (z['nid'], z['np'])), shape=(g.n, P))
    Bn = sp.diags(1 / np.maximum(np.sqrt(np.asarray(B.multiply(B).sum(1)).ravel()), 1e-9)) @ B
    C = (Bn.T @ Bn).toarray(); np.fill_diagonal(C, 0)
    dom = g.nodes.dominant_neuropil.values
    rem = g.removable
    dom_r = pd.Series(dom[rem]).fillna('NONE').values
    counts = pd.Series(dom_r).value_counts()
    elig = [p for p in counts.index if p != 'NONE' and counts[p] >= 0.005 * rem.size]
    if core is None:
        core = rng.choice(elig)
    pidx = {p: i for i, p in enumerate(names)}
    order_np = [core]; S = np.zeros(P, bool); S[pidx[core]] = True
    score = C[pidx[core]].copy()
    while S.sum() < P:
        sc = np.where(S, -1, score)
        j = int(np.argmax(sc)); S[j] = True; order_np.append(names[j]); score += C[j]
    rank = {p: i for i, p in enumerate(order_np)}
    rk = np.array([rank.get(p, P + 1) for p in dom_r], float)
    rk = rk + rng.random(rk.size) * 0.5  # random within neuropil
    return rem[np.argsort(rk)], core


def removal_order(g, regime, rng, demand=None):
    rem = g.removable
    if regime == 'random':
        return rng.permutation(rem)
    if regime.startswith('metabolic'):
        w = np.asarray(demand[rem], float)
        w = np.where(np.isfinite(w) & (w > 0), w, np.nanmedian(w[w > 0]))
        keys = rng.random(rem.size) ** (1.0 / w)   # Efraimidis-Spirakis weighted sampling w/o replacement
        return rem[np.argsort(-keys)]
    if regime == 'degree':
        d = g.deg[rem] + rng.random(rem.size) * 0.1
        return rem[np.argsort(-d)]
    if regime == 'degree_matched_random':
        # random neurons drawn from the same log-degree bins as the degree-targeted sequence
        d = g.deg[rem]
        b = np.digitize(np.log1p(d), np.linspace(0, np.log1p(d.max()) + 1e-9, 41))
        tgt = rem[np.argsort(-(d + rng.random(rem.size) * 0.1))]
        pools = {k: list(rng.permutation(rem[b == k])) for k in np.unique(b)}
        bt = b[np.searchsorted(rem, tgt)] if np.all(np.diff(rem) > 0) else None
        out = []
        for k in bt:
            out.append(pools[k].pop())
        return np.array(out)
    raise ValueError(regime)


def run_curve(g, order, fgrid, **kw):
    kw.setdefault('maxL', 30)
    rows = []
    N = g.removable.size
    for f in fgrid:
        alive = np.ones(g.n, bool); alive[g.excluded] = False
        alive[order[:int(round(f * N))]] = False
        m = g.evaluate(alive, **kw); m['f'] = f
        m.setdefault('gscc', np.nan)
        rows.append(m)
    return rows
