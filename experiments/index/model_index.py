"""Recoverability-index track, step 2: predict dynamics outcomes from connectome-only indices of the surviving graph.
numpy/scipy only (scikit-learn unavailable: no disk space to install). Inputs: data/run_indices.csv (+ data/run_perclib.csv),
dynamics/tables/all_runs.pkl. Outputs: tables/*.csv, tables/model_results.json.
CV: grouped by condition name (all seeds of a condition in the same fold), 5 folds x 20 repeats; cluster bootstrap over
conditions (1000 resamples) for CIs; leave-arm-out transfer (train on B+C, test on D / E / H)."""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src",):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import os, json, numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr, rankdata

D = os.environ.get('RESEARCH_ROOT', _P.OUT); OUT = os.path.join(D, 'followup', 'index')
TAB = os.path.join(OUT, 'tables'); os.makedirs(TAB, exist_ok=True)
RNG = np.random.default_rng(20261018)

RESET = {'reset', 'revive_reset', 'clr3_reset', 'noise_reset', 'adapt_reset', 'sdblock_reset', 'block10', 'block20', 'block40',
         'block20_del5', 'block20_del20', 'staged_fwd', 'staged_rev'}
NORESET = {'naive', 'naive_noinj', 'adapt', 'voff-10', 'voff-20', 'sdblock', 'clr3', 'ramp10', 'ramp30', 'ramp30_noinj', 'ramp60',
           'gaba2', 'gaba4', 'noise', 'block2', 'block5', 'pulsed3x10', 'pulsed6x5', 'cool12', 'cool17'}
# expected direction: +1 = larger value -> more recoverable (fixed a priori, before looking at outcomes)
FEATS = {'alive': +1, 'Reff_mf': +1, 'Reff_mf_known': +1, 'Reff_signed': +1, 'ipr_exc': -1, 'EI_syn': +1, 'exc_in_frac': +1,
         'hub1': +1, 'hub5': +1, 'hub10': +1, 'rich1_n': +1, 'recip_n': +1, 'loop2_super': +1, 'cyc3_n': +1, 'loop3_super': +1,
         'gscc_exc': +1, 'ecore100': +1, 'ecore200': +1, 'score200': +1, 'reach_sim': +1, 'reachN_L6': +1, 'reachG_L6': +1,
         'influence_n': +1, 'redund_L6': +1}


def load():
    import glob
    X = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(OUT, 'data', 'run_indices*.csv')))], ignore_index=True)
    X = X.drop_duplicates('id', keep='last')
    pc = os.path.join(OUT, 'data', 'run_perclib.csv')
    if os.path.exists(pc):
        X = X.merge(pd.read_csv(pc), on='id', how='left')
    ref = json.load(open(os.path.join(OUT, 'data', 'fafb_intact_index.json')))
    R = pd.read_pickle(os.path.join(D, 'dynamics', 'tables', 'all_runs.pkl'))
    R.index.name = None; R = R.reset_index(drop=True)
    R = R[['id', 'arm', 'name', 'seed', 'insult', 'D', 'reperf', 'T_ins', 'dead_frac', 'reach', 'Reach_rel', 'FII_norm', 'outcome',
           'protect', 'age', 'age_mode', 'syn_lost_frac', 'obs_late_Hz', 'FSA']]
    df = R.merge(X, on='id', how='inner')
    df['alive'] = 1 - df['dead']
    df['rich1_n'] = df.rich1 / ref['rich1']; df['recip_n'] = df.recip_exc / ref['recip_exc']; df['cyc3_n'] = df.cyc3_exc / ref['cyc3_exc']
    df['reach_sim'] = df['reach']
    if 'influence_L6' in df: df['influence_n'] = df.influence_L6 / df.loc[df.dead == 0, 'influence_L6'].median()
    df['reperf'] = df.reperf.fillna('')
    df['proto'] = np.where(df.reperf.isin(RESET), 'reset', np.where(df.reperf.isin(NORESET), 'noreset', 'other'))
    df['runaway'] = (df.outcome == 'runaway').astype(int)
    df['recovered'] = (df.outcome == 'recovered').astype(int)
    df['cls'] = df.outcome.map({'recovered': 'R', 'runaway': 'A', 'structural_loss': 'S', 'brain_death': 'X', 'impaired': 'I'})
    df['cond'] = df.arm + ':' + df.name
    return df, ref


def auc(y, s):
    y = np.asarray(y).astype(bool); s = np.asarray(s, float)
    n1, n0 = y.sum(), (~y).sum()
    if n1 == 0 or n0 == 0: return np.nan
    r = rankdata(s); return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def cluster_boot(df, fn, B=1000):
    conds = df.cond.unique(); idx = {c: np.flatnonzero(df.cond.values == c) for c in conds}
    out = []
    for b in range(B):
        pick = RNG.choice(conds, conds.size, replace=True)
        ii = np.concatenate([idx[c] for c in pick])
        v = fn(df.iloc[ii])
        if np.isfinite(v): out.append(v)
    return (np.percentile(out, 2.5), np.percentile(out, 97.5)) if out else (np.nan, np.nan)


def folds(df, k=5, rep=20):
    conds = df.cond.unique()
    for r in range(rep):
        perm = RNG.permutation(conds); grp = {c: i % k for i, c in enumerate(perm)}
        g = df.cond.map(grp).values
        for f in range(k):
            yield g != f, g == f


# ---------------- models
def std_fit(Xtr):
    mu = np.nanmean(Xtr, 0); sd = np.nanstd(Xtr, 0); sd[sd == 0] = 1; return mu, sd


def logit_fit(X, y, lam=1.0, w=None):
    n, p = X.shape; Xb = np.c_[np.ones(n), X]; w = np.ones(n) if w is None else w
    def f(b):
        z = Xb @ b; l = np.logaddexp(0, z) - y * z
        g = Xb.T @ (w * (1 / (1 + np.exp(-z)) - y)) + lam * np.r_[0, b[1:]]
        return (w * l).sum() + 0.5 * lam * (b[1:] ** 2).sum(), g
    return minimize(f, np.zeros(p + 1), jac=True, method='L-BFGS-B').x


def logit_pred(b, X): return 1 / (1 + np.exp(-(np.c_[np.ones(len(X)), X] @ b)))


def softmax_fit(X, Y, lam=1.0):
    n, p = X.shape; K = Y.shape[1]; Xb = np.c_[np.ones(n), X]
    def f(bf):
        B = bf.reshape(p + 1, K); Z = Xb @ B; Z -= Z.max(1, keepdims=True); P = np.exp(Z); P /= P.sum(1, keepdims=True)
        l = -(Y * np.log(P + 1e-12)).sum() + 0.5 * lam * (B[1:] ** 2).sum()
        G = Xb.T @ (P - Y); G[1:] += lam * B[1:]
        return l, G.ravel()
    return minimize(f, np.zeros((p + 1) * K), jac=True, method='L-BFGS-B').x.reshape(p + 1, K)


def softmax_pred(B, X):
    Z = np.c_[np.ones(len(X)), X] @ B; Z -= Z.max(1, keepdims=True); P = np.exp(Z); return P / P.sum(1, keepdims=True)


class Tree:
    """CART (gini for classification, MSE for regression) with random feature subsets."""
    def __init__(self, depth=6, min_leaf=5, mtry=None, reg=False, rng=None):
        self.depth, self.min_leaf, self.mtry, self.reg, self.rng = depth, min_leaf, mtry, reg, rng

    def _leaf(self, y): return y.mean(0)

    def _imp(self, y):
        if self.reg: return ((y - y.mean()) ** 2).sum()
        p = y.mean(0); p = np.array([p, 1 - p]) if np.ndim(p) == 0 else p; return len(y) * (1 - (p ** 2).sum())

    def fit(self, X, y):
        self.nodes = []; self._grow(X, y, 0); return self

    def _grow(self, X, y, d):
        i = len(self.nodes); self.nodes.append(None)
        if d >= self.depth or len(y) < 2 * self.min_leaf or self._imp(y) <= 1e-12:
            self.nodes[i] = ('L', self._leaf(y)); return i
        p = X.shape[1]; feats = self.rng.choice(p, self.mtry or p, replace=False)
        best = (self._imp(y), None, None)
        for j in feats:
            xs = X[:, j]; qs = np.unique(np.quantile(xs, np.linspace(0.05, 0.95, 19)))
            for t in qs:
                m = xs <= t; nl = m.sum()
                if nl < self.min_leaf or len(y) - nl < self.min_leaf: continue
                v = self._imp(y[m]) + self._imp(y[~m])
                if v < best[0] - 1e-12: best = (v, j, t)
        if best[1] is None:
            self.nodes[i] = ('L', self._leaf(y)); return i
        j, t = best[1], best[2]; m = X[:, j] <= t
        l = self._grow(X[m], y[m], d + 1); r = self._grow(X[~m], y[~m], d + 1)
        self.nodes[i] = ('S', j, t, l, r); return i

    def predict(self, X):
        out = []
        for x in X:
            n = self.nodes[0]
            while n[0] == 'S': n = self.nodes[n[3]] if x[n[1]] <= n[2] else self.nodes[n[4]]
            out.append(n[1])
        return np.array(out)


def rf_fit(X, y, reg=False, n=150, seed=0):
    rng = np.random.default_rng(seed); trees = []; mtry = max(1, int(np.sqrt(X.shape[1])))
    for b in range(n):
        ii = rng.integers(0, len(y), len(y))
        trees.append(Tree(6, 5, mtry, reg, rng).fit(X[ii], y[ii]))
    return trees


def rf_pred(trees, X): return np.mean([t.predict(X) for t in trees], 0)


def r2(y, p): return 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()


# ---------------- analyses
def univariate(df, target, feats, kind):
    rows = []
    sub = df.dropna(subset=[target])
    for f, sgn in feats.items():
        if f not in sub or sub[f].isna().all(): continue
        s = sub.dropna(subset=[f])
        if kind == 'auc':
            a = auc(s[target], sgn * s[f]); ci = cluster_boot(s, lambda q: auc(q[target], sgn * q[f]), 500)
            # CV of a one-feature logistic (orientation + calibration learnt in-fold)
            ps = np.full(len(s), np.nan); accs = []
            aucs = []
            for tr, te in folds(s, 5, 10):
                x = s[f].values[:, None]; mu, sd = std_fit(x[tr])
                b = logit_fit((x[tr] - mu) / sd, s[target].values[tr], 1e-3)
                aucs.append(auc(s[target].values[te], logit_pred(b, (x[te] - mu) / sd)))
            rows.append(dict(feature=f, n=len(s), n_pos=int(s[target].sum()), AUC_apriori=a, AUC_lo=ci[0], AUC_hi=ci[1],
                             AUC_cv_mean=np.nanmean(aucs), spearman=spearmanr(s[f], s[target]).correlation))
        else:
            y = s[target].values; x = s[f].values
            def cvr2(q):
                yy = q[target].values; xx = q[f].values; pr = np.zeros(len(q)); cnt = np.zeros(len(q))
                for tr, te in folds(q, 5, 4):
                    if tr.sum() < 5: continue
                    A = np.c_[np.ones(tr.sum()), xx[tr], xx[tr] ** 2]; c = np.linalg.lstsq(A, yy[tr], rcond=None)[0]
                    pr[te] += np.c_[np.ones(te.sum()), xx[te], xx[te] ** 2] @ c; cnt[te] += 1
                return r2(yy, pr / np.maximum(cnt, 1))
            v = cvr2(s); ci = cluster_boot(s, cvr2, 200)
            rows.append(dict(feature=f, n=len(s), R2_cv=v, R2_lo=ci[0], R2_hi=ci[1], spearman=spearmanr(x, y).correlation))
    return pd.DataFrame(rows)


def combined(df, target, feats, kind, extra=()):
    s = df.dropna(subset=[target] + list(feats)).reset_index(drop=True)
    X = s[list(feats) + list(extra)].values.astype(float); y = s[target].values.astype(float)
    res = {}
    for name in ['logit', 'rf']:
        pr = np.zeros(len(s)); cnt = np.zeros(len(s))
        for k, (tr, te) in enumerate(folds(s, 5, 4 if name == 'rf' else 20)):
            mu, sd = std_fit(X[tr]); Xtr = (X[tr] - mu) / sd; Xte = (X[te] - mu) / sd
            if name == 'logit':
                if kind == 'auc': p = logit_pred(logit_fit(Xtr, y[tr], 1.0), Xte)
                else:
                    A = np.c_[np.ones(tr.sum()), Xtr]; lam = 1.0 * np.eye(A.shape[1]); lam[0, 0] = 0
                    c = np.linalg.solve(A.T @ A + lam, A.T @ y[tr]); p = np.c_[np.ones(te.sum()), Xte] @ c
            else:
                p = rf_pred(rf_fit(Xtr, y[tr], reg=(kind != 'auc'), n=60, seed=k), Xte)
            pr[te] += p; cnt[te] += 1
        pr /= np.maximum(cnt, 1)
        if kind == 'auc':
            v = auc(y, pr); s['_p'] = pr; ci = cluster_boot(s, lambda q: auc(q[target], q['_p']), 500)
            res[name] = dict(AUC_cv=v, lo=ci[0], hi=ci[1])
        else:
            v = r2(y, pr); s['_p'] = pr; ci = cluster_boot(s, lambda q: r2(q[target].values, q['_p'].values), 500)
            res[name] = dict(R2_cv=v, lo=ci[0], hi=ci[1])
    return res, s


def leave_arm_out(df, target, feats, kind, train_arms=('B', 'C'), test_arms=('D', 'E', 'H')):
    s = df.dropna(subset=[target] + list(feats))
    tr = s[s.arm.isin(train_arms)]; out = {}
    X = tr[list(feats)].values.astype(float); mu, sd = std_fit(X)
    if kind == 'auc': b = logit_fit((X - mu) / sd, tr[target].values, 1.0)
    else:
        A = np.c_[np.ones(len(tr)), (X - mu) / sd]; lam = np.eye(A.shape[1]); lam[0, 0] = 0
        c = np.linalg.solve(A.T @ A + lam, A.T @ tr[target].values)
    for a in test_arms:
        te = s[s.arm == a]
        if len(te) == 0: continue
        Xt = (te[list(feats)].values - mu) / sd
        if kind == 'auc':
            out[a] = dict(n=len(te), AUC=auc(te[target], logit_pred(b, Xt)))
        else:
            p = np.c_[np.ones(len(te)), Xt] @ c
            out[a] = dict(n=len(te), R2=r2(te[target].values, p), MAE=float(np.abs(te[target].values - p).mean()))
    return out
