"""Step 2 runner: univariate + combined prediction of dynamics outcomes from connectome-only indices.
Writes tables/univariate_<target>.csv, tables/model_results.json, tables/index_dataset.csv."""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src",):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import os, sys, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model_index as MI

df, ref = MI.load()
df = df[df.arm.isin(['B', 'C', 'D', 'E', 'H'])].reset_index(drop=True)
df.to_csv(os.path.join(MI.TAB, 'index_dataset.csv'), index=False)
feats = {k: v for k, v in MI.FEATS.items() if k in df.columns and df[k].notna().mean() > 0.95}
res = {'n_runs': int(len(df)), 'n_conditions': int(df.cond.nunique()), 'arms': df.arm.value_counts().to_dict(), 'features': list(feats)}

# targets
T = {}
reset = df[df.proto == 'reset']
T['FII_reset_nonrunaway'] = (reset[reset.outcome != 'runaway'].copy(), 'FII_norm', 'r2')
T['recovered_reset'] = (reset.copy(), 'recovered', 'auc')
nr = df[df.proto == 'noreset']
T['runaway_noreset'] = (nr.copy(), 'runaway', 'auc')
COMB = [f for f in ['alive', 'Reff_mf_known', 'Reff_mf', 'hub1', 'reach_sim', 'gscc_exc', 'EI_syn', 'ecore100'] if f in feats]
for name, (d, y, kind) in T.items():
    u = MI.univariate(d, y, feats, kind)
    u = u.sort_values('AUC_cv_mean' if kind == 'auc' else 'R2_cv', ascending=False)
    u.to_csv(os.path.join(MI.TAB, f'univariate_{name}.csv'), index=False)
    comb, s = MI.combined(d, y, COMB, kind)
    lao = MI.leave_arm_out(d, y, COMB, kind)
    lao1 = {f: MI.leave_arm_out(d, y, [f], kind) for f in u.feature.head(3)}
    res[name] = dict(n=int(len(d)), n_cond=int(d.cond.nunique()), kind=kind, top=u.head(6).to_dict('records'),
                     combined=comb, combined_features=COMB, leave_arm_out_combined=lao, leave_arm_out_top3=lao1)
    print(name, len(d), u.head(6).to_string(), comb, flush=True)

# 5-class outcome (R/A/S/X/I) with protocol indicator; balanced accuracy, grouped CV
d = df.dropna(subset=COMB).reset_index(drop=True)
cls = sorted(d.cls.dropna().unique()); Y = np.stack([(d.cls == c).values for c in cls], 1).astype(float)
X = np.c_[d[COMB].values.astype(float), (d.proto == 'reset').values.astype(float)]
P = np.zeros_like(Y); cnt = np.zeros(len(d))
for tr, te in MI.folds(d, 5, 10):
    mu, sd = MI.std_fit(X[tr]); B = MI.softmax_fit((X[tr] - mu) / sd, Y[tr], 1.0)
    P[te] += MI.softmax_pred(B, (X[te] - mu) / sd); cnt[te] += 1
P /= np.maximum(cnt, 1)[:, None]
pred = P.argmax(1); true = Y.argmax(1)
def bacc(t, p): return float(np.mean([np.mean(p[t == k] == k) for k in np.unique(t)]))
d['_t'] = true; d['_p'] = pred
ci = MI.cluster_boot(d, lambda q: bacc(q['_t'].values, q['_p'].values), 500)
# protocol-only baseline
Xp = (d.proto == 'reset').values.astype(float)[:, None]; P0 = np.zeros_like(Y); c0 = np.zeros(len(d))
for tr, te in MI.folds(d, 5, 10):
    B = MI.softmax_fit(Xp[tr], Y[tr], 1.0); P0[te] += MI.softmax_pred(B, Xp[te]); c0[te] += 1
pred0 = (P0 / np.maximum(c0, 1)[:, None]).argmax(1)
conf = pd.crosstab(pd.Series(np.array(cls)[true], name='true'), pd.Series(np.array(cls)[pred], name='pred'))
conf.to_csv(os.path.join(MI.TAB, 'outcome_class_confusion.csv'))
res['outcome_class'] = dict(classes=cls, n=int(len(d)), balanced_acc=bacc(true, pred), lo=ci[0], hi=ci[1],
                            balanced_acc_protocol_only=bacc(true, pred0),
                            per_class_auc={c: MI.auc(Y[:, k], P[:, k]) for k, c in enumerate(cls)})
print(res['outcome_class'], flush=True)
json.dump(res, open(os.path.join(MI.TAB, 'model_results.json'), 'w'), indent=1, default=float)
print('done')
