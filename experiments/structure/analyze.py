"""Fit Reach(f) curves (logistic-4 vs linear vs exponential-3), D50, susceptibility, bootstrap ratios.
python analyze.py  -> tables/fits_all.csv, tables/d50_all.csv, tables/ratios_all.csv, tables/susceptibility_all.csv
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
import os, sys, glob, warnings
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
warnings.filterwarnings('ignore')
ROOT = os.path.join(_P.OUT, 'structure')
T = os.path.join(ROOT, 'tables')
METRICS = ['reachN_L6', 'reachN_L4', 'reachN_Linf', 'reachN_L6_DN', 'reachN_L6_MN', 'reachG_L6', 'reachGloose_L6',
           'redund_L6', 'influence_L6', 'gscc']
PRIMARY = 'reachN_L6'
NBOOT = int(os.environ.get('NBOOT', 2000))
W_FACTOR = 2 * np.log(9)  # 10-90% width = 4.394 * s


def logistic(f, a, d, c, s):
    return d + (a - d) / (1 + np.exp(np.clip((f - c) / s, -50, 50)))


def expo(f, a, b, c):
    return a * np.exp(-b * f) + c


def aic(rss, n, k):
    return n * np.log(max(rss, 1e-300) / n) + 2 * k


def fit_logistic(f, y, p0=None):
    best = None
    inits = [p0] if p0 is not None else []
    inits += [(1, 0, c0, s0) for c0 in (0.05, 0.2, 0.4, 0.6, 0.8, 1.0) for s0 in (0.02, 0.08, 0.2)]
    for p in inits:
        try:
            popt, _ = curve_fit(logistic, f, y, p0=p, bounds=([0.5, -0.2, -0.5, 1e-3], [1.5, 1.0, 2.0, 2.0]), maxfev=4000)
            rss = float(((y - logistic(f, *popt)) ** 2).sum())
            if best is None or rss < best[1]: best = (popt, rss)
        except Exception:
            pass
        if p0 is not None and best is not None: break
    return best


def fit_all(f, y):
    n = y.size
    lo = fit_logistic(f, y)
    b = np.polyfit(f, y, 1); rss_lin = float(((y - np.polyval(b, f)) ** 2).sum())
    best_e = None
    for b0 in (0.5, 2, 5, 15):
        try:
            pe, _ = curve_fit(expo, f, y, p0=(1, b0, 0), bounds=([0, 0, -1], [2, 100, 1]), maxfev=4000)
            r = float(((y - expo(f, *pe)) ** 2).sum())
            if best_e is None or r < best_e[1]: best_e = (pe, r)
        except Exception:
            pass
    a_log = aic(lo[1], n, 4); a_lin = aic(rss_lin, n, 2); a_exp = aic(best_e[1], n, 3) if best_e else np.inf
    a, d, c, s = lo[0]
    return dict(n=n, f_c=c, s=s, width_10_90=W_FACTOR * s, upper=a, lower=d, AIC_log=a_log, AIC_lin=a_lin, AIC_exp=a_exp,
                dAIC=a_log - min(a_lin, a_exp), best_alt='linear' if a_lin <= a_exp else 'exponential',
                lin_slope=b[0], exp_rate=best_e[0][1] if best_e else np.nan, rss_log=lo[1]), lo[0]


def d50(fv, yv):
    """first f at which normalised curve falls to <=0.5 (linear interpolation); nan if never (censored)."""
    idx = np.where(yv <= 0.5)[0]
    if idx.size == 0: return np.nan
    i = idx[0]
    if i == 0: return fv[0]
    f0, f1, y0, y1 = fv[i - 1], fv[i], yv[i - 1], yv[i]
    return f0 + (0.5 - y0) * (f1 - f0) / (y1 - y0) if y1 != y0 else f1


def h1_category(dAIC, w, fc):
    if not (0 <= fc <= 0.9):
        return 'no transition within f<=0.9 (logistic midpoint extrapolated)'
    if dAIC <= -10 and w <= 0.10: return 'supported (structural criteria)'
    if w > 0.25 or dAIC > -2: return 'falsified'
    return 'indeterminate (graded transition)'


def main():
    files = [p for p in glob.glob(os.path.join(T, '*_curves*.csv'))]
    if not files: return
    cur = pd.concat([pd.read_csv(p) for p in files], ignore_index=True)
    cur = cur.drop_duplicates(['dataset', 'graph', 'regime', 'order', 'f'])
    # normalise by f=0 of each curve
    for m in METRICS:
        b = cur[cur.f == 0].set_index(['dataset', 'graph', 'regime', 'order'])[m]
        key = pd.MultiIndex.from_frame(cur[['dataset', 'graph', 'regime', 'order']])
        cur[m + '_n'] = cur[m].values / b.reindex(key).values
    cur.to_csv(os.path.join(T, 'curves_normalised_all.csv.gz'), index=False)
    fits, d50s, sus, params = [], [], [], {}
    for (ds, gr, rg), sub in cur.groupby(['dataset', 'graph', 'regime']):
        for m in METRICS:
            if m + '_n' not in sub or sub[m + '_n'].isna().all(): continue
            s2 = sub.dropna(subset=[m + '_n'])
            f = s2.f.values.astype(float); y = s2[m + '_n'].values.astype(float)
            r, p = fit_all(f, y); params[(ds, gr, rg, m)] = p
            r.update(dataset=ds, graph=gr, regime=rg, metric=m, n_orders=s2.order.nunique(),
                     baseline_raw=float(s2[s2.f == 0][m].mean()))
            r['H1_category'] = h1_category(r['dAIC'], r['width_10_90'], r['f_c'])
            fits.append(r)
            dd = [d50(g_.f.values, g_[m + '_n'].values) for _, g_ in s2.sort_values('f').groupby('order')]
            d50s.append(dict(dataset=ds, graph=gr, regime=rg, metric=m, D50_mean=np.nanmean(dd) if np.isfinite(dd).any() else np.nan,
                             D50_sd=np.nanstd(dd), n_censored=int(np.isnan(dd).sum()), n_orders=len(dd)))
            if m == PRIMARY and s2.order.nunique() > 1:
                v = s2.groupby('f')[m + '_n'].agg(['mean', 'var', 'std']).reset_index()
                fpk = float(v.loc[v['var'].idxmax(), 'f'])
                sus.append(dict(dataset=ds, graph=gr, regime=rg, f_peak_var=fpk, max_var=float(v['var'].max()),
                                f_c=r['f_c'], within_0p05=abs(fpk - r['f_c']) <= 0.05))
    fits = pd.DataFrame(fits); d50s = pd.DataFrame(d50s); sus = pd.DataFrame(sus)
    fits.to_csv(os.path.join(T, 'fits_all.csv'), index=False)
    d50s.to_csv(os.path.join(T, 'd50_all.csv'), index=False)
    sus.to_csv(os.path.join(T, 'susceptibility_all.csv'), index=False)

    # bootstrap ratios targeted / random on primary metric (thr5 graph and null graphs)
    rng = np.random.default_rng(7)
    rat = []
    for (ds, gr), sub in cur.groupby(['dataset', 'graph']):
        if 'random' not in set(sub.regime): continue
        R = sub[sub.regime == 'random']
        for rg in sorted(set(sub.regime) - {'random'}):
            Tt = sub[sub.regime == rg]
            for m in [PRIMARY, 'reachN_Linf', 'reachG_L6', 'redund_L6', 'gscc']:
                pr = params.get((ds, gr, 'random', m)); pt = params.get((ds, gr, rg, m))
                if pr is None or pt is None: continue
                fr = fits.query('dataset==@ds and graph==@gr and regime=="random" and metric==@m').iloc[0]
                ft = fits.query('dataset==@ds and graph==@gr and regime==@rg and metric==@m').iloc[0]
                ro = R.order.unique(); to = Tt.order.unique()
                Rg = {o: g_ for o, g_ in R.groupby('order')}; Tg = {o: g_ for o, g_ in Tt.groupby('order')}
                dR = np.array([d50(g_.sort_values('f').f.values, g_.sort_values('f')[m + '_n'].values) for o, g_ in Rg.items()])
                dT = np.array([d50(g_.sort_values('f').f.values, g_.sort_values('f')[m + '_n'].values) for o, g_ in Tg.items()])
                nb = NBOOT if (m == PRIMARY and gr == 'thr5') else 100
                bc, bd = [], []
                for b in range(nb):
                    sr = rng.choice(ro, ro.size); st = rng.choice(to, to.size)
                    if True:
                        rr = pd.concat([Rg[o] for o in sr]); tt = pd.concat([Tg[o] for o in st])
                        a1 = fit_logistic(rr.f.values, rr[m + '_n'].values, p0=tuple(pr))
                        a2 = fit_logistic(tt.f.values, tt[m + '_n'].values, p0=tuple(pt))
                        if a1 and a2 and a1[0][2] > 0: bc.append(a2[0][2] / a1[0][2])
                    ir = rng.integers(0, dR.size, dR.size); it = rng.integers(0, dT.size, dT.size)
                    with np.errstate(all='ignore'):
                        bd.append(np.nanmean(dT[it]) / np.nanmean(dR[ir]))
                bc = np.array(bc); bd = np.array(bd)
                with np.errstate(all='ignore'):
                    d50ratio = np.nanmean(dT) / np.nanmean(dR)
                rat.append(dict(dataset=ds, graph=gr, regime=rg, metric=m,
                                fc_random=fr.f_c, fc_target=ft.f_c, ratio_fc=ft.f_c / fr.f_c,
                                ratio_fc_lo=np.nanpercentile(bc, 2.5) if bc.size else np.nan,
                                ratio_fc_hi=np.nanpercentile(bc, 97.5) if bc.size else np.nan, n_boot_fc=bc.size,
                                D50_random=np.nanmean(dR), D50_target=np.nanmean(dT), ratio_D50=d50ratio,
                                ratio_D50_lo=np.nanpercentile(bd, 2.5) if np.isfinite(bd).any() else np.nan,
                                ratio_D50_hi=np.nanpercentile(bd, 97.5) if np.isfinite(bd).any() else np.nan,
                                n_random_orders=ro.size, n_target_orders=to.size,
                                n_random_censored=int(np.isnan(dR).sum()), n_target_censored=int(np.isnan(dT).sum()),
                                ratio_D50_upper_bound_if_random_censored=(np.nanmean(dT) / 0.9) if np.isnan(dR).all() else np.nan,
                                random_fc_in_range=bool(0 <= fr.f_c <= 0.9), target_fc_in_range=bool(0 <= ft.f_c <= 0.9)))
    pd.DataFrame(rat).to_csv(os.path.join(T, 'ratios_all.csv'), index=False)
    print('fits', len(fits), 'ratios', len(rat))


if __name__ == '__main__':
    main()
