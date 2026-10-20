"""Dynamics-track analysis: load harvested results (hpc/<job>/out/results/<arm>/*.json|npz), recompute FII uniformly
from the stored readout response matrices, classify outcomes, test H1/H3/H4/H5, write tables.
usage: python analyze.py <results_dir> [<results_dir> ...]   (writes into $RESEARCH_ROOT\\dynamics\\tables)"""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src",):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import os, sys, json, glob, math
import numpy as np, pandas as pd
from scipy.optimize import curve_fit

OUTD = os.environ.get("DYN_OUT", os.path.join(_P.OUT, 'dynamics'))
TAB = os.path.join(OUTD, "tables"); os.makedirs(TAB, exist_ok=True)
RNG = np.random.default_rng(0)


def fii_matrix(Rb, Rc):
    U = (Rb > 0) | (Rc > 0)
    if U.sum() < 3 or not (Rc[U] != 0).any():
        return 0.0
    x, y = Rb[U].ravel(), Rc[U].ravel()
    if x.std() == 0 or y.std() == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def load(dirs):
    rows, Rs, TR = [], {}, {}
    for d in dirs:
        for f in glob.glob(os.path.join(d, "*", "*.json")):
            j = json.load(open(f))
            cid = j["id"]
            if cid in Rs:
                continue
            z = np.load(f[:-5] + ".npz")
            Rs[cid] = z["R"].astype(np.float32)
            TR[cid] = {k[6:]: z[k] for k in z.files if k.startswith("trace_")}
            if "probe_rate" in z.files:
                TR[cid]["probe_rate"] = z["probe_rate"]
            rows.append(j)
    df = pd.DataFrame(rows).set_index("id", drop=False)
    return df, Rs, TR


def boot_ci(x, n=2000, f=np.mean):
    x = np.asarray([v for v in x if v is not None and np.isfinite(v)], float)
    if x.size == 0:
        return (np.nan, np.nan, np.nan)
    if x.size == 1:
        return (float(x[0]), np.nan, np.nan)
    b = [f(RNG.choice(x, x.size, replace=True)) for _ in range(n)]
    return (float(f(x)), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)))


def main(dirs):
    df, Rs, TR = load(dirs)
    base_ids = sorted(i for i in df.index if df.loc[i, "arm"] == "base")
    Rb = np.mean([Rs[i] for i in base_ids], 0)
    # noise ceiling: leave-one-out single-seed vs mean of remaining baseline seeds
    loo = [fii_matrix(np.mean([Rs[k] for k in base_ids if k != i], 0), Rs[i]) for i in base_ids]
    ceil = float(np.median(loo))
    base_evoked = float(df.loc[base_ids, "n_evoked"].mean()); base_reach = float(df.loc[base_ids, "reach"].mean())
    df["FII"] = [fii_matrix(Rb, Rs[i]) if i not in base_ids else loo[base_ids.index(i)] for i in df.index]
    df["FII_norm"] = df.FII / ceil
    # rank-based variant (Spearman)
    from scipy.stats import spearmanr
    def sp_(i):
        U = (Rb > 0) | (Rs[i] > 0)
        if U.sum() < 3 or not (Rs[i][U] != 0).any():
            return 0.0
        r = spearmanr(Rb[U].ravel(), Rs[i][U].ravel()).correlation
        return float(r) if np.isfinite(r) else 0.0
    df["FII_spearman"] = [sp_(i) for i in df.index]
    df["Act_rel"] = df.n_evoked / base_evoked
    df["Reach_rel"] = df.reach / base_reach
    df["mag_ratio"] = [float(np.abs(Rs[i]).sum() / np.abs(Rb).sum()) for i in df.index]
    # shuffled-label negative control & gain control (H5 controls)
    perm = RNG.permutation(Rb.shape[1])
    ctrl = dict(noise_ceiling=ceil, loo=loo, shuffled_readouts=fii_matrix(Rb, Rb[:, perm]), half_gain=fii_matrix(Rb, 0.5 * Rb))

    def outcome(r):
        late = r.get("obs_late_Hz"); late = 0.0 if late is None or not np.isfinite(late) else late
        if late > 0.5 or r.FSA > 0.05:
            return "runaway"
        if r.FSA < 0.01 and r.n_readout_resp == 0:
            return "brain_death"
        if r.FII_norm >= 0.8 and r.Reach_rel >= 0.8:
            return "recovered"
        if r.dead_frac >= 0.10:
            return "structural_loss"
        return "impaired"
    df["outcome"] = [outcome(df.loc[i]) for i in df.index]
    df.drop(columns=[c for c in ["named_resp", "reach_per_probe", "n_evoked_per_probe", "probe_names", "FII_neuron_per_probe"] if c in df], errors="ignore").to_csv(os.path.join(TAB, "all_runs.csv"), index=False)

    # per-condition summary
    g = []
    for (arm, name), sub in df.groupby(["arm", "name"]):
        row = dict(arm=arm, name=name, n=len(sub))
        for k in ["FII_norm", "FII_neuron", "dead_frac", "SD_onset_s", "K_peak", "surge_peak_Hz", "obs_late_Hz", "FSA", "Act_rel", "Reach_rel", "mag_ratio"]:
            if k in sub:
                m, lo, hi = boot_ci(sub[k].values.astype(float) if sub[k].dtype != object else [v for v in sub[k].values if v is not None])
                row[k] = m; row[k + "_lo"] = lo; row[k + "_hi"] = hi
        oc = sub.outcome.value_counts()
        row["outcome_mode"] = oc.index[0]; row["outcomes"] = ";".join(f"{a}:{b}" for a, b in oc.items())
        for k in ["D", "T_ins", "reperf", "insult"]:
            if k in sub:
                row[k] = sub[k].iloc[0]
        g.append(row)
    S = pd.DataFrame(g); S.to_csv(os.path.join(TAB, "condition_summary.csv"), index=False)
    res = dict(controls=ctrl, base=dict(n=len(base_ids), evoked=base_evoked, reach=base_reach), n_runs=int(len(df)),
               n_conditions=int(len(S)), outcomes=df.outcome.value_counts().to_dict())

    # ---------------- A: calibration
    A = df[df.arm == "A"]
    res["A"] = {n: dict(SD_onset=boot_ci(s.SD_onset_s.dropna().values.astype(float)), dead=boot_ci(s.dead_frac.values),
                        FII_norm=boot_ci(s.FII_norm.values), outcomes=s.outcome.value_counts().to_dict())
                for n, s in A.groupby("name")}
    # LT50 extrapolation from J slope for arrest runs
    lt = {}
    for n, s in A[A.name.str.startswith("arrest")].groupby("name"):
        vals = []
        for i in s.index:
            tr = TR[i]; t = tr["t"]; J = tr["J"]; mk = df.loc[i, "marks"]
            w = (t >= mk["insult_start"] + 60) & (t < mk["reperf_start"])
            if w.sum() > 5:
                slope = np.polyfit(t[w] - mk["insult_start"], J[w], 1)[0]
                j50 = float(df.loc[i, "mx"].get("J50_s", 600.0)); hmax = 1 / 30.0
                tt = np.linspace(0, 200000, 200001); JJ = np.maximum(slope * tt, 0)
                H = np.cumsum(hmax * JJ ** 4 / (JJ ** 4 + j50 ** 4)) * (tt[1] - tt[0])
                k = np.searchsorted(H, math.log(2))
                vals.append(dict(J_slope_per_s=float(slope), LT50_s=float(tt[k]) if k < tt.size else None,
                                 dead_end=float(df.loc[i, "dead_frac"]), e_end_insult=float(tr["e"][w][-1])))
        lt[n] = vals
    res["A_LT50_extrapolation"] = lt

    # ---------------- B: brain-death map and H4
    B = df[(df.arm == "B") & (df.insult == "anoxia") & ~df.name.str.startswith("q10")]
    mapt = B.groupby(["reperf", "T_ins", "D"]).agg(FII_norm=("FII_norm", "mean"), dead=("dead_frac", "mean"),
                                                    SD=("SD_onset_s", lambda x: np.nanmean([v for v in x if v is not None]) if any(v is not None for v in x) else np.nan),
                                                    outcome=("outcome", lambda x: x.value_counts().index[0]),
                                                    rec_frac=("outcome", lambda x: float((x == "recovered").mean()))).reset_index()
    mapt.to_csv(os.path.join(TAB, "B_brain_death_map.csv"), index=False)

    def dmax(sub):
        """largest D (log-interp) at which mean FII_norm >= 0.8 and Reach >=0.8 (recovered); per-seed bootstrap."""
        Ds = np.array(sorted(sub.D.unique()), float)
        def crit(ss):
            m = np.array([ss[ss.D == d].FII_norm.mean() if ((ss.D == d).any()) else np.nan for d in Ds])
            ok = m >= 0.8
            if not ok[0]:
                return Ds[0] / 2
            if ok.all():
                return np.inf
            k = np.argmax(~ok)          # first failing D
            d0, d1 = Ds[k - 1], Ds[k]; m0, m1 = m[k - 1], m[k]
            fr = (m0 - 0.8) / max(m0 - m1, 1e-9)
            return float(np.exp(np.log(d0) + fr * (np.log(d1) - np.log(d0))))
        est = crit(sub)
        bs = []
        seeds = sub.seed.unique()
        for _ in range(1000):
            pick = []
            for d in Ds:
                ss = sub[sub.D == d]
                pick.append(ss.iloc[RNG.integers(0, len(ss), len(ss))])
            bs.append(crit(pd.concat(pick)))
        bs = np.array(bs)
        return est, bs
    H4 = {}
    for (lab, q, sel) in [("cal_q2.5", 2.5, B[B.reperf == "reset"])] + \
            [(f"q{q}", q, df[(df.arm == "B") & df.name.str.startswith(f"q10_{q}_")]) for q in (2.0, 3.0)]:
        out = {}
        for T in (17.0, 22.0, 29.0, 33.0):
            sub = sel[sel.T_ins == T]
            if len(sub):
                e, bs = dmax(sub); out[T] = (e, bs)
        if 17.0 in out and 29.0 in out:
            e17, b17 = out[17.0]; e29, b29 = out[29.0]
            def q10w(a, b):
                return (a / b) ** (10.0 / 12.0) if np.isfinite(a) and np.isfinite(b) and b > 0 else np.nan
            qs = np.array([q10w(a, b) for a, b in zip(b17, b29)])
            ratio = np.array([a / b if np.isfinite(a) and np.isfinite(b) else np.nan for a, b in zip(b17, b29)])
            H4[lab] = dict(model_q10=q, Dmax={str(T): out[T][0] for T in out},
                           Dmax_CI={str(T): [float(np.nanpercentile(np.where(np.isfinite(out[T][1]), out[T][1], np.nan), 2.5)) if np.isfinite(out[T][1]).any() else None,
                                             float(np.nanpercentile(np.where(np.isfinite(out[T][1]), out[T][1], np.nan), 97.5)) if np.isfinite(out[T][1]).any() else None] for T in out},
                           Q10_window=q10w(e17, e29), Q10_window_CI=[float(np.nanpercentile(qs, 2.5)), float(np.nanpercentile(qs, 97.5))] if np.isfinite(qs).any() else None,
                           ratio_CI=[float(np.nanpercentile(ratio, 2.5)), float(np.nanpercentile(ratio, 97.5))] if np.isfinite(ratio).any() else None,
                           censored_frac_17=float(np.mean(~np.isfinite(b17))), censored_frac_29=float(np.mean(~np.isfinite(b29))))
    res["H4"] = H4
    # other B insults
    res["B_other"] = {n: dict(FII_norm=boot_ci(s.FII_norm.values), dead=boot_ci(s.dead_frac.values), K_peak=boot_ci(s.K_peak.values),
                              SD_onset=boot_ci([v for v in s.SD_onset_s.values if v is not None]), outcomes=s.outcome.value_counts().to_dict())
                      for n, s in df[(df.arm == "B") & (df.insult != "anoxia")].groupby("name")}

    # ---------------- H1: FII vs damage fraction (logistic vs linear vs exponential)
    pool = df[(df.arm.isin(["B", "C", "H"])) & (df.reperf.isin(["reset"]) | df.name.str.contains("reset")) & (df.outcome != "runaway")]
    f = pool.dead_frac.values.astype(float); y = pool.FII_norm.values.astype(float)
    def aic(rss, n, k):
        return n * math.log(max(rss, 1e-12) / n) + 2 * k
    H1 = dict(n=int(len(f)), f_range=[float(f.min()), float(f.max())] if len(f) else None)
    try:
        lin = np.polyfit(f, y, 1); rss_lin = float(((np.polyval(lin, f) - y) ** 2).sum())
        ex = lambda x, a, b: a * np.exp(-b * x)
        pe, _ = curve_fit(ex, f, y, p0=[1, 5], maxfev=20000); rss_ex = float(((ex(f, *pe) - y) ** 2).sum())
        lg = lambda x, top, bot, fc, s: bot + (top - bot) / (1 + np.exp((x - fc) / s))
        pl, _ = curve_fit(lg, f, y, p0=[1, 0, np.median(f) if f.max() > 0 else 0.1, 0.02], maxfev=50000,
                          bounds=([0, -0.5, 0, 1e-4], [1.5, 1, 1, 1]))
        rss_lg = float(((lg(f, *pl) - y) ** 2).sum()); n = len(f)
        A_ = dict(linear=aic(rss_lin, n, 2), exponential=aic(rss_ex, n, 2), logistic=aic(rss_lg, n, 4))
        w = 2 * math.log(9) * pl[3]
        H1.update(AIC=A_, dAIC_logistic_vs_best_alt=A_["logistic"] - min(A_["linear"], A_["exponential"]),
                  logistic=dict(top=pl[0], bottom=pl[1], f_c=pl[2], scale=pl[3], width_10_90=w))
        # recovered at f >= f_c+0.10 under any protocol?
        allp = df[(df.arm.isin(["B", "C", "H"]))]
        H1["n_recovered_beyond_fc_plus_0.10"] = int(((allp.dead_frac >= pl[2] + 0.10) & (allp.outcome == "recovered")).sum())
        H1["max_dead_frac_recovered"] = float(allp[allp.outcome == "recovered"].dead_frac.max()) if (allp.outcome == "recovered").any() else None
    except Exception as e:
        H1["error"] = repr(e)
    res["H1"] = H1

    # ---------------- C: resuscitation ranking and H3
    C = df[df.arm == "C"].copy(); C["prot"] = C.name.str.split("_", n=1).str[1]
    rank = []
    for p, s in C.groupby("prot"):
        row = dict(protocol=p)
        for D in (30, 180, 360):
            ss = s[s.D == D]
            m, lo, hi = boot_ci(ss.FII_norm.values); row[f"FII_norm_D{D}"] = m; row[f"lo_D{D}"] = lo; row[f"hi_D{D}"] = hi
            row[f"dead_D{D}"] = float(ss.dead_frac.mean()); row[f"rec_frac_D{D}"] = float((ss.outcome == "recovered").mean())
            row[f"surge_D{D}"] = float(np.nanmean([v for v in ss.surge_peak_Hz.values if v is not None])) if len(ss) else np.nan
            row[f"outcome_D{D}"] = ss.outcome.value_counts().index[0] if len(ss) else None
        row["FII_norm_mean"] = np.nanmean([row[f"FII_norm_D{D}"] for D in (30, 180, 360)])
        rank.append(row)
    RK = pd.DataFrame(rank).sort_values("FII_norm_mean", ascending=False); RK.to_csv(os.path.join(TAB, "C_resuscitation_ranking.csv"), index=False)
    res["C_ranking"] = RK[["protocol", "FII_norm_mean", "FII_norm_D30", "FII_norm_D180", "FII_norm_D360", "dead_D360", "outcome_D30", "outcome_D180", "outcome_D360"]].to_dict("records")

    def delta(pa, pb, D, key="FII_norm"):
        a = C[(C.prot == pa) & (C.D == D)].sort_values("seed")[key].values.astype(float)
        b = C[(C.prot == pb) & (C.D == D)].sort_values("seed")[key].values.astype(float)
        if len(a) == 0 or len(b) == 0:
            return None
        d = [np.mean(RNG.choice(a, len(a))) - np.mean(RNG.choice(b, len(b))) for _ in range(2000)]
        return dict(delta=float(a.mean() - b.mean()), lo=float(np.percentile(d, 2.5)), hi=float(np.percentile(d, 97.5)))
    H3 = {}
    for staged in ["ramp10", "ramp30", "ramp60", "pulsed3x10", "pulsed6x5"]:
        H3[staged] = {D: dict(dFII=delta(staged, "naive", D), dSurv=delta(staged, "naive", D, "dead_frac"),
                              surge=(float(np.nanmean([v for v in C[(C.prot == staged) & (C.D == D)].surge_peak_Hz if v is not None])),
                                     float(np.nanmean([v for v in C[(C.prot == "naive") & (C.D == D)].surge_peak_Hz if v is not None]))))
                      for D in (30, 180, 360)}
    H3["negative_control_noinj"] = {D: delta("ramp30_noinj", "naive_noinj", D) for D in (30, 180, 360)}
    H3["reset_vs_naive"] = {D: delta("block20", "naive", D) for D in (30, 180, 360)}
    res["H3"] = H3

    # ---------------- H5: D50 ordering across insult duration (B anoxia, 22C, reset and naive pooled per protocol)
    H5 = {}
    for rp in ("reset", "naive"):
        sub = B[(B.T_ins == 22.0) & (B.reperf == rp)]
        Ds = np.array(sorted(sub.D.unique()), float)
        def d50(col, ss):
            m = np.array([ss[ss.D == d][col].mean() for d in Ds]); top, bot = m[0], np.nanmin(m)
            half = bot + 0.5 * (top - bot)
            if top - bot < 0.05:
                return np.inf
            k = np.argmax(m <= half)
            if k == 0:
                return Ds[0]
            d0, d1, m0, m1 = Ds[k - 1], Ds[k], m[k - 1], m[k]
            fr = (m0 - half) / max(m0 - m1, 1e-9)
            return float(np.exp(np.log(d0) + fr * (np.log(d1) - np.log(d0))))
        out = {}
        for col in ["FII_norm", "Act_rel", "Reach_rel"]:
            bs = []
            for _ in range(500):
                pick = pd.concat([sub[sub.D == d].sample(frac=1, replace=True, random_state=int(RNG.integers(1e9))) for d in Ds])
                bs.append(d50(col, pick))
            bs = np.array(bs)
            out[col] = dict(D50=d50(col, sub), CI=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))], censored=float(np.mean(~np.isfinite(bs))))
        H5[rp] = out
    grid_runs = df[df.arm.isin(["B", "C", "D", "E", "H"])]
    rbd = (grid_runs.Act_rel >= 0.8) & (grid_runs.Reach_rel >= 0.8) & (grid_runs.FII_norm <= 0.5)
    H5["resuscitated_but_different_frac_runs"] = float(rbd.mean())
    H5["rbd_frac_runs_excluding_runaway"] = float((rbd & (grid_runs.outcome != "runaway")).mean())
    cs = S[S.arm.isin(["B", "C", "D", "E", "H"])]
    H5["rbd_frac_conditions"] = float(((cs.Act_rel >= 0.8) & (cs.Reach_rel >= 0.8) & (cs.FII_norm <= 0.5)).mean())
    res["H5"] = H5

    # ---------------- D, E, H
    res["D"] = {n: dict(FII_norm=boot_ci(s.FII_norm.values), dead=boot_ci(s.dead_frac.values), SD_onset=boot_ci([v for v in s.SD_onset_s if v is not None]),
                        outcomes=s.outcome.value_counts().to_dict()) for n, s in df[df.arm == "D"].groupby("name")}
    res["E"] = {n: dict(FII_norm=boot_ci(s.FII_norm.values), dead=boot_ci(s.dead_frac.values), SD_onset=boot_ci([v for v in s.SD_onset_s if v is not None]),
                        syn_lost=float(s.syn_lost_frac.mean()) if "syn_lost_frac" in s else None, Act_rel=boot_ci(s.Act_rel.values),
                        dep_max=boot_ci([v for v in s.dep_max if v is not None]) if "dep_max" in s else None,
                        outcomes=s.outcome.value_counts().to_dict()) for n, s in df[df.arm == "E"].groupby("name")}
    Hh = {}
    for n, s in df[df.arm == "H"].groupby("name"):
        Hh[n] = dict(FII_norm=boot_ci(s.FII_norm.values), dead=boot_ci(s.dead_frac.values), n_prot=int(s.n_protected.iloc[0]),
                     outcomes=s.outcome.value_counts().to_dict())
    for rp in ("naive", "reset"):
        ref = C[(C.prot == ("naive" if rp == "naive" else "block20")) & (C.D == 360)]
        r0 = float(ref.FII_norm.mean()) if len(ref) else np.nan; d0 = float(ref.dead_frac.mean()) if len(ref) else np.nan
        for n in list(Hh):
            if n.endswith(rp):
                Hh[n]["ref_unprotected_FII_norm"] = r0; Hh[n]["ref_unprotected_dead"] = d0
                Hh[n]["FII_gain_per_1000_protected"] = (Hh[n]["FII_norm"][0] - r0) / Hh[n]["n_prot"] * 1000
                Hh[n]["deaths_prevented_per_protected"] = (d0 - Hh[n]["dead"][0]) * 139255 / Hh[n]["n_prot"]
    res["H"] = Hh
    json.dump(res, open(os.path.join(TAB, "analysis_results.json"), "w"), indent=1, default=lambda x: (x.tolist() if hasattr(x, "tolist") else (None if x is None else str(x))))
    df.to_pickle(os.path.join(TAB, "all_runs.pkl"))
    import pickle; pickle.dump(dict(Rb=Rb, TR=TR), open(os.path.join(TAB, "traces.pkl"), "wb"))
    print(json.dumps(dict(n_runs=res["n_runs"], n_conditions=res["n_conditions"], outcomes=res["outcomes"], ceil=ceil), default=str))
    return df, S, res


if __name__ == "__main__":
    main(sys.argv[1:])
