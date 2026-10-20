"""Analysis of the local MCNS core contrast (+ FAFB reference recomputed the same way).

FII (preregistered matrix version) = Pearson r of vec(R_condition) vs vec(R_baseline-mean) over entries >0 in either.
Variants: FII_resp = brain readouts x RESPONSIVE probes only (a probe is responsive if >=3 brain readouts respond >8 Hz
in the baseline mean); FII_all = all 7 probes x brain readouts; FII_vnc = VNC motor/efferent readouts, responsive probes.
Noise ceiling = median leave-one-out single-baseline-seed r; FII_norm = FII / ceiling.
Outcome classes as in dynamics/analyze.py: runaway (late obs rate >0.5 Hz/neuron or FSA >0.05) > brain_death
(FSA<0.01 and no readout >8 Hz) > recovered (FII_norm>=0.8 and Reach>=0.8 x baseline) > structural_loss (dead>=10 %) > impaired.
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
import os, json, glob
import numpy as np, pandas as pd

R = os.path.join(_P.OUT, 'followup', 'replication')
FAFB_RAW = os.path.join(_P.OUT, 'dynamics', 'raw')
PROBES = ["sugar", "bitter", "JO_aud", "JO_wg", "LC4", "LC9", "turnL"]


def fii(Rb, Rc):
    U = (Rb > 0) | (Rc > 0)
    if U.sum() < 3 or not (Rc[U] != 0).any():
        return 0.0
    x, y = Rb[U].ravel(), Rc[U].ravel()
    if x.std() == 0 or y.std() == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def ro_masks(ds):
    t = pd.read_csv(os.path.join(R, "cache", f"cns_{ds}_neurons.csv.gz"), low_memory=False)
    ty = t.primary_type.astype(str); sc = t.super_class.astype(str); reg = t.region.fillna("NA").astype(str)
    mn9 = np.flatnonzero(ty == "MN9")
    brain = np.union1d(np.union1d(np.flatnonzero(sc == "descending"), np.flatnonzero((sc == "motor") & (reg != "VNC"))), mn9)
    vnc = np.flatnonzero(sc.isin(["motor", "efferent"]) & (reg == "VNC"))
    allr = np.union1d(brain, vnc)
    return np.isin(allr, brain), np.isin(allr, vnc)


def load_cns(res_dir, label):
    js, Rm = {}, {}
    for f in glob.glob(os.path.join(res_dir, label, "*", "*.json")):
        o = json.load(open(f)); js[o["id"]] = o
        Rm[o["id"]] = np.load(f[:-5] + ".npz")["R"].astype(np.float32)
    return js, Rm


def load_fafb(ids_prefix=("base_", "C_D30_naive_s", "C_D30_block20_s", "C_D180_naive_s", "C_D180_block20_s")):
    rows = pd.read_csv(os.path.join(_P.OUT, 'dynamics', 'tables', 'all_runs.csv'))
    js, Rm = {}, {}
    for _, r in rows.iterrows():
        if not r["id"].startswith(ids_prefix):
            continue
        f = os.path.join(FAFB_RAW, r["arm"], r["id"] + ".npz")
        if os.path.exists(f):
            Rm[r["id"]] = np.load(f)["R"].astype(np.float32); js[r["id"]] = r.to_dict()
    return js, Rm


def responsive(Rbase_mean, mask, thr=8.0, k=3):
    return np.array([(Rbase_mean[j, mask] > thr).sum() >= k for j in range(Rbase_mean.shape[0])])


def score(js, Rm, isb, isv, ds):
    base = sorted(k for k in js if k.startswith("base"))
    Rb = np.mean([Rm[k] for k in base], 0)
    resp_b = responsive(Rb, isb); resp_v = responsive(Rb, isv) if isv is not None and isv.any() else None
    sets = {"resp": (resp_b, isb), "all": (np.ones(len(PROBES), bool), isb)}
    if resp_v is not None and resp_v.any():
        sets["vnc"] = (resp_v, isv)
    ceil = {}
    for nm, (pm, cm) in sets.items():
        rs = [fii(np.mean([Rm[j] for j in base if j != k], 0)[np.ix_(pm, cm)], Rm[k][np.ix_(pm, cm)]) for k in base]
        ceil[nm] = float(np.median(rs))
    reach_key = "reach_brain" if ds != "FAFB" else "reach"
    rb = np.mean([js[k][reach_key] for k in base])
    out = []
    for k, o in js.items():
        d = dict(ds=ds, id=k, name=("base" if k.startswith("base") else o["name"]), seed=o["seed"])
        for nm, (pm, cm) in sets.items():
            d[f"FII_{nm}"] = fii(Rb[np.ix_(pm, cm)], Rm[k][np.ix_(pm, cm)]); d[f"FIIn_{nm}"] = d[f"FII_{nm}"] / ceil[nm]
        for key in ["dead_frac", "FSA", "obs_late_Hz", "obs_late_Hz_OL", "obs_late_Hz_CB", "obs_late_Hz_VNC", "surge_peak_Hz",
                    "SD_onset_s", "K_peak", "repol_s_after_reperf", "ignite_s_OL", "ignite_s_CB", "ignite_s_VNC", "wall_s"]:
            v = o.get(key); d[key] = None if (v is None or (isinstance(v, float) and np.isnan(v))) else v
        d["reach"] = o[reach_key]
        for rg in ["OL", "CB", "VNC"]:
            d[f"FSA_{rg}"] = (o.get("FSA_region") or {}).get(rg) if isinstance(o.get("FSA_region"), dict) else None
        ol = d["obs_late_Hz"] or 0.0
        runaway = (ol > 0.5) or (d["FSA"] > 0.05)
        rmax = max([v for v in [o.get("readout_resp_max_brain"), o.get("readout_resp_max_vnc"), o.get("readout_resp_max")] if v is not None and v == v] or [0])
        if runaway:
            oc = "runaway"
        elif d["FSA"] < 0.01 and rmax <= 8:
            oc = "brain_death"
        elif d["FIIn_resp"] >= 0.8 and d["reach"] >= 0.8 * rb:
            oc = "recovered"
        elif d["dead_frac"] >= 0.10:
            oc = "structural_loss"
        else:
            oc = "impaired"
        d["outcome"] = oc
        out.append(d)
    return pd.DataFrame(out), ceil, dict(resp_probes=[p for p, b in zip(PROBES, resp_b) if b],
                                          vnc_resp_probes=[p for p, b in zip(PROBES, resp_v)] if resp_v is not None else [])


def boot_ci(x, n=2000, seed=0):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed); b = rng.choice(x, (n, len(x))).mean(1)
    return (float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)))


if __name__ == "__main__":
    isb, isv = ro_masks("MCNS")
    jsM, RmM = load_cns(os.path.join(R, "raw", "local_v3", "res"), "MCNS")
    dM, cM, pM = score(jsM, RmM, isb, isv, "MCNS")
    jsF, RmF = load_fafb()
    nF = RmF[next(iter(RmF))].shape[1]
    dF, cF, pF = score(jsF, RmF, np.ones(nF, bool), None, "FAFB")
    dF["name"] = dF["name"].str.replace("block20", "reset")
    df = pd.concat([dM, dF], ignore_index=True)
    os.makedirs(os.path.join(R, "tables"), exist_ok=True)
    df.to_csv(os.path.join(R, "tables", "local_core_runs.csv"), index=False)
    rows = []
    for (ds, nm), g in df.groupby(["ds", "name"]):
        r = dict(ds=ds, condition=nm, n=len(g))
        for c in ["FIIn_resp", "FIIn_all", "FIIn_vnc", "dead_frac", "FSA", "obs_late_Hz", "obs_late_Hz_OL", "obs_late_Hz_CB",
                  "obs_late_Hz_VNC", "surge_peak_Hz", "SD_onset_s", "K_peak"]:
            if c in g and g[c].notna().any():
                r[c] = float(g[c].astype(float).mean()); lo, hi = boot_ci(g[c].astype(float))
                r[c + "_lo"] = lo; r[c + "_hi"] = hi
        r["outcomes"] = json.dumps(g.outcome.value_counts().to_dict())
        rows.append(r)
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(R, "tables", "local_core_summary.csv"), index=False)
    meta = dict(ceiling_MCNS=cM, ceiling_FAFB=cF, probes_MCNS=pM, probes_FAFB=pF)
    json.dump(meta, open(os.path.join(R, "tables", "local_core_meta.json"), "w"), indent=1)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print(json.dumps(meta))
    print(S[["ds", "condition", "n", "FIIn_resp", "FIIn_resp_lo", "FIIn_resp_hi", "FIIn_all", "dead_frac", "FSA", "obs_late_Hz",
             "obs_late_Hz_OL", "obs_late_Hz_CB", "obs_late_Hz_VNC", "outcomes"]].round(3).to_string())
