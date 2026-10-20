"""Dynamics-track sweep driver (FAFB v783 princeton5, flysim_dyn LIF + metabolic layer).

python sweep.py --shard K --nshard N --workers W [--arms base,A,B,...] [--ids id1,id2] [--list]
Resumable: a condition is skipped if RES_DIR/<arm>/<id>.json exists.  One writer per path.
Per condition: <id>.json (scalars + named readouts + inline FII), <id>.npz (full response matrix, traces,
sparse neuron spike counts per eval epoch, alive mask).  Copies of json/npz-light under ./out/.
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
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import json, time, copy, math, argparse, traceback, glob
import numpy as np
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flysim_dyn as fs

RES = os.environ.get("RES_DIR", os.path.join(_P.OUT, "dynamics", "raw"))
OUT = os.environ.get("OUT_DIR", "out")

# ---------------------------------------------------------------- calibrated parameter set (arm A validates)
SP = 50.0                       # time compression during protocol (slow s per bio s)
SP_COMA = 100.0                 # compression during long anoxic coma plateau (arm A arrest runs only)
CAL = dict(tau_atp_s=30.0, J50_s=600.0, q10_injury=2.5)
REPERF_S = 60.0                 # reperfusion window (slow s)
OBS_S = 20.0                    # free observation after the protocol (slow s), S=1, no interventions
PRE_S = 2.0
PROBE_MS = 250.0; SPONT_MS = 500.0

# ---------------------------------------------------------------- model + populations (loaded once, fork-shared)
M = None; POP = {}; BASE = None; AT = None; DEG = None
OL_R = ["ME_R", "LO_R", "LOP_R", "AME_R", "LA_R"]


def setup():
    global M, POP, AT, DEG
    M = fs.load_model("princeton5")
    m = M; t = m.neurons
    P = {}
    P["sugar"] = (m.idx(fs.SHIU_SUGAR_GRN), 150.0)
    P["bitter"] = (m.select(None, **{"class": "gustatory", "sub_class": "bitter"}), 150.0)
    P["JO_aud"] = (m.select(None, **{"class": "mechanosensory", "sub_class": "auditory"}), 150.0)
    P["JO_wg"] = (m.select(None, **{"class": "mechanosensory", "sub_class": "wind_gravity"}), 150.0)
    P["LC4"] = (m.select("LC4"), 100.0)
    P["LC9"] = (m.select("LC9"), 100.0)
    P["turnL"] = (np.union1d(m.select("LC9", side="left"), m.select("DNa02", side="left")), 180.0)
    POP["probes"] = P
    mn9 = m.idx([fs.SHIU_MN9])
    POP["readouts"] = np.union1d(m.select(None, super_class=["descending", "motor"]), mn9)
    named = {"MN9": mn9}
    for ty in ["DNp01", "DNg100", "DNp09", "DNa02", "DNa01", "DNg97"]:
        for sd in ["left", "right"]:
            named[f"{ty}_{sd[0].upper()}"] = m.select(ty, side=sd)
    named["MDN"] = m.select("MDN")
    POP["named"] = named
    sc = t.super_class.fillna("NA").values
    POP["G_sens"] = np.isin(sc, ["sensory", "sensory_ascending", "visual_projection", "optic"])
    POP["G_cent"] = np.isin(sc, ["central", "visual_centrifugal", "endocrine", "ascending", "NA"])
    POP["G_out"] = np.isin(sc, ["descending", "motor"])
    A = sp.csr_matrix((np.ones(len(m.indices), np.float32), m.indices, m.indptr), shape=(m.N, m.N))
    AT = A.T.tocsr()
    DEG = np.diff(m.indptr) + np.bincount(m.indices, minlength=m.N)
    np_ = m.neuropils
    POP["core_OL_R"] = np.array([n in OL_R for n in np_])
    POP["core_hemi_R"] = np.array([n.endswith("_R") for n in np_])
    POP["base_demand"] = fs.demand_base(m, dict(fs.METAB_DEFAULTS))


def load_base():
    """Baseline (intact) mean probe rates for inline FII (from RES/base/*.npz)."""
    global BASE
    fl = sorted(glob.glob(os.path.join(RES, "base", "*.npz")))
    if len(fl) < 3:
        BASE = None; return
    Rs, rates = [], []
    for f in fl:
        z = np.load(f)
        Rs.append(z["R"]); rates.append(z["probe_rate"].astype(np.float32))
    BASE = dict(R=np.mean(Rs, 0), rates=np.mean(rates, 0), n=len(fl))

# ---------------------------------------------------------------- helpers

def reach(alive):
    """fraction of (probe, readout) pairs with a directed path <=6 hops through surviving neurons."""
    ro = POP["readouts"]; tot = 0; hit = 0
    per = {}
    for name, (idx, _) in POP["probes"].items():
        x = np.zeros(M.N, np.float32); src = idx[alive[idx]]; x[src] = 1.0
        seen = x > 0
        for _ in range(6):
            x = (AT @ x > 0) & alive & ~seen
            seen |= x
            x = x.astype(np.float32)
            if not x.any():
                break
        ok = seen[ro] & alive[ro]
        per[name] = float(ok.mean()); hit += ok.sum(); tot += ro.size
    return float(hit / tot), per


def wrap(f, ts0):
    return (lambda ts, _f=f, _o=ts0: _f(ts - _o)) if callable(f) else f


def phase(st, dur_s, sp_, metab_over, supply=1.0, temp=22.0, iv=None, stim=None, seed=0, inh=1.0, exc=1.0,
          model=None, rec=True):
    m = model or M
    ts0 = (st.t_ms * 1e-3 * sp_) if st is not None else 0.0
    mx = dict(metab_over); mx["speedup"] = sp_
    ivw = {k: (wrap(v, ts0) if callable(v) else v) for k, v in (iv or {}).items()}
    r = fs.run(m, dur_s * 1000.0 / sp_, stim=stim, metab=mx, supply=wrap(supply, ts0), temperature=wrap(temp, ts0),
               interventions=ivw, state=st, groups=None, rec_every_ms=(500.0 / sp_ if rec else 1e9), seed=seed,
               inh_gain=inh, exc_gain=exc)
    return r


class Trace:
    def __init__(self):
        self.t0 = 0.0; self.d = {k: [] for k in ["t", "e", "dep", "blk", "dead", "Km", "Kx", "syn", "J", "pop_t", "pop"]}
        self.marks = {}

    def add(self, r, dur_s, sp_):
        sl = r["slow"]; ts = np.asarray(sl["t_s"]);
        if ts.size:
            tt = self.t0 + (ts - ts[0])
            K = np.asarray(sl["K"])
            for k, v in [("t", tt), ("e", sl["e_mean"]), ("dep", sl["dep_frac"]), ("blk", sl["block_frac"]),
                         ("dead", sl["dead_frac"]), ("Km", K.mean(1)), ("Kx", K.max(1)), ("syn", sl["syn_frac"]),
                         ("J", sl["J_mean"])]:
                self.d[k].append(np.asarray(v, np.float32))
        pop = r["pop_ms"].sum(1).astype(np.float64)
        b = max(1, int(round(500.0 / sp_)))           # 0.5 slow-s bins
        nb = len(pop) // b
        if nb:
            hz = pop[: nb * b].reshape(nb, b).sum(1) / (b * 1e-3) / M.N
            self.d["pop_t"].append((self.t0 + 0.5 * np.arange(nb)).astype(np.float32)); self.d["pop"].append(hz.astype(np.float32))
        self.t0 += dur_s

    def arrays(self):
        return {k: (np.concatenate(v) if v else np.zeros(0, np.float32)) for k, v in self.d.items()}


def piece(pts):
    return fs.piecewise(pts)

# ---------------------------------------------------------------- reperfusion protocols -> list of sub-phases

def reperf_protocol(name, T_rep):
    """returns list of dict(dur, supply, temp, iv, inh, stim_fn, mx, revive)."""
    S = lambda dur, **kw: dict(dict(dur=dur, supply=1.0, temp=T_rep, iv={}, inh=1.0, stim=None, mx={}, revive=False), **kw)
    R = REPERF_S
    blk = lambda on, off: {"block_all": piece([(0, False), (on, True), (on + off, False)])}
    if name == "naive":
        return [S(R)]
    if name.startswith("block"):              # block{d}_del{x}
        d, x = name[5:].split("_del") if "_del" in name else (name[5:], "0")
        return [S(R, iv=blk(float(x), float(d)))]
    if name == "reset":
        return [S(R, iv=blk(0, 20))]
    if name.startswith("voff"):
        v = float(name[4:]); return [S(R, iv={"v_offset": piece([(0, v), (20.0, 0.0)])})]
    if name.startswith("gaba"):
        g = float(name[4:]); return [S(20.0, inh=g), S(R - 20.0)]
    if name == "pulsed3x10":
        return [S(R, supply=piece([(0, 1.0), (10, 0.0), (20, 1.0), (30, 0.0), (40, 1.0)]))]
    if name == "pulsed6x5":
        return [S(R, supply=piece([(i * 5.0, float(i % 2 == 0)) for i in range(11)] + [(55.0, 1.0)]))]
    if name.startswith("ramp"):
        nm = name.split("_")[0]; Tr = float(nm[4:])
        iv = {"injury_scale": 0.0} if name.endswith("noinj") else {}
        return [S(R, supply=fs.ramp([(0, 0.0), (Tr, 1.0), (1e6, 1.0)]), iv=iv)]
    if name == "naive_noinj":
        return [S(R, iv={"injury_scale": 0.0})]
    if name.startswith("cool"):
        Tc = float(name[4:]); return [S(R, temp=piece([(0, Tc), (40.0, 22.0)]))]
    if name == "sdblock":
        return [S(R, iv={"k_rel_scale": 0.2})]
    if name == "sdblock_reset":
        return [S(R, iv=dict(blk(0, 20), k_rel_scale=0.2))]
    if name == "clr3":
        return [S(R, iv={"clr_scale": 3.0})]
    if name == "clr3_reset":
        return [S(R, iv=dict(blk(0, 20), clr_scale=3.0))]
    if name in ("staged_fwd", "staged_rev"):
        allm = np.ones(M.N, bool)
        g1, g2, g3 = POP["G_sens"], POP["G_cent"], POP["G_out"]
        if name == "staged_fwd":      # release sensory -> central -> descending/motor
            pts = [(0, allm), (20.0, ~g1), (30.0, g3 | (~g1 & ~g2 & ~g3)), (40.0, None)]
        else:
            pts = [(0, allm), (20.0, ~g3), (30.0, g1 | (~g1 & ~g2 & ~g3)), (40.0, None)]
        pts[2] = (30.0, pts[2][1] & ~(g3 if name == "staged_fwd" else g1)) if False else pts[2]
        return [S(R, iv={"block_mask": piece(pts)})]
    if name in ("noise", "noise_reset"):
        rng = np.random.default_rng(123); idx = np.sort(rng.choice(M.N, int(0.02 * M.N), replace=False))
        st = [dict(idx=idx, rate=2.0)]
        if name == "noise":
            return [S(20.0, stim=st), S(R - 20.0)]
        return [S(20.0, iv=blk(0, 20)), S(20.0, stim=st), S(R - 40.0)]
    if name == "revive_reset":
        return [S(R, iv=blk(0, 20), revive=True)]
    if name == "adapt":
        return [S(R, mx={"k_adapt": 0.2})]
    if name == "adapt_reset":
        return [S(R, mx={"k_adapt": 0.2}, iv=blk(0, 20))]
    raise ValueError(name)

# ---------------------------------------------------------------- one condition

def evaluate(st, seed, mx, iv_eval, model=None):
    m = model or M
    mxe = dict(mx); mxe["speedup"] = 1.0
    r = fs.run(m, SPONT_MS, metab=mxe, supply=1.0, temperature=22.0, state=st, groups=None, seed=seed + 1000,
               interventions=iv_eval, rec_every_ms=1e9)
    spont = r["counts"][0].copy(); st = r["state"]
    alive = st.alive.copy()
    pc = []
    for j, (name, (idx, rate)) in enumerate(POP["probes"].items()):
        s2 = copy.deepcopy(st)
        rp = fs.run(m, PROBE_MS, stim=[dict(idx=idx, rate=rate)], metab=mxe, supply=1.0, temperature=22.0, state=s2,
                    groups=None, seed=seed * 100 + 5000 + j, interventions=iv_eval, rec_every_ms=1e9)
        pc.append(rp["counts"][0].copy())
    pc = np.array(pc)
    spont_rate = spont / (SPONT_MS * 1e-3); probe_rate = pc / (PROBE_MS * 1e-3)
    ro = POP["readouts"]
    Rm = probe_rate[:, ro] - spont_rate[ro][None, :]
    return dict(spont=spont, pc=pc, spont_rate=spont_rate, probe_rate=probe_rate, R=Rm.astype(np.float32), alive=alive)


def fii_matrix(Rb, Rc):
    U = (Rb > 0) | (Rc > 0)
    if U.sum() < 3 or not (Rc[U] != 0).any():
        return 0.0
    x, y = Rb[U].ravel(), Rc[U].ravel()
    if x.std() == 0 or y.std() == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def make_aged_model(age, mode, seed):
    m = copy.copy(M)
    p = 0.3 * age
    rng = np.random.default_rng(10_000 + seed)
    ne = len(M.indices)
    if mode == "random":
        pe = np.full(ne, p)
    else:   # demand-proportional loss of output synapses (high-demand presynaptic neurons lose more)
        b = np.repeat(POP["base_demand"], np.diff(M.indptr))
        pe = np.clip(p * b / b.mean(), 0, 0.95)
    keep = rng.random(ne) >= pe
    m.weight = (M.weight * keep).astype(np.float32)
    return m, float(1 - keep.mean())


def run_condition(c):
    t0 = time.time()
    seed = c["seed"]
    spc = float(c.get("sp", SP))
    mx = dict(CAL); mx.update(c.get("mx", {}))
    T_ins = c.get("T_ins", 22.0); T_rep = c.get("T_rep", T_ins); T_pre = c.get("T_pre", T_ins)
    model = None; extra = {}
    iv_all = dict(c.get("iv_all", {}))          # interventions applied in every phase (ageing clearance etc.)
    if c.get("age") is not None:
        model, lost = make_aged_model(c["age"], c["age_mode"], seed); extra["syn_lost_frac"] = lost
        mx["reserve"] = 4.0 * (1 - 0.5 * c["age"]) if not c.get("reserve_boost") else 4.0
        iv_all.setdefault("clr_scale", (1 - 0.5 * c["age"]) * c.get("clr_boost", 1.0))
    prot = {}
    if c.get("protect"):
        kind, k = c["protect"]
        n = M.N if kind == "whole" else int(round(k / 100.0 * M.N))
        if kind == "hub":
            sel = np.argsort(-DEG, kind="stable")[:n]
        elif kind == "rand":
            sel = np.random.default_rng(777 + seed).choice(M.N, n, replace=False)
        else:
            sel = np.arange(M.N)
        sc = np.ones(M.N, np.float32); sc[sel] = 0.0
        prot = {"injury_scale": sc}; extra["n_protected"] = int(n)
    tr = Trace(); st = None
    # --- baseline-only conditions
    if c.get("insult") == "none":
        st = phase(None, PRE_S, spc, mx, temp=22.0, iv=iv_all, seed=seed, model=model)["state"] if c.get("age") is not None else fs.init_state(M, dict(fs.METAB_DEFAULTS, **mx), seed)
        ev = evaluate(st, seed, mx, iv_all or None, model)
        return finish(c, ev, tr, {}, t0, extra, st)
    # --- pre phase (+ preconditioning)
    r = phase(None, PRE_S, spc, mx, temp=T_pre, iv=iv_all, seed=seed, model=model); tr.add(r, PRE_S, spc); st = r["state"]
    if c.get("precond"):
        for dur, sup, iv in [(15.0, 0.0, {}), (20.0, 1.0, {"block_all": True}), (60.0, 1.0, {})]:
            r = phase(st, dur, spc, mx, supply=sup, temp=T_pre, iv=dict(iv_all, **iv), seed=seed, model=model)
            tr.add(r, dur, spc); st = r["state"]
    # --- insult
    ins = c["insult"]; D = float(c["D"])
    if ins == "anoxia":
        S_ins = 0.0
    elif ins.startswith("focal_"):
        S_ins = np.where(POP["core_" + ins[6:]], 0.0, 1.0)
    elif ins.startswith("partial_"):
        S_ins = float(ins.split("_")[1])
    elif ins == "heat":
        S_ins = 1.0
    else:
        raise ValueError(ins)
    iv_ins = dict(iv_all); iv_ins.update(prot); iv_ins.update(c.get("iv_ins", {}))
    mark = {"insult_start": tr.t0}
    segs = [(D, spc)] if not c.get("coma_split") or D <= 60 else [(60.0, spc), (D - 60.0, SP_COMA)]
    for dur, s_ in segs:
        r = phase(st, dur, s_, mx, supply=S_ins, temp=T_ins, iv=iv_ins, seed=seed, model=model); tr.add(r, dur, s_); st = r["state"]
    mark["reperf_start"] = tr.t0
    dead_after_insult = float(1 - st.alive.mean())
    # --- reperfusion
    for sub in reperf_protocol(c["reperf"], T_rep):
        iv = dict(iv_all); iv.update(prot); iv.update(sub["iv"])
        if sub["revive"]:
            iv["revive"] = True
        mxs = dict(mx); mxs.update(sub["mx"])
        stim = None
        if sub["stim"]:
            stim = [dict(s_, t_on=0.0, t_off=1e18) for s_ in sub["stim"]]
        r = phase(st, sub["dur"], spc, mxs, supply=sub["supply"], temp=sub["temp"], iv=iv, stim=stim, seed=seed,
                  inh=sub["inh"], model=model)
        tr.add(r, sub["dur"], spc); st = r["state"]
    mark["obs_start"] = tr.t0
    r = phase(st, OBS_S, spc, mx, supply=1.0, temp=T_rep, iv=iv_all, seed=seed, model=model); tr.add(r, OBS_S, spc); st = r["state"]
    mark["end"] = tr.t0
    ev = evaluate(st, seed, mx, iv_all or None, model)
    extra["dead_after_insult"] = dead_after_insult
    return finish(c, ev, tr, mark, t0, extra, st)


def finish(c, ev, tr, mark, t0, extra, st):
    A = tr.arrays(); alive = ev["alive"]
    out = dict(c); out.update(extra)
    out["dead_frac"] = float(1 - alive.mean())
    rch, rper = reach(alive); out["reach"] = rch; out["reach_per_probe"] = rper
    sr = ev["spont_rate"]
    out["FSA"] = float((ev["spont"] > 0).mean()); out["FSA_alive"] = float((ev["spont"][alive] > 0).mean()) if alive.any() else 0.0
    out["spont_Hz_per_neuron"] = float(sr.mean())
    evoked = (ev["probe_rate"] >= sr[None, :] + 8.0)
    out["n_evoked"] = int(evoked.any(0).sum())
    out["n_evoked_per_probe"] = evoked.sum(1).astype(int).tolist()
    R = ev["R"]; ro = POP["readouts"]
    out["readout_resp_max"] = float(R.max()); out["n_readout_resp"] = int((R > 8.0).any(0).sum())
    pos = {int(i): k for k, i in enumerate(ro)}
    nm = {}
    for name, idx in POP["named"].items():
        nm[name] = [round(float(np.mean([R[j, pos[int(i)]] for i in idx if int(i) in pos])) if any(int(i) in pos for i in idx) else float("nan"), 2) for j in range(R.shape[0])]
    out["named_resp"] = nm; out["probe_names"] = list(POP["probes"].keys())
    if BASE is not None:
        out["FII"] = fii_matrix(BASE["R"], R)
        nf = []
        for j in range(R.shape[0]):
            d = fs.functional_identity_index(BASE["rates"][j], ev["probe_rate"][j]); nf.append(d["FII"] if d["FII"] == d["FII"] else 0.0)
        out["FII_neuron"] = float(np.mean(nf)); out["FII_neuron_per_probe"] = [round(x, 4) for x in nf]
        out["resp_mag_ratio"] = float(np.abs(R).sum() / max(np.abs(BASE["R"]).sum(), 1e-9))
    if A["t"].size and mark:
        t = A["t"]; ins0 = mark["insult_start"]; rp0 = mark["reperf_start"]
        w = (t >= ins0) & (t < rp0)
        dep = A["dep"]
        on = np.flatnonzero(w & (dep > 0.5))
        out["SD_onset_s"] = float(t[on[0]] - ins0) if on.size else None
        on1 = np.flatnonzero(w & (dep > 0.1))
        out["SD10_onset_s"] = float(t[on1[0]] - ins0) if on1.size else None
        out["K_peak"] = float(A["Kx"].max()); out["Kmean_peak"] = float(A["Km"].max())
        out["e_min"] = float(A["e"].min()); out["dep_max"] = float(dep.max())
        wr = t >= rp0
        rep = np.flatnonzero(wr & (dep < 0.01))
        out["repol_s_after_reperf"] = float(t[rep[0]] - rp0) if rep.size else None
        out["K_mean_end"] = float(A["Km"][-1]); out["e_end"] = float(A["e"][-1])
        pt, pp = A["pop_t"], A["pop"]
        wr2 = (pt >= rp0) & (pt < mark["obs_start"])
        out["surge_peak_Hz"] = float(pp[wr2].max()) if wr2.any() else None
        out["reperf_mean_Hz"] = float(pp[wr2].mean()) if wr2.any() else None
        wo = pt >= mark["end"] - 10.0
        out["obs_late_Hz"] = float(pp[wo].mean()) if wo.any() else None
        out["marks"] = mark
    out["wall_s"] = round(time.time() - t0, 1)
    try:
        import resource; out["maxrss_MB"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
    except Exception:
        pass
    # sparse spike counts per eval epoch (spont, probes)
    cnt = np.vstack([ev["spont"][None, :], ev["pc"]])
    sp_i = [np.flatnonzero(x).astype(np.int32) for x in cnt]
    sp_c = [np.minimum(x[i], 65535).astype(np.uint16) for x, i in zip(cnt, sp_i)]
    heavy = dict(R=R, alive=np.packbits(alive), probe_rate=(ev["probe_rate"] if c["arm"] == "base" else np.zeros(0)),
                 **{f"trace_{k}": v for k, v in A.items()},
                 **{f"cnt_i{j}": a for j, a in enumerate(sp_i)}, **{f"cnt_c{j}": a for j, a in enumerate(sp_c)})
    light = dict(R=R.astype(np.float16), **{f"trace_{k}": v for k, v in A.items()})
    if c["arm"] == "base":
        light["probe_rate"] = ev["probe_rate"].astype(np.float32)
    return out, heavy, light


def jdump(o, path):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(o, f, default=lambda x: (x.tolist() if hasattr(x, "tolist") else str(x)))
    os.replace(tmp, path)


def worker(c):
    arm, cid = c["arm"], c["id"]
    dj = os.path.join(RES, arm, cid + ".json")
    if os.path.exists(dj):
        return cid, "skip", 0
    try:
        out, heavy, light = run_condition(c)
        os.makedirs(os.path.join(RES, arm), exist_ok=True); os.makedirs(os.path.join(OUT, "results", arm), exist_ok=True)
        np.savez_compressed(os.path.join(RES, arm, cid + ".npz"), **heavy)
        np.savez_compressed(os.path.join(OUT, "results", arm, cid + ".npz"), **light)
        jdump(out, os.path.join(OUT, "results", arm, cid + ".json"))
        jdump(out, dj)                       # written last = completion marker
        return cid, "ok", out["wall_s"]
    except Exception:
        os.makedirs(os.path.join(OUT, "errors"), exist_ok=True)
        tb = traceback.format_exc()
        open(os.path.join(OUT, "errors", cid + ".txt"), "w").write(tb)
        try:
            os.makedirs(os.path.join(RES, "_errors"), exist_ok=True)
            open(os.path.join(RES, "_errors", cid + ".txt"), "w").write(tb)
        except Exception:
            pass
        return cid, "error", 0

# ---------------------------------------------------------------- grid

def grid():
    G = []
    def add(arm, name, seeds=None, **kw):
        if seeds is None:
            seeds = (0, 1, 2, 3, 4) if arm in ("B", "C") else (0, 1, 2)
        for s in seeds:
            G.append(dict(arm=arm, id=f"{arm}_{name}_s{s}", name=name, seed=s, **kw))
    for s in range(10):
        G.append(dict(arm="base", id=f"base_s{s}", name="base", seed=s, insult="none"))
    # A. calibration
    for tau in (20.0, 30.0, 40.0):
        for spd in (25.0, 50.0, 100.0):
            add("A", f"tau{int(tau)}_sp{int(spd)}_D60_reset", insult="anoxia", D=60, reperf="reset",
                mx=dict(tau_atp_s=tau), sp=spd)
    for spd in (25.0, 100.0):
        add("A", f"sp{int(spd)}_D30_naive", insult="anoxia", D=30, reperf="naive", sp=spd)
    add("A", "sham_naive", insult="anoxia", D=0.0, reperf="naive")
    add("A", "sham_reset", insult="anoxia", D=0.0, reperf="reset")
    for j50 in (60.0, 300.0):
        add("A", f"J50_{int(j50)}_D120_reset", insult="anoxia", D=120, reperf="reset", mx=dict(J50_s=j50))
    for arr in (1.0, 0.05):
        for j50 in (600.0, 3600.0, 14400.0):
            for D in (960, 3600):
                add("A", f"arrest{arr}_J50_{int(j50)}_D{D}_reset", insult="anoxia", D=D, reperf="reset",
                    mx=dict(coma_arrest=arr, J50_s=j50), coma_split=True)
    # B. brain-death map
    for D in (7, 30, 60, 120, 240, 480, 960):
        for T in (17.0, 22.0, 29.0, 33.0):
            for rp in ("naive", "reset"):
                add("B", f"anox_D{D}_T{int(T)}_{rp}", insult="anoxia", D=D, T_ins=T, reperf=rp)
    for q in (2.0, 3.0):
        for T in (17.0, 29.0):
            for D in (120, 240, 480, 960):
                add("B", f"q10_{q}_anox_D{D}_T{int(T)}_reset", insult="anoxia", D=D, T_ins=T, reperf="reset",
                    mx=dict(q10_demand=q, q10_injury=q))
    for core in ("OL_R", "hemi_R"):
        for D in (60, 240, 960):
            for rp in ("naive", "reset"):
                add("B", f"focal_{core}_D{D}_{rp}", insult=f"focal_{core}", D=D, reperf=rp)
    for S in (0.1, 0.25):
        for D in (240, 960):
            add("B", f"partial{S}_D{D}_reset", insult=f"partial_{S}", D=D, reperf="reset")
    for T in (33.0, 37.0, 40.0):
        for rp in ("naive", "reset"):
            add("B", f"heat_T{int(T)}_D300_{rp}", insult="heat", D=300, T_ins=T, T_rep=22.0, reperf=rp)
    # C. resuscitation protocols
    PROT = (["naive"] + [f"block{d}" for d in (2, 5, 10, 20, 40)] + ["block20_del5", "block20_del20", "voff-10", "voff-20",
            "gaba2", "gaba4", "pulsed3x10", "pulsed6x5", "ramp10", "ramp30", "ramp60", "cool17", "cool12", "sdblock",
            "sdblock_reset", "clr3", "clr3_reset", "staged_fwd", "staged_rev", "noise", "noise_reset", "revive_reset",
            "adapt", "adapt_reset", "naive_noinj", "ramp30_noinj"])
    for D in (30, 180, 360):
        for p in PROT:
            add("C", f"D{D}_{p}", insult="anoxia", D=D, reperf=p)
    # D. tolerance / preconditioning
    for D in (180, 360):
        for ds in (0.2, 0.5, 0.8):
            add("D", f"hypomet{ds}_D{D}_reset", insult="anoxia", D=D, reperf="reset", iv_ins=dict(demand_scale=ds), pre_ds=ds)
        for Tc in (12.0, 17.0):
            add("D", f"precool{int(Tc)}_D{D}_reset", insult="anoxia", D=D, T_ins=Tc, T_pre=Tc, T_rep=22.0, reperf="reset")
        add("D", f"precond15_D{D}_reset", insult="anoxia", D=D, reperf="reset", precond=True)
        add("D", f"arrest0.05_D{D}_reset", insult="anoxia", D=D, reperf="reset", mx=dict(coma_arrest=0.05))
    # E. ageing
    for a in (0.25, 0.5, 0.75, 1.0):
        for mode in ("demand", "random"):
            add("E", f"age{a}_{mode}_none", insult="none", age=a, age_mode=mode)
            add("E", f"age{a}_{mode}_D60_reset", insult="anoxia", D=60, reperf="reset", age=a, age_mode=mode)
    for a in (0.0, 0.25, 0.5, 0.75, 1.0):
        add("E", f"age{a}_demand_partial0.3_D120_reset", insult="partial_0.3", D=120, reperf="reset", age=a, age_mode="demand")
    add("E", "age1.0_demand_D60_reset_reserveboost", insult="anoxia", D=60, reperf="reset", age=1.0, age_mode="demand", reserve_boost=True)
    add("E", "age1.0_demand_D60_reset_clrboost", insult="anoxia", D=60, reperf="reset", age=1.0, age_mode="demand", clr_boost=2.0)
    add("E", "age1.0_demand_D60_reset_hypomet0.8", insult="anoxia", D=60, reperf="reset", age=1.0, age_mode="demand", iv_ins=dict(demand_scale=0.8))
    # H. targeted (hub) neuroprotection (lead request)
    for rp in ("naive", "reset"):
        for k in (1, 2, 5, 10):
            add("H", f"hub{k}_D360_{rp}", insult="anoxia", D=360, reperf=rp, protect=("hub", k))
            add("H", f"rand{k}_D360_{rp}", insult="anoxia", D=360, reperf=rp, protect=("rand", k))
        add("H", f"whole_D360_{rp}", insult="anoxia", D=360, reperf=rp, protect=("whole", 100))
    return G


def _fix(c):
    # arm-A speedup variants: protocol compression set per condition
    return c


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--arms", default="")
    ap.add_argument("--ids", default=""); ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    G = grid()
    if a.list:
        from collections import Counter
        print(len(G), Counter(c["arm"] for c in G)); sys.exit()
    if a.arms:
        G = [c for c in G if c["arm"] in a.arms.split(",")]
    if a.ids:
        G = [c for c in G if c["id"] in a.ids.split(",")]
    G = G[a.shard::a.nshard]
    os.makedirs(OUT, exist_ok=True)
    todo = [c for c in G if not os.path.exists(os.path.join(RES, c["arm"], c["id"] + ".json"))]
    print(f"shard {a.shard}/{a.nshard}: {len(G)} conditions, {len(todo)} to run, workers {a.workers}", flush=True)
    t0 = time.time(); setup(); load_base()
    print(f"model loaded {time.time()-t0:.1f}s; baseline for inline FII: {None if BASE is None else BASE['n']}", flush=True)
    log = open(os.path.join(OUT, f"progress_shard{a.shard}.log"), "a")
    if a.workers <= 1:
        it = map(worker, todo)
    else:
        import multiprocessing as mp_
        pool = mp_.get_context("fork").Pool(a.workers, maxtasksperchild=4)
        it = pool.imap_unordered(worker, todo, chunksize=1)
    n = 0
    for cid, stt, w in it:
        n += 1
        msg = f"[{time.time()-t0:8.0f}s] {n}/{len(todo)} {cid} {stt} {w}s"
        print(msg, flush=True); log.write(msg + "\n"); log.flush()
    print("DONE", time.time() - t0, flush=True)
