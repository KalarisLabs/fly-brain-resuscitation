"""Mechanism & robustness follow-up driver (FAFB v783 princeton5; flysim_mech = flysim_dyn + optional excitotoxic term).

python mech_sweep.py --shard A|B --workers W [--arms L,X,N,T,M,Y] [--smoke] [--list]

Arms
  L  Latin-hypercube robustness sweep: 1 calibrated + 48 LHS parameter sets x {D30,D180} x {naive,reset} x 3 seeds,
     plus a per-set intact baseline (insult none, seeds 100-103) for set-specific FII.
  X  Excitotoxic reperfusion variant (k_exc): H3 re-test, 12 reperfusion protocols x D{30,180,360} x k{1,3} x 5 seeds,
     reset-combination protocols at k=0, calibration k{0.3,10} at D180.
  N  Attractor anatomy: naive reperfusion D{30,180} x seeds 0-4, neuron-level epoch spike counts, first-spike latency,
     end-state effective resting potential.  -> driver sets (built on Modal, shard A) -> T.
  T  Targeted 20-s reset (block_mask) of driver sets (sustain / ignite), top-degree hubs, random, k in
     {0.1,0.5,1,2,5,10,20} % vs whole-brain block and naive, D{30,180}, HELD-OUT seeds 10-14.
  M  Mechanism dissection: naive for 20 s then a manipulation for 40 s (synapses off, excitation off, inhibition off,
     K+ clamp 5 mM, k_spk 0, k_rel 0), then 20 s free observation; D{30,180} x 5 seeds.
  Y  Bistability / hysteresis in the intact network (no insult): K+ clamp kicks x v_offset, x excitatory gain,
     spike kicks, quench by brief total block, kicks after a D30 + reset recovery, slow v_offset ramp up/down.
Resumable: skip if RES/<arm>/<id>.json exists.  One writer per path.
"""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src", _os.path.join("experiments", "dynamics")):
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import flysim_mech as fs
import sweep as SW
SW.fs = fs                                   # all template helpers now use the mech variant

RES = os.environ.get("RES_DIR", os.path.join(_P.OUT, "followup", "mechanism", "raw", "results"))
OUT = os.environ.get("OUT_DIR", "out")
BASE_RES = os.environ.get("BASE_RES", os.path.join(_P.OUT, "dynamics", "raw"))
SP = SW.SP
CAL = dict(SW.CAL)
M2 = None
DRVFILE = None
MV = {}; MV_N = {}

# ------------------------------------------------------------------------------------------------ setup

def setup():
    global M2, DRVFILE
    SW.setup()
    M = SW.M
    deg = np.diff(M.indptr)
    s2 = fs._sign_vector(M.neurons.nt.values, glut_inhibitory=False)
    M2 = copy.copy(M)
    M2.weight = (np.abs(M.weight) * np.repeat(s2, deg).astype(np.float32)).astype(np.float32)
    M2.sign = s2
    t = M.neurons
    cls = t["class"].fillna("").values; pt = t.primary_type.fillna("")
    alln = (cls == "ALLN") | pt.str.contains(r"^(?:l|v|v2|il|il2)LN", regex=True).values
    nocall = t.nt_neuron.isna().values
    low = nocall | (t.nt_type_score.fillna(0).values < 0.5)
    for var, mask, val in (("A", alln, -1), ("B", low, -1), ("C", nocall, 0)):
        sv = M.sign.copy(); sv[mask] = val
        mv = copy.copy(M); mv.sign = sv
        mv.weight = (np.abs(M.weight) * np.repeat(sv, deg).astype(np.float32)).astype(np.float32)
        MV[var] = mv; MV_N[var] = int(mask.sum())
    SW.RES = BASE_RES
    SW.load_base()
    DRVFILE = os.path.join(RES, "drivers.npz")
    if not os.path.exists(DRVFILE) and os.path.exists(os.path.join(HERE, "drivers.npz")):
        import shutil; os.makedirs(RES, exist_ok=True); shutil.copy(os.path.join(HERE, "drivers.npz"), DRVFILE)


def model_for(c):
    if c.get("signvar"):
        base = MV[c["signvar"]]
        if c.get("w_syn") is not None:
            m = copy.copy(base); m.w_syn = float(c["w_syn"]); return m
        return base
    base = M2 if c.get("glut_exc") else SW.M
    if c.get("w_syn") is not None and abs(float(c["w_syn"]) - SW.M.w_syn) > 1e-9:
        m = copy.copy(base); m.w_syn = float(c["w_syn"]); return m
    return None if base is SW.M else base


class Trace2(SW.Trace):
    def __init__(self):
        super().__init__(); self.d["Jx"] = []

    def add(self, r, dur_s, sp_):
        sl = r["slow"]
        if np.asarray(sl["t_s"]).size:
            self.d["Jx"].append(np.asarray(sl["Jx_mean"], np.float32))
        super().add(r, dur_s, sp_)


def phase2(st, dur_s, sp_, mx, supply=1.0, temp=22.0, iv=None, stim=None, seed=0, inh=1.0, exc=1.0, model=None,
           epochs=None, rec=True):
    m = model or SW.M
    ts0 = (st.t_ms * 1e-3 * sp_) if st is not None else 0.0
    mxx = dict(mx); mxx["speedup"] = sp_
    ivw = {k: (SW.wrap(v, ts0) if callable(v) else v) for k, v in (iv or {}).items()}
    return fs.run(m, dur_s * 1000.0 / sp_, stim=stim, metab=mxx, supply=SW.wrap(supply, ts0), temperature=SW.wrap(temp, ts0),
                  interventions=ivw, state=st, groups=None, rec_every_ms=(500.0 / sp_ if rec else 1e9), seed=seed,
                  inh_gain=inh, exc_gain=exc, epochs=epochs)

# ------------------------------------------------------------------------------------------------ protocols

def blk(on, off):
    return {"block_all": SW.piece([(0, False), (on, True), (on + off, False)])}


def target_set(kind, k, seed, var=None):
    N = SW.M.N
    n = max(1, int(round(k / 100.0 * N)))
    if kind == "hub":
        idx = np.argsort(-SW.DEG, kind="stable")[:n]
    elif kind == "rand":
        idx = np.random.default_rng(4242 + 97 * seed + int(round(k * 100))).choice(N, n, replace=False)
    elif kind.endswith("180"):
        z = np.load(os.path.join(HERE, "drivers_extra.npz"))      # D180-only rankings (computed locally, same rule)
        idx = z["rank_" + kind][:n]
    else:
        z = np.load(DRVFILE if not var else os.path.join(RES, f"drivers_{var}.npz"))
        idx = z["rank_" + kind][:n]
    return np.sort(np.asarray(idx, np.int64))


def protocol(name, seed=0, var=None):
    S = lambda dur, **kw: dict(dict(dur=dur, supply=1.0, iv={}, inh=1.0, exc=1.0, stim=None, mx={}), **kw)
    R = SW.REPERF_S
    ramp = lambda Tr: fs.ramp([(0, 0.0), (Tr, 1.0), (1e6, 1.0)])
    pulsed = SW.piece([(0, 1.0), (10, 0.0), (20, 1.0), (30, 0.0), (40, 1.0)])
    if name == "naive":
        return [S(R)]
    if name == "reset":
        return [S(R, iv=blk(0, 20))]
    if name == "block20_del20":
        return [S(R, iv=blk(20, 20))]
    if name in ("block40", "block60"):
        return [S(R, iv=blk(0, float(name[5:])))]
    if name in ("ramp10", "ramp30", "ramp60"):
        return [S(R, supply=ramp(float(name[4:])))]
    if name == "pulsed3x10":
        return [S(R, supply=pulsed)]
    if name in ("ramp10_reset", "ramp30_reset"):
        return [S(R, supply=ramp(float(name[4:6])), iv=blk(0, 20))]
    if name == "ramp30_del20":
        return [S(R, supply=ramp(30.0), iv=blk(20, 20))]
    if name == "ramp60_del40":
        return [S(R, supply=ramp(60.0), iv=blk(40, 20))]
    if name == "pulsed3x10_del40":
        return [S(R, supply=pulsed, iv=blk(40, 20))]
    if name.startswith("tblk:"):
        _, kind, k = name.split(":")
        idx = target_set(kind, float(k), seed, var)
        return [S(R, iv={"block_mask": SW.piece([(0, idx), (20.0, None)])})]
    if name.startswith("dis:"):
        mnp = name[4:]
        kw = {}
        if mnp == "nosyn":
            kw = dict(exc=0.0, inh=0.0)
        elif mnp == "noexc":
            kw = dict(exc=0.0)
        elif mnp == "noinh":
            kw = dict(inh=0.0)
        elif mnp == "Kclamp5":
            kw = dict(iv={"K_clamp": 5.0})
        elif mnp == "kspk0":
            kw = dict(mx={"k_spk": 0.0})
        elif mnp == "krel0":
            kw = dict(mx={"k_rel": 0.0})
        else:
            raise ValueError(name)
        return [S(20.0), S(R - 20.0, **kw)]
    raise ValueError(name)

# ------------------------------------------------------------------------------------------------ standard condition

def run_cond(c):
    t0 = time.time(); seed = c["seed"]; spc = SP
    mx = dict(CAL); mx.update(c.get("mx", {}))
    model = model_for(c)
    tr = Trace2(); extra = {}
    if c.get("insult") == "none":
        st = fs.init_state(SW.M, dict(fs.METAB_DEFAULTS, **mx), seed)
        ev = SW.evaluate(st, seed, mx, None, model)
        out, heavy, light = SW.finish(c, ev, tr, {}, t0, extra, st)
        light["probe_rate"] = ev["probe_rate"].astype(np.float32)
        return out, heavy, light
    T = float(c.get("T", 22.0))
    piv = {}
    if c.get("protect"):
        pk = float(c["protect"][1]); npr = int(round(pk / 100.0 * SW.M.N))
        sc_ = np.ones(SW.M.N, np.float32); sc_[np.argsort(-SW.DEG, kind="stable")[:npr]] = 0.0
        piv = {"injury_scale": sc_}; extra["n_protected"] = npr
    r = phase2(None, SW.PRE_S, spc, mx, temp=T, iv=piv, seed=seed, model=model); tr.add(r, SW.PRE_S, spc); st = r["state"]
    D = float(c["D"]); mark = {"insult_start": tr.t0}
    r = phase2(st, D, spc, mx, supply=0.0, temp=T, iv=piv, seed=seed, model=model); tr.add(r, D, spc); st = r["state"]
    mark["reperf_start"] = tr.t0
    extra["dead_after_insult"] = float(1 - st.alive.mean())
    extra["J_mean_after_insult"] = float(st.J.mean())
    win = {}
    rp = c["reperf"]
    if rp.startswith("CL"):
        # closed-loop adaptive reset: block until max compartment K+ < K_th - margin AND mean ATP e > e_crit
        idx = None
        if rp != "CL":
            _, kind, k = rp.split(":"); idx = target_set(kind, float(k), seed, c.get("signvar"))
        mp_ = dict(fs.METAB_DEFAULTS); mp_.update(mx)
        Kth = mp_["K0"] * math.exp((fs.LIF_DEFAULTS["v_th"] - fs.LIF_DEFAULTS["v_rest"]) / (mp_["g_K"] * mp_["RT_F"]))
        Kcrit = Kth - float(c.get("CL_margin", 1.0)); Ecrit = float(c.get("CL_e", 0.9)); Tmax = float(c.get("CL_max", 120.0))
        ivb = {"block_all": True} if idx is None else {"block_mask": idx}
        Tb = 0.0; met = False; em = 0.0
        while Tb < Tmax:
            r = phase2(st, 1.0, spc, mx, iv=dict(piv, **ivb), temp=T, seed=seed, model=model); tr.add(r, 1.0, spc); st = r["state"]; Tb += 1.0
            e_ = st.a / st.a0; em = float(e_[st.alive].mean()) if st.alive.any() else 0.0
            if Tb >= 2 and float(st.K.max()) < Kcrit and em > Ecrit:
                met = True; break
        extra.update(CL_block_s=Tb, CL_met=met, CL_Kcrit=Kcrit, CL_Kmax_release=float(st.K.max()), CL_e_release=em)
        r = phase2(st, 40.0, spc, mx, iv=piv, temp=T, seed=seed, model=model)
        pop = r["pop_ms"].sum(1); win["post_release_mean_Hz"] = float(pop.sum() / (len(pop) * 1e-3) / SW.M.N)
        tr.add(r, 40.0, spc); st = r["state"]
    else:
        for j, sub in enumerate(protocol(rp, seed, c.get("signvar"))):
            mxs = dict(mx); mxs.update(sub["mx"])
            stim = [dict(s_, t_on=0.0, t_off=1e18) for s_ in sub["stim"]] if sub["stim"] else None
            r = phase2(st, sub["dur"], spc, mxs, supply=sub["supply"], iv=dict(piv, **sub["iv"]), temp=T, stim=stim, seed=seed,
                       inh=sub["inh"], exc=sub["exc"], model=model)
            pop = r["pop_ms"].sum(1)
            win[f"sub{j}_mean_Hz"] = float(pop.sum() / (len(pop) * 1e-3) / SW.M.N)
            b = max(1, len(pop) // 4)
            win[f"sub{j}_lastq_Hz"] = float(pop[-b:].sum() / (b * 1e-3) / SW.M.N)
            tr.add(r, sub["dur"], spc); st = r["state"]
    extra.update(win)
    mark["obs_start"] = tr.t0
    r = phase2(st, SW.OBS_S, spc, mx, temp=T, iv=piv, seed=seed, model=model); tr.add(r, SW.OBS_S, spc); st = r["state"]
    mark["end"] = tr.t0
    extra["dead_before_eval"] = float(1 - st.alive.mean())
    Jx = getattr(st, "Jx", None)
    if Jx is not None:
        extra["Jx_mean_end"] = float(Jx.mean())
        dead = ~st.alive
        extra["Jx_mean_dead"] = float(Jx[dead].mean()) if dead.any() else 0.0
        extra["frac_alive_Jx_gt100"] = float((Jx[st.alive] > 100).mean()) if st.alive.any() else 0.0
    extra["J_mean_end"] = float(st.J.mean())
    ev = SW.evaluate(st, seed, mx, None, model)
    out, heavy, light = SW.finish(c, ev, tr, mark, t0, extra, st)
    A = tr.arrays()
    light["trace_Jx"] = A.get("Jx", np.zeros(0, np.float32))
    return out, heavy, light

# ------------------------------------------------------------------------------------------------ anatomy (N)

def vr_eff_of(st, mx, model=None):
    m = model or SW.M
    mp = dict(fs.METAB_DEFAULTS); mp.update(mx)
    e = (st.a / st.a0).astype(np.float64)
    dVa = np.clip((mp["e_ok"] - e) / mp["e_ok"], 0, 1) ** mp["n_atp"] * mp["dV_atp_max"]
    dVk_c = mp["g_K"] * mp["RT_F"] * np.log(np.maximum(st.K, 0.1) / mp["K0"])
    dVk_i = m.Mn @ dVk_c
    vr = fs.LIF_DEFAULTS["v_rest"] + mp["V_span"] * np.tanh((dVa + dVk_i) / mp["V_span"])
    vth = fs.LIF_DEFAULTS["v_th"] + mp["theta_block"] * (1 - st.h)
    return vr.astype(np.float32), vth.astype(np.float32), e.astype(np.float32)


def run_anat(c):
    t0 = time.time(); seed = c["seed"]; spc = SP; mx = dict(CAL)
    mdl = model_for(c)
    tr = Trace2()
    r = phase2(None, SW.PRE_S, spc, mx, seed=seed, model=mdl); tr.add(r, SW.PRE_S, spc); st = r["state"]
    D = float(c["D"])
    r = phase2(st, D, spc, mx, supply=0.0, seed=seed, model=mdl); tr.add(r, D, spc); st = r["state"]
    vr0, vth0, e0 = vr_eff_of(st, mx)
    K0c = st.K.astype(np.float32).copy()
    rep_ms = SW.REPERF_S * 1000.0 / spc
    ep = [(i * 5.0, (i + 1) * 5.0) for i in range(60)] + [(300.0 + i * 50.0, 350.0 + i * 50.0) for i in range(int((rep_ms - 300) / 50))]
    r = phase2(st, SW.REPERF_S, spc, mx, seed=seed, epochs=ep, model=mdl); tr.add(r, SW.REPERF_S, spc); st = r["state"]
    C1 = r["counts"]
    vr1, vth1, e1 = vr_eff_of(st, mx)
    obs_ms = SW.OBS_S * 1000.0 / spc
    ep2 = [(i * 100.0, (i + 1) * 100.0) for i in range(int(obs_ms / 100))]
    r = phase2(st, SW.OBS_S, spc, mx, seed=seed, epochs=ep2, model=mdl); tr.add(r, SW.OBS_S, spc); st = r["state"]
    C2 = r["counts"]
    vr2, vth2, e2 = vr_eff_of(st, mx)
    any_ = C1 > 0
    first_ep = np.where(any_.any(0), any_.argmax(0), -1)
    ep_start = np.array([a for a, b in ep], np.float32)
    first_ms = np.where(first_ep >= 0, ep_start[np.maximum(first_ep, 0)], np.inf).astype(np.float32)
    light = dict(first_ms=first_ms, cnt_rep=C1.sum(0).astype(np.int32), cnt_rep_early=C1[:60].sum(0).astype(np.int32),
                 cnt_obs=C2.sum(0).astype(np.int32), cnt_obs_last=C2[-1].astype(np.int32),
                 vr_rep0=vr0, vth_rep0=vth0, e_rep0=e0, vr_rep_end=vr1, vth_rep_end=vth1, vr_end=vr2, vth_end=vth2, e_end=e2,
                 K_rep0=K0c, K_rep_end=st.K.astype(np.float32), alive=st.alive.copy(), pop5ms=C1[:60].sum(1).astype(np.int32),
                 ep_start=ep_start, **{f"trace_{k}": v for k, v in tr.arrays().items()})
    heavy = dict(C1=np.minimum(C1, 65535).astype(np.uint16), C2=np.minimum(C2, 65535).astype(np.uint16), ep=np.array(ep, np.float32))
    out = dict(c); out["wall_s"] = round(time.time() - t0, 1)
    al = st.alive
    rate_obs = C2.sum(0) / (obs_ms * 1e-3)
    out["dead_frac"] = float(1 - al.mean()); out["obs_mean_Hz"] = float(rate_obs.mean())
    out["frac_active20"] = float((rate_obs >= 20).mean())
    out["frac_supra_end"] = float(((vr2 > vth2) & al).mean())
    out["frac_supra_rep_end"] = float(((vr1 > vth1) & al).mean())
    out["median_first_ms"] = float(np.median(first_ms[np.isfinite(first_ms)])) if np.isfinite(first_ms).any() else None
    return out, heavy, light


def build_drivers(var=None):
    """Driver rankings from the N runs (excitatory neurons only).  sustain: r_i * sum of excitatory synapses onto
    neurons active (>=20 Hz) in the late runaway; ignite: earliest first spike after reperfusion (ties -> sustain)."""
    dfile = DRVFILE if not var else os.path.join(RES, f"drivers_{var}.npz")
    if os.path.exists(dfile):
        return
    M = SW.M if not var else MV[var]; N = M.N
    fl = sorted(glob.glob(os.path.join(RES, "N", "*.anat.npz"))) if not var else sorted(glob.glob(os.path.join(RES, "SN", f"SN_{var}_*.anat.npz")))
    S_ = []; F_ = []
    deg = np.diff(M.indptr); pre = np.repeat(np.arange(N), deg)
    wpos = np.maximum(M.weight, 0).astype(np.float64)
    for f in fl:
        z = np.load(f)
        rate = z["cnt_obs"] / (SW.OBS_S / SP)
        act = (rate >= 20) & z["alive"]
        outw = np.bincount(pre, weights=wpos * act[M.indices], minlength=N)
        S_.append(rate * outw)
        F_.append(np.where(np.isfinite(z["first_ms"]), z["first_ms"], 5000.0))
    s = np.mean(S_, 0); fm = np.mean(F_, 0)
    exc = M.sign > 0
    s_ex = np.where(exc, s, -1.0)
    rank_s = np.argsort(-s_ex, kind="stable")
    rank_i = np.lexsort((-s_ex, np.where(exc, fm, 1e9)))
    tmp = dfile + ".tmp.npz"
    np.savez_compressed(tmp, rank_sustain=rank_s.astype(np.int32), rank_ignite=rank_i.astype(np.int32),
                        s_mean=s.astype(np.float32), first_mean=fm.astype(np.float32), n_runs=len(fl))
    os.replace(tmp, dfile)
    import shutil; shutil.copy(dfile, os.path.join(OUT, os.path.basename(dfile)))

# ------------------------------------------------------------------------------------------------ bistability (Y)

def run_Y(c):
    t0 = time.time(); seed = c["seed"]; spc = SP; mx = dict(CAL)
    mdl = model_for(c)
    _p2 = phase2
    phase2_ = lambda *a, **k: _p2(*a, model=mdl, **k)
    vo = float(c.get("voff", 0.0)); g = float(c.get("gain", 1.0))
    ivb = {"v_offset": vo} if vo else {}
    tr = Trace2(); marks = {}
    noise = None
    if c.get("noise"):
        rng = np.random.default_rng(123); idx = np.sort(rng.choice(SW.M.N, int(0.02 * SW.M.N), replace=False))
        noise = [dict(idx=idx, rate=2.0, t_on=0.0, t_off=1e18)]
    st = None
    if c.get("prereset"):
        r = phase2_(None, SW.PRE_S, spc, mx, seed=seed); tr.add(r, SW.PRE_S, spc); st = r["state"]
        r = phase2_(st, 30.0, spc, mx, supply=0.0, seed=seed); tr.add(r, 30.0, spc); st = r["state"]
        r = phase2_(st, SW.REPERF_S, spc, mx, iv=blk(0, 20), seed=seed); tr.add(r, SW.REPERF_S, spc); st = r["state"]
        r = phase2_(st, SW.OBS_S, spc, mx, seed=seed); tr.add(r, SW.OBS_S, spc); st = r["state"]
    r = phase2_(st, SW.PRE_S, spc, mx, iv=ivb, seed=seed, exc=g, stim=noise); tr.add(r, SW.PRE_S, spc); st = r["state"]
    kind = c["kind"]
    if kind == "hyst":
        Tr = float(c.get("Tramp", 160.0)); lo, hi = float(c.get("vlo", -6.0)), float(c.get("vhi", 10.0))
        f = fs.ramp([(0, lo), (Tr, hi), (2 * Tr, lo), (1e6, lo)])
        marks["ramp_start"] = tr.t0
        r = phase2_(st, 2 * Tr, spc, mx, iv={"v_offset": f}, seed=seed, exc=g, stim=noise); tr.add(r, 2 * Tr, spc); st = r["state"]
    else:
        marks["kick_start"] = tr.t0
        if kind == "K":
            Kk = float(c["Kk"])
            ivk = dict(ivb, **({"K_clamp": Kk} if Kk > 0 else {}))
            r = phase2_(st, 2.0, spc, mx, iv=ivk, seed=seed, exc=g, stim=noise)
        elif kind == "spk":
            rng = np.random.default_rng(555 + seed)
            idx = np.sort(rng.choice(SW.M.N, int(c["frac"] * SW.M.N), replace=False))
            r = phase2_(st, 2.0, spc, mx, iv=ivb, seed=seed, exc=g, stim=[dict(idx=idx, rate=100.0, t_on=0.0, t_off=1e18)])
        else:
            raise ValueError(kind)
        tr.add(r, 2.0, spc); st = r["state"]
        marks["free_start"] = tr.t0
        if c.get("quench") is not None:
            r = phase2_(st, 20.0, spc, mx, iv=ivb, seed=seed, exc=g); tr.add(r, 20.0, spc); st = r["state"]
            marks["pre_quench_end"] = tr.t0
            B = float(c["quench"])
            r = phase2_(st, B, spc, mx, iv=dict(ivb, block_all=True), seed=seed, exc=g); tr.add(r, B, spc); st = r["state"]
            marks["quench_end"] = tr.t0
            r = phase2_(st, 30.0, spc, mx, iv=ivb, seed=seed, exc=g); tr.add(r, 30.0, spc); st = r["state"]
        else:
            r = phase2_(st, float(c.get("free_s", 40.0)), spc, mx, iv=ivb, seed=seed, exc=g, stim=noise)
            tr.add(r, float(c.get("free_s", 40.0)), spc); st = r["state"]
    marks["end"] = tr.t0
    A = tr.arrays()
    out = dict(c); out["marks"] = marks; out["wall_s"] = round(time.time() - t0, 1)
    pt, pp = A["pop_t"], A["pop"]
    late = pt >= marks["end"] - 10.0
    out["late_Hz"] = float(pp[late].mean()) if late.any() else None
    out["peak_Hz"] = float(pp.max()) if pp.size else None
    out["K_mean_end"] = float(A["Km"][-1]) if A["Km"].size else None
    out["K_max_end"] = float(A["Kx"][-1]) if A["Kx"].size else None
    out["dead_frac"] = float(1 - st.alive.mean())
    if "pre_quench_end" in marks:
        w = (pt >= marks["pre_quench_end"] - 5.0) & (pt < marks["pre_quench_end"])
        out["prequench_Hz"] = float(pp[w].mean()) if w.any() else None
    vr, vth, e = vr_eff_of(st, mx)
    out["frac_supra_end"] = float(((vr + vo > vth) & st.alive).mean())
    out["runaway"] = bool(out["late_Hz"] is not None and out["late_Hz"] > 0.5)
    light = {f"trace_{k}": v for k, v in A.items()}
    return out, light, light

# ------------------------------------------------------------------------------------------------ grid

LHS_KEYS = ["w_syn", "tau_atp_s", "J50_s", "g_K", "k_rel", "k_clr", "reserve", "demand_beta", "K_max", "glut_inh"]
LHS_RANGE = dict(w_syn=(0.15, 0.25, "lin"), tau_atp_s=(15.0, 60.0, "log"), J50_s=(300.0, 1200.0, "log"),
                 g_K=(0.3, 0.9, "lin"), k_rel=(20.0, 120.0, "log"), k_clr=(0.5, 2.0, "log"), reserve=(2.0, 8.0, "log"),
                 demand_beta=(0.25, 1.0, "lin"), K_max=(50.0, 100.0, "lin"), glut_inh=(0, 1, "bin"))
DEFAULT_SET = dict(w_syn=0.185, tau_atp_s=30.0, J50_s=600.0, g_K=0.6, k_rel=60.0, k_clr=1.0, reserve=4.0,
                   demand_beta=0.5, K_max=80.0, glut_inh=1)


def lhs_sets(n=48, seed=2026):
    from scipy.stats import qmc
    U = qmc.LatinHypercube(d=len(LHS_KEYS), seed=seed).random(n)
    sets = [dict(DEFAULT_SET, set=0)]
    for i in range(n):
        p = {}
        for j, k in enumerate(LHS_KEYS):
            lo, hi, kind = LHS_RANGE[k]; u = U[i, j]
            if kind == "lin":
                p[k] = lo + u * (hi - lo)
            elif kind == "log":
                p[k] = float(np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo))))
            else:
                p[k] = int(u < 0.5)
            if k != "glut_inh":
                p[k] = round(float(p[k]), 4)
        p["set"] = i + 1
        sets.append(p)
    return sets


def lhs_cond(p):
    mx = {k: p[k] for k in ["tau_atp_s", "J50_s", "g_K", "k_rel", "k_clr", "reserve", "demand_beta", "K_max"]}
    return dict(mx=mx, w_syn=p["w_syn"], glut_exc=(p["glut_inh"] == 0), lhs=p)


X_PROT = ["naive", "ramp10", "ramp30", "ramp60", "pulsed3x10", "reset", "block20_del20", "ramp10_reset", "ramp30_reset",
          "ramp30_del20", "ramp60_del40", "pulsed3x10_del40"]
T_KINDS = ["sustain", "ignite", "hub", "rand"]
T_KS = [0.02, 0.05, 0.1, 0.5, 1, 2, 5, 10, 20]
M_PROT = ["naive", "dis:nosyn", "dis:noexc", "dis:noinh", "dis:Kclamp5", "dis:kspk0", "dis:krel0"]


def grid():
    G = []
    def add(arm, name, seeds, **kw):
        for s in seeds:
            G.append(dict(arm=arm, id=f"{arm}_{name}_s{s}", name=name, seed=s, **kw))
    # L
    for p in lhs_sets():
        lc = lhs_cond(p)
        add("L", f"set{p['set']:02d}_base", (100, 101, 102, 103), insult="none", **lc)
        for D in (30, 180):
            for rp in ("naive", "reset"):
                add("L", f"set{p['set']:02d}_D{D}_{rp}", (0, 1, 2), insult="anoxia", D=D, reperf=rp, **lc)
    # X
    S5 = (0, 1, 2, 3, 4)
    for kx in (1.0, 3.0):
        for D in (30, 180, 360):
            for p in X_PROT:
                add("X", f"k{kx:g}_D{D}_{p}", S5, insult="anoxia", D=D, reperf=p, mx=dict(k_exc=kx), k_exc=kx)
    for D in (30, 180, 360):
        for p in ["ramp10_reset", "ramp30_reset", "ramp30_del20", "ramp60_del40", "pulsed3x10_del40"]:
            add("X", f"k0_D{D}_{p}", S5, insult="anoxia", D=D, reperf=p, k_exc=0.0)
    for D in (30, 180, 360):                     # follow-up: strong excitotoxicity H3 (k=10)
        for p in X_PROT:
            add("X", f"k10_D{D}_{p}", S5, insult="anoxia", D=D, reperf=p, mx=dict(k_exc=10.0), k_exc=10.0)
    for kx in (0.3, 10.0):
        for p in ("naive", "reset"):
            add("X", f"k{kx:g}_D180_{p}", S5, insult="anoxia", D=180, reperf=p, mx=dict(k_exc=kx), k_exc=kx)
    # N
    for D in (30, 180):
        add("N", f"anat_D{D}", S5, insult="anoxia", D=D, reperf="naive")
    # T (held-out seeds)
    S10 = (10, 11, 12, 13, 14)
    for D in (30, 180):
        add("T", f"D{D}_naive", S10, insult="anoxia", D=D, reperf="naive")
        add("T", f"D{D}_reset", S10, insult="anoxia", D=D, reperf="reset")
        for kind in T_KINDS:
            for k in T_KS:
                add("T", f"D{D}_tblk_{kind}_{k:g}", S10, insult="anoxia", D=D, reperf=f"tblk:{kind}:{k:g}", tkind=kind, tk=k)
    for kind, ks, D in [("sustain", (12.5, 15, 17.5), 180), ("ignite", (12.5, 15, 17.5), 180), ("hub", (30, 50), 180),
                        ("rand", (30, 50), 180), ("sustain180", (2, 5, 10, 15, 20), 180), ("ignite180", (2, 5, 10, 15, 20), 180),
                        ("sustain", (0.2, 0.3), 30), ("ignite", (0.2, 0.3), 30), ("hub", (1.5,), 30)]:
        for k in ks:                             # follow-up: bracket the minimal rescuing set
            add("T", f"D{D}_tblk_{kind}_{k:g}", S10, insult="anoxia", D=D, reperf=f"tblk:{kind}:{k:g}", tkind=kind, tk=k)
    # F (follow-up): longer total block for every LHS set at D180
    for p in lhs_sets():
        lc = lhs_cond(p)
        for rp in ("block40", "block60"):
            add("F", f"set{p['set']:02d}_D180_{rp}", (0, 1, 2), insult="anoxia", D=180, reperf=rp, **lc)
    # M
    for D in (30, 180):
        for p in M_PROT:
            add("M", f"D{D}_{p.replace(':', '_')}", S5, insult="anoxia", D=D, reperf=p)
    # Y
    for vo in (-6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0):
        for Kk in (0.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0, 30.0):
            add("Y", f"kick_v{vo:g}_K{Kk:g}", S5 if vo == 0.0 else (0, 1, 2), kind="K", voff=vo, Kk=Kk)
    for g in (0.5, 0.75, 1.25, 1.5):
        for Kk in (0.0, 8.0, 10.0, 12.0, 15.0, 20.0):
            add("Y", f"kick_g{g:g}_K{Kk:g}", (0, 1, 2), kind="K", gain=g, Kk=Kk)
    for fr in (0.01, 0.1, 0.5):
        add("Y", f"spk_f{fr:g}", S5, kind="spk", frac=fr)
    for vo in (0.0, 4.0):
        for B in (1.0, 2.0, 5.0, 10.0, 20.0):
            add("Y", f"quench_v{vo:g}_B{B:g}", (0, 1, 2), kind="K", voff=vo, Kk=15.0, quench=B)
    for Kk in (0.0, 10.0, 15.0):
        add("Y", f"postreset_K{Kk:g}", (0, 1, 2), kind="K", Kk=Kk, prereset=True)
    # S (sign-robustness): variants A (AL local neurons inhibitory), B (no-call / low-confidence NT inhibitory),
    # C (no-call NT sign 0)
    for var in ("A", "B", "C"):
        add("S", f"{var}_base", (100, 101, 102, 103, 104), insult="none", signvar=var)
        for D in (30, 180):
            for rp in ("naive", "reset"):
                add("S", f"{var}_D{D}_{rp}", S5, insult="anoxia", D=D, reperf=rp, signvar=var)
    for var in ("A", "B"):
        for D in (30, 180):
            add("SN", f"{var}_anat_D{D}", S5, insult="anoxia", D=D, reperf="naive", signvar=var)
            add("ST", f"{var}_D{D}_naive", S10, insult="anoxia", D=D, reperf="naive", signvar=var)
            add("ST", f"{var}_D{D}_reset", S10, insult="anoxia", D=D, reperf="reset", signvar=var)
            for kind, ks in (("sustain", T_KS[2:]), ("hub", T_KS[2:]), ("rand", T_KS[2:]), ("ignite", (0.5, 2, 20))):
                for k in ks:
                    add("ST", f"{var}_D{D}_tblk_{kind}_{k:g}", S10, insult="anoxia", D=D, reperf=f"tblk:{kind}:{k:g}", tkind=kind, tk=k, signvar=var)
        for Kk in (0.0, 9.0, 10.0, 12.0, 15.0):
            add("SY", f"{var}_kick_v0_K{Kk:g}", S5, kind="K", voff=0.0, Kk=Kk, signvar=var)
    for nz in (0, 1):
        add("Y", f"hyst_noise{nz}", (0, 1, 2), kind="hyst", noise=nz)
        add("Y", f"hystwide_noise{nz}", (0, 1, 2), kind="hyst", noise=nz, vlo=-30.0, vhi=10.0, Tramp=200.0)
    return G


_DONE = None


def cost(c):
    global _DONE
    if _DONE is None:
        fp = os.path.join(HERE, "done_ids.txt")
        _DONE = set(open(fp).read().split()) if os.path.exists(fp) else set()
    if c["id"] in _DONE:
        return 0.0
    a = c["arm"]
    if a in ("Y", "SY"):
        if c["kind"] == "hyst":
            return 3000
        return 150 if (c.get("Kk", 0) >= 8 or c["kind"] == "spk") else 10
    if c.get("insult") == "none":
        return 30
    bio = (82 + c["D"]) / SP
    rp = c["reperf"]
    run_prone = rp in ("naive", "ramp10", "ramp30", "ramp60", "pulsed3x10") or rp.startswith("dis:") or rp.startswith("tblk:")
    return 60 * (bio + 2.25) if run_prone else 25 * bio + 40


def shard_split(G, shard):
    fixedA = [c for c in G if c["arm"] in ("N", "T", "SN", "ST")]
    rest = sorted([c for c in G if c["arm"] not in ("N", "T", "SN", "ST")], key=lambda c: -cost(c))
    tot = {"A": sum(cost(c) for c in fixedA), "B": 0.0}
    asg = {"A": [], "B": []}
    for c in rest:
        s = "A" if tot["A"] <= tot["B"] else "B"
        asg[s].append(c); tot[s] += cost(c)
    if shard == "A":
        N_ = [c for c in fixedA if c["arm"] in ("N", "SN")]; T_ = [c for c in fixedA if c["arm"] in ("T", "ST")]
        return N_ + asg["A"], T_, tot
    return asg["B"], [], tot

# ------------------------------------------------------------------------------------------------ worker

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
        if arm in ("N", "SN"):
            out, heavy, light = run_anat(c); ext = ".anat.npz"
        elif arm in ("Y", "SY"):
            out, heavy, light = run_Y(c); ext = ".npz"
        else:
            out, heavy, light = run_cond(c); ext = ".npz"
        os.makedirs(os.path.join(RES, arm), exist_ok=True); os.makedirs(os.path.join(OUT, "results", arm), exist_ok=True)
        if arm in ("N", "SN"):
            np.savez_compressed(os.path.join(RES, arm, cid + ".heavy.npz"), **heavy)
            np.savez_compressed(os.path.join(RES, arm, cid + ext), **light)
        else:
            np.savez_compressed(os.path.join(RES, arm, cid + ext), **heavy)
        np.savez_compressed(os.path.join(OUT, "results", arm, cid + ext), **light)
        jdump(out, os.path.join(OUT, "results", arm, cid + ".json"))
        jdump(out, dj)
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


SMOKE_IDS = ["S_A_base_s100", "S_B_D30_naive_s0", "S_C_D30_reset_s0", "SN_A_anat_D30_s0", "SN_B_anat_D30_s0",
             "SY_A_kick_v0_K12_s0", "F_set07_D180_block40_s0", "X_k10_D30_ramp30_s0"]
SMOKE_T = ["ST_A_D30_tblk_sustain_0.5_s10", "ST_B_D30_tblk_sustain_0.5_s10", "T_D180_tblk_sustain180_15_s10"]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", default="A"); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--arms", default=""); ap.add_argument("--smoke", action="store_true"); ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    G = grid()
    if a.smoke:
        RES = RES + "_smoke"
        allc = {c["id"]: c for c in G}
        if "X_k0_D30_naive_s0" not in allc:
            allc["X_k0_D30_naive_s0"] = dict(arm="X", id="X_k0_D30_naive_s0", name="k0_D30_naive", seed=0, insult="anoxia", D=30, reperf="naive", k_exc=0.0)
        first = [allc[i] for i in SMOKE_IDS]; second = [allc[i] for i in SMOKE_T]
        tot = {}
    else:
        cost({"id": "", "arm": "x"}) if False else None
        fp = os.path.join(HERE, "done_ids.txt")
        if os.path.exists(fp):
            dn = set(open(fp).read().split()); G = [c for c in G if c["id"] not in dn]
        seen = set(); G = [c for c in G if not (c["id"] in seen or seen.add(c["id"]))]
        if a.arms:
            G = [c for c in G if c["arm"] in a.arms.split(",")]
        first, second, tot = shard_split(G, a.shard)
    if a.list:
        from collections import Counter
        print(len(G), Counter(c["arm"] for c in G), "shard", a.shard, len(first), len(second), {k: round(v / 3600, 1) for k, v in tot.items()}); sys.exit()
    os.makedirs(OUT, exist_ok=True); os.makedirs(RES, exist_ok=True)
    t0 = time.time(); setup()
    json.dump(lhs_sets(), open(os.path.join(OUT, "lhs_sets.json"), "w"))
    print("sign variants n_changed:", MV_N, flush=True)
    print(f"shard {a.shard}: first {len(first)} second {len(second)} est CPU-h {tot}; model {time.time()-t0:.1f}s; base n={None if SW.BASE is None else SW.BASE['n']}", flush=True)
    log = open(os.path.join(OUT, f"progress_{a.shard}.log"), "a")
    import multiprocessing as mp_
    pool = mp_.get_context("fork").Pool(a.workers, maxtasksperchild=4)
    def report(res, n=[0]):
        cid, stt, w = res; n[0] += 1
        msg = f"[{time.time()-t0:8.0f}s] {n[0]} {cid} {stt} {w}s"; print(msg, flush=True); log.write(msg + "\n"); log.flush()
    asy_first = [pool.apply_async(worker, (c,), callback=report) for c in first]
    if second:
        for x, c in zip(asy_first, first):
            if c["arm"] in ("N", "SN"):
                x.wait()
        for var in (None, "A", "B"):
            try:
                build_drivers(var); print("drivers built", var, flush=True)
            except Exception:
                print("DRIVER BUILD FAILED", var, "\n" + traceback.format_exc(), flush=True)
        asy_second = [pool.apply_async(worker, (c,), callback=report) for c in second]
    else:
        asy_second = []
    for x in asy_first + asy_second:
        x.wait()
    pool.close(); pool.join()
    print("DONE", round(time.time() - t0), flush=True)
