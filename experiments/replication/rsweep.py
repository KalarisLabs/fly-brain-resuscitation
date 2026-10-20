"""Replication sweep (MCNS v1.0 / BANC v888 / FAFB-proxy control) — same protocol as dynamics/sweep.py.

python rsweep.py --ds MCNS --wsyn 0.123 --workers 32 [--shard k --nshard n] [--arms ...] [--list]
Protocol (slow-layer s, speedup 50): 2 s normoxia -> insult D s -> 60 s reperfusion protocol -> 20 s free observation
-> evaluation at speedup 1 (500-ms spontaneous window + 7 probes x 250 ms from copies of the end state).
Calibrated metabolic set as FAFB: tau_atp 30 s, J50 600 s, q10_injury 2.5; demand proxy = total synapse count.
Results: RES_DIR/<DS>/<arm>/<id>.json (+ heavy .npz) on the Volume; light copies under OUT_DIR/results/<DS>/<arm>/.
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
import json, time, copy, argparse, traceback, re
import numpy as np, pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flysim_cns as fc
fs = fc.fs

RES = os.environ.get("RES_DIR", os.path.join(_P.OUT, "followup", "replication", "raw", "results"))
OUT = os.environ.get("OUT_DIR", "out")
SP = 50.0
CAL = dict(tau_atp_s=30.0, J50_s=600.0, q10_injury=2.5)
REPERF_S = 60.0; OBS_S = 20.0; PRE_S = 2.0
PROBE_MS = 250.0; SPONT_MS = 500.0
DS = None; LABEL = None; M = None; POP = {}; AT = None; DEG = None; WSYN = None; REGNAMES = None


def setup(ds, wsyn, alln=False):
    global DS, M, POP, AT, DEG, WSYN, REGNAMES
    DS, WSYN = ds, wsyn
    if ds == "FAFBproxy":
        M = fc.fafb_proxy(w_syn=wsyn, alln_inhib=alln)
        dn = M.neurons.dominant_neuropil.astype(str)
        M.neurons["region"] = np.where(dn.str.match(r"^(ME|LO|LOP|AME|LA)(_|$)"), "OL", "CB")
        bat = fc.battery(M, "FAFB")
    else:
        M = fc.load_cns(ds, w_syn=wsyn, alln_inhib=alln)
        M.neurons["region"] = M.neurons.region.fillna("NA")
        bat = fc.battery(M, ds)
    CAL["demand_mode"] = "syn_total"
    P, ro, named, rep = bat
    POP["probes"] = P; POP["ro"] = ro; POP["named"] = named; POP["report"] = rep
    POP["readouts"] = ro["all"]
    POP["ro_is_brain"] = np.isin(ro["all"], ro["brain"]); POP["ro_is_vnc"] = np.isin(ro["all"], ro["vnc"])
    A = sp.csr_matrix((np.ones(len(M.indices), np.float32), M.indices, M.indptr), shape=(M.N, M.N))
    AT = A.T.tocsr()
    DEG = np.diff(M.indptr) + np.bincount(M.indices, minlength=M.N)
    reg = M.neurons.region.astype(str).values
    REGNAMES = ["OL", "CB", "VNC"]
    POP["reg"] = {r: reg == r for r in REGNAMES}


def reach(alive, ro):
    tot = 0; hit = 0
    for name, (idx, _) in POP["probes"].items():
        x = np.zeros(M.N, np.float32); x[idx[alive[idx]]] = 1.0
        seen = x > 0
        for _ in range(6):
            x = (AT @ x > 0) & alive & ~seen
            seen |= x
            x = x.astype(np.float32)
            if not x.any():
                break
        ok = seen[ro] & alive[ro]
        hit += ok.sum(); tot += ro.size
    return float(hit / max(tot, 1))


def wrap(f, ts0):
    return (lambda ts, _f=f, _o=ts0: _f(ts - _o)) if callable(f) else f


def phase(st, dur_s, sp_, metab_over, supply=1.0, temp=22.0, iv=None, stim=None, seed=0):
    ts0 = (st.t_ms * 1e-3 * sp_) if st is not None else 0.0
    mx = dict(metab_over); mx["speedup"] = sp_
    ivw = {k: (wrap(v, ts0) if callable(v) else v) for k, v in (iv or {}).items()}
    return fs.run(M, dur_s * 1000.0 / sp_, stim=stim, metab=mx, supply=wrap(supply, ts0), temperature=wrap(temp, ts0),
                  interventions=ivw, state=st, groups="region", rec_every_ms=500.0 / sp_, seed=seed)


class Trace:
    def __init__(self):
        self.t0 = 0.0
        self.d = {k: [] for k in ["t", "e", "dep", "blk", "dead", "Km", "Kx", "syn", "J", "pop_t", "pop"] +
                  [f"pop_{r}" for r in ["OL", "CB", "VNC"]]}

    def add(self, r, dur_s, sp_):
        sl = r["slow"]; ts = np.asarray(sl["t_s"])
        if ts.size:
            tt = self.t0 + (ts - ts[0]); K = np.asarray(sl["K"])
            for k, v in [("t", tt), ("e", sl["e_mean"]), ("dep", sl["dep_frac"]), ("blk", sl["block_frac"]),
                         ("dead", sl["dead_frac"]), ("Km", K.mean(1)), ("Kx", K.max(1)), ("syn", sl["syn_frac"]),
                         ("J", sl["J_mean"])]:
                self.d[k].append(np.asarray(v, np.float32))
        pm = r["pop_ms"].astype(np.float64); gn = list(r["group_names"])
        b = max(1, int(round(500.0 / sp_))); nb = pm.shape[0] // b
        if nb:
            tot = pm.sum(1)[: nb * b].reshape(nb, b).sum(1) / (b * 1e-3) / M.N
            self.d["pop_t"].append((self.t0 + 0.5 * np.arange(nb)).astype(np.float32)); self.d["pop"].append(tot.astype(np.float32))
            for rg in ["OL", "CB", "VNC"]:
                nreg = max(int(POP["reg"][rg].sum()), 1)
                x = pm[:, gn.index(rg)] if rg in gn else np.zeros(pm.shape[0])
                self.d[f"pop_{rg}"].append((x[: nb * b].reshape(nb, b).sum(1) / (b * 1e-3) / nreg).astype(np.float32))
        self.t0 += dur_s

    def arrays(self):
        return {k: (np.concatenate(v) if v else np.zeros(0, np.float32)) for k, v in self.d.items()}


def piece(pts):
    return fs.piecewise(pts)


def reperf_protocol(name, T_rep):
    S = lambda dur, **kw: dict(dict(dur=dur, supply=1.0, temp=T_rep, iv={}, revive=False), **kw)
    blk = lambda on, off: {"block_all": piece([(0, False), (on, True), (on + off, False)])}
    if name == "naive":
        return [S(REPERF_S)]
    if name == "reset":
        return [S(REPERF_S, iv=blk(0, 20))]
    raise ValueError(name)


def evaluate(st, seed, mx):
    mxe = dict(mx); mxe["speedup"] = 1.0
    r = fs.run(M, SPONT_MS, metab=mxe, supply=1.0, temperature=22.0, state=st, groups=None, seed=seed + 1000, rec_every_ms=1e9)
    spont = r["counts"][0].copy(); st = r["state"]; alive = st.alive.copy()
    pc = []
    for j, (name, (idx, rate)) in enumerate(POP["probes"].items()):
        s2 = copy.deepcopy(st)
        rp = fs.run(M, PROBE_MS, stim=[dict(idx=idx, rate=rate)], metab=mxe, supply=1.0, temperature=22.0, state=s2,
                    groups=None, seed=seed * 100 + 5000 + j, rec_every_ms=1e9)
        pc.append(rp["counts"][0].copy())
    pc = np.array(pc)
    spont_rate = spont / (SPONT_MS * 1e-3); probe_rate = pc / (PROBE_MS * 1e-3)
    ro = POP["readouts"]
    Rm = probe_rate[:, ro] - spont_rate[ro][None, :]
    return dict(spont=spont, pc=pc, spont_rate=spont_rate, probe_rate=probe_rate, R=Rm.astype(np.float32), alive=alive)


def run_condition(c):
    t0 = time.time(); seed = c["seed"]; spc = SP
    mx = dict(CAL); mx.update(c.get("mx", {}))
    T_ins = c.get("T_ins", 22.0); T_rep = c.get("T_rep", T_ins)
    extra = {}; prot = {}
    if c.get("wscale"):
        mx = dict(mx)
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
        extra["protected_region_counts"] = {r: int(POP["reg"][r][sel].sum()) for r in POP["reg"]}
    tr = Trace(); mark = {}
    if c.get("insult") == "none":
        st = fs.init_state(M, dict(fs.METAB_DEFAULTS, **mx), seed)
        ev = evaluate(st, seed, mx)
        return finish(c, ev, tr, mark, t0, extra)
    r = phase(None, PRE_S, spc, mx, temp=T_ins, seed=seed); tr.add(r, PRE_S, spc); st = r["state"]
    D = float(c["D"])
    mark["insult_start"] = tr.t0
    r = phase(st, D, spc, mx, supply=0.0, temp=T_ins, iv=prot, seed=seed); tr.add(r, D, spc); st = r["state"]
    mark["reperf_start"] = tr.t0
    extra["dead_after_insult"] = float(1 - st.alive.mean())
    if c["reperf"] == "closedloop":
        # closed-loop reset: total block in 1-s slow steps until max compartment K+ < K_thr and mean ATP (alive) > E_thr
        K_thr = float(c.get("K_thr", 6.0)); E_thr = float(c.get("E_thr", 0.9)); tb = 0.0
        while tb < REPERF_S:
            r = phase(st, 1.0, spc, mx, supply=1.0, temp=T_rep, iv=dict(prot, block_all=True), seed=seed)
            tr.add(r, 1.0, spc); st = r["state"]; tb += 1.0
            e_ = st.a / st.a0
            if float(np.max(st.K)) < K_thr and float(e_[st.alive].mean()) > E_thr:
                break
        extra["cl_block_s"] = tb; extra["cl_Kmax_release"] = float(np.max(st.K))
        extra["cl_e_release"] = float((st.a / st.a0)[st.alive].mean())
        if REPERF_S - tb > 0:
            r = phase(st, REPERF_S - tb, spc, mx, supply=1.0, temp=T_rep, iv=dict(prot), seed=seed)
            tr.add(r, REPERF_S - tb, spc); st = r["state"]
    else:
      for sub in reperf_protocol(c["reperf"], T_rep):
        iv = dict(prot); iv.update(sub["iv"])
        r = phase(st, sub["dur"], spc, mx, supply=sub["supply"], temp=sub["temp"], iv=iv, seed=seed)
        tr.add(r, sub["dur"], spc); st = r["state"]
    mark["obs_start"] = tr.t0
    r = phase(st, OBS_S, spc, mx, supply=1.0, temp=T_rep, seed=seed); tr.add(r, OBS_S, spc); st = r["state"]
    mark["end"] = tr.t0
    ev = evaluate(st, seed, mx)
    return finish(c, ev, tr, mark, t0, extra)


def finish(c, ev, tr, mark, t0, extra):
    A = tr.arrays(); alive = ev["alive"]
    out = dict(c); out.update(extra); out["ds"] = DS; out["label"] = LABEL; out["w_syn"] = WSYN
    out["dead_frac"] = float(1 - alive.mean())
    out["dead_frac_region"] = {r: float(1 - alive[m].mean()) for r, m in POP["reg"].items() if m.any()}
    out["reach_brain"] = reach(alive, POP["ro"]["brain"])
    out["reach_vnc"] = reach(alive, POP["ro"]["vnc"]) if len(POP["ro"]["vnc"]) else None
    sr = ev["spont_rate"]
    out["FSA"] = float((ev["spont"] > 0).mean())
    out["FSA_region"] = {r: float((ev["spont"][m] > 0).mean()) for r, m in POP["reg"].items() if m.any()}
    out["spont_Hz_per_neuron"] = float(sr.mean())
    out["spont_Hz_region"] = {r: float(sr[m].mean()) for r, m in POP["reg"].items() if m.any()}
    evoked = (ev["probe_rate"] >= sr[None, :] + 8.0)
    out["n_evoked"] = int(evoked.any(0).sum())
    out["n_evoked_region"] = {r: int(evoked[:, m].any(0).sum()) for r, m in POP["reg"].items() if m.any()}
    R = ev["R"]
    out["readout_resp_max_brain"] = float(R[:, POP["ro_is_brain"]].max()) if POP["ro_is_brain"].any() else None
    out["readout_resp_max_vnc"] = float(R[:, POP["ro_is_vnc"]].max()) if POP["ro_is_vnc"].any() else None
    ro = POP["readouts"]; pos = {int(i): k for k, i in enumerate(ro)}
    nm = {}
    for name, idx in POP["named"].items():
        ii = [pos[int(i)] for i in idx if int(i) in pos]
        nm[name] = [round(float(np.mean(R[j, ii])), 2) if ii else None for j in range(R.shape[0])]
    out["named_resp"] = nm; out["probe_names"] = list(POP["probes"].keys())
    if A["t"].size and mark:
        t = A["t"]; ins0 = mark["insult_start"]; rp0 = mark["reperf_start"]
        w = (t >= ins0) & (t < rp0); dep = A["dep"]
        on = np.flatnonzero(w & (dep > 0.5)); out["SD_onset_s"] = float(t[on[0]] - ins0) if on.size else None
        out["K_peak"] = float(A["Kx"].max()); out["Kmean_peak"] = float(A["Km"].max()); out["e_min"] = float(A["e"].min())
        rep = np.flatnonzero((t >= rp0) & (dep < 0.01))
        out["repol_s_after_reperf"] = float(t[rep[0]] - rp0) if rep.size else None
        out["K_mean_end"] = float(A["Km"][-1])
        pt = A["pop_t"]; wr2 = (pt >= rp0) & (pt < mark["obs_start"]); wo = pt >= mark["end"] - 10.0
        for key, arr in [("", A["pop"])] + [(f"_{r}", A[f"pop_{r}"]) for r in ["OL", "CB", "VNC"]]:
            if arr.size == pt.size and arr.size:
                out["surge_peak_Hz" + key] = float(arr[wr2].max()) if wr2.any() else None
                out["obs_late_Hz" + key] = float(arr[wo].mean()) if wo.any() else None
        # region ignition latency after reperfusion / block release: first bin > 1 Hz/neuron
        for r in ["OL", "CB", "VNC"]:
            arr = A[f"pop_{r}"]
            if arr.size == pt.size:
                k = np.flatnonzero((pt >= rp0) & (arr > 1.0))
                out[f"ignite_s_{r}"] = float(pt[k[0]] - rp0) if k.size else None
        out["marks"] = mark
    out["wall_s"] = round(time.time() - t0, 1)
    try:
        import resource; out["maxrss_MB"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
    except Exception:
        pass
    cnt = np.vstack([ev["spont"][None, :], ev["pc"]])
    sp_i = [np.flatnonzero(x).astype(np.int32) for x in cnt]
    sp_c = [np.minimum(x[i], 65535).astype(np.uint16) for x, i in zip(cnt, sp_i)]
    heavy = dict(R=R, alive=np.packbits(alive), probe_rate=(ev["probe_rate"].astype(np.float32) if c["arm"] == "base" else np.zeros(0)),
                 **{f"trace_{k}": v for k, v in A.items()},
                 **{f"cnt_i{j}": a for j, a in enumerate(sp_i)}, **{f"cnt_c{j}": a for j, a in enumerate(sp_c)})
    light = dict(R=R.astype(np.float16), alive=np.packbits(alive), ro=POP["readouts"].astype(np.int32),
                 ro_is_brain=POP["ro_is_brain"], ro_is_vnc=POP["ro_is_vnc"],
                 **{f"trace_{k}": v for k, v in A.items()})
    if c["arm"] == "base":
        light["probe_rate"] = ev["probe_rate"].astype(np.float16)
    return out, heavy, light


def jdump(o, path):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(o, f, default=lambda x: (x.tolist() if hasattr(x, "tolist") else str(x)))
    os.replace(tmp, path)


def worker(c):
    arm, cid = c["arm"], c["id"]
    dj = os.path.join(RES, LABEL, arm, cid + ".json")
    if os.path.exists(dj):
        return cid, "skip", 0
    try:
        out, heavy, light = run_condition(c)
        os.makedirs(os.path.join(RES, LABEL, arm), exist_ok=True); os.makedirs(os.path.join(OUT, "results", LABEL, arm), exist_ok=True)
        np.savez_compressed(os.path.join(RES, LABEL, arm, cid + ".npz"), **heavy)
        np.savez_compressed(os.path.join(OUT, "results", LABEL, arm, cid + ".npz"), **light)
        jdump(out, os.path.join(OUT, "results", LABEL, arm, cid + ".json"))
        jdump(out, dj)
        return cid, "ok", out["wall_s"]
    except Exception:
        tb = traceback.format_exc()
        os.makedirs(os.path.join(OUT, "errors"), exist_ok=True)
        open(os.path.join(OUT, "errors", f"{LABEL}_{cid}.txt"), "w").write(tb)
        return cid, "error", 0


def grid(ds, variant="main"):
    G = []
    S5 = (0, 1, 2, 3, 4)

    def add(arm, name, seeds=S5, **kw):
        for s in seeds:
            G.append(dict(arm=arm, id=f"{arm}_{name}_s{s}", name=name, seed=s, **kw))
    if variant == "alln":       # sign-corrected (antennal-lobe local neurons inhibitory): one core contrast
        for s in range(5):
            G.append(dict(arm="base", id=f"base_s{s}", name="base", seed=s, insult="none"))
        for D in (30, 180):
            for rp in ("naive", "reset"):
                add("core", f"D{D}_T22_{rp}", insult="anoxia", D=D, reperf=rp)
        return G
    if variant == "wlo":        # w_syn sensitivity: 0.7 x w_crit_any, 3 seeds
        for s in range(3):
            G.append(dict(arm="base", id=f"base_s{s}", name="base", seed=s, insult="none"))
        for D in (30, 180):
            for rp in ("naive", "reset"):
                add("core", f"D{D}_T22_{rp}", seeds=(0, 1, 2), insult="anoxia", D=D, reperf=rp)
        return G
    nb = 10 if ds != "FAFBproxy" else 5
    for s in range(nb):
        G.append(dict(arm="base", id=f"base_s{s}", name="base", seed=s, insult="none"))
    if ds == "FAFBproxy":
        for D in (30, 180):
            for rp in ("naive", "reset"):
                add("core", f"D{D}_T22_{rp}", insult="anoxia", D=D, reperf=rp)
        add("hub", "D360_none_reset", insult="anoxia", D=360, reperf="reset")
        add("hub", "D360_hub1_reset", insult="anoxia", D=360, reperf="reset", protect=("hub", 1))
        add("cl", "D180_T22_closedloop", insult="anoxia", D=180, reperf="closedloop")
        return G
    add("core", "sham_reset", seeds=(0, 1, 2), insult="anoxia", D=0.0, reperf="reset")
    add("core", "sham_naive", seeds=(0, 1, 2), insult="anoxia", D=0.0, reperf="naive")
    for D in (30, 180):
        for rp in ("naive", "reset"):
            add("core", f"D{D}_T22_{rp}", insult="anoxia", D=D, reperf=rp)
    for T in (17.0, 29.0):
        for rp in ("naive", "reset"):
            add("temp", f"D180_T{int(T)}_{rp}", insult="anoxia", D=180, T_ins=T, reperf=rp)
    for rp in ("naive", "reset"):
        add("hub", f"D360_none_{rp}", insult="anoxia", D=360, reperf=rp)
        for k in (1, 10):
            add("hub", f"D360_hub{k}_{rp}", insult="anoxia", D=360, reperf=rp, protect=("hub", k))
    for k in (1, 10):
        add("hub", f"D360_rand{k}_reset", insult="anoxia", D=360, reperf="reset", protect=("rand", k))
    add("hub", "D360_whole_reset", insult="anoxia", D=360, reperf="reset", protect=("whole", 100))
    add("cl", "D180_T22_closedloop", insult="anoxia", D=180, reperf="closedloop")
    return G


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ds", required=True); ap.add_argument("--wsyn", type=float, required=True)
    ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--arms", default="")
    ap.add_argument("--ids", default=""); ap.add_argument("--list", action="store_true")
    ap.add_argument("--variant", default="main")
    a = ap.parse_args()
    G = grid(a.ds, a.variant)
    LBL = a.ds if a.variant == "main" else f"{a.ds}_{a.variant}"
    if a.list:
        from collections import Counter
        print(len(G), Counter(c["arm"] for c in G)); sys.exit()
    if a.arms:
        G = [c for c in G if c["arm"] in a.arms.split(",")]
    if a.ids:
        G = [c for c in G if c["id"] in a.ids.split(",")]
    G = G[a.shard::a.nshard]
    # longest first (D360, then D180) for better packing
    G.sort(key=lambda c: -float(c.get("D", 0) or 0))
    os.makedirs(OUT, exist_ok=True)
    DS = a.ds; LABEL = LBL
    todo = [c for c in G if not os.path.exists(os.path.join(RES, LABEL, c["arm"], c["id"] + ".json"))]
    print(f"{a.ds} w_syn={a.wsyn} shard {a.shard}/{a.nshard}: {len(G)} conditions, {len(todo)} to run, workers {a.workers}", flush=True)
    t0 = time.time(); setup(a.ds, a.wsyn, alln=(a.variant == "alln"))
    os.makedirs(os.path.join(OUT, "results", LABEL), exist_ok=True)
    json.dump(dict(report=POP["report"], w_syn=a.wsyn, N=M.N, info=M.info, CAL=CAL, deg_top1_min=int(np.sort(DEG)[-int(round(0.01 * M.N))])),
              open(os.path.join(OUT, "results", LABEL, "setup.json"), "w"), default=str, indent=1)
    print(f"model loaded {time.time()-t0:.1f}s", flush=True)
    log = open(os.path.join(OUT, f"progress_{LABEL}_shard{a.shard}.log"), "a")
    if a.workers <= 1:
        it = map(worker, todo)
    else:
        import multiprocessing as mp_
        pool = mp_.get_context("fork").Pool(a.workers, maxtasksperchild=4)
        it = pool.imap_unordered(worker, todo, chunksize=1)
    n = 0
    for cid, stt, w in it:
        n += 1
        msg = f"[{time.time()-t0:8.0f}s] {n}/{len(todo)} {a.ds} {cid} {stt} {w}s"
        print(msg, flush=True); log.write(msg + "\n"); log.flush()
    print("DONE", a.ds, time.time() - t0, flush=True)
