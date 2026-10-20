"""Round 3 (Modal, adaptive): sign robustness (A/B/C) incl. driver recomputation + targeted-reset bisection,
605-run follow-up (sequential seeding), closed-loop adaptive reset (CL).  One container per --group.

python mech_round3.py --group A|B|smoke --workers W
  group A: sign arms (S, SN, SY, BI chains for def/A/B/C), CL arms (CL, CLL, CLH, CLT)
  group B: F (LHS 40/60-s blocks), X k_exc=10 H3 grid, Y wide hysteresis
Priority queue keeps <= W tasks in flight so sequential bisection steps are not stuck behind bulk work.
Sequential seeding: 3 seeds; +2 seeds (3,4) when the 3 seeds disagree on the decision (rescued / runaway).
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import json, time, glob, math, argparse, traceback, heapq, queue, shutil
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mech_sweep as MS

RES = MS.RES; OUT = MS.OUT
STAGED = os.path.join(os.path.dirname(HERE), "staged")
LOGF = None
FAILING_SETS = [7, 8, 10, 12, 14, 23, 26, 30, 38, 39, 43]          # re-ignition at D180 reset (main sweep)
STRUCT_SETS = []                                                       # filled from staged list
CEIL_DEF = 0.9940


def log(*a):
    msg = " ".join(str(x) for x in a)
    print(msg, flush=True)
    if LOGF:
        LOGF.write(msg + "\n"); LOGF.flush()

# ------------------------------------------------------------------------------------------------ staging / baselines

def stage():
    if not os.path.isdir(STAGED):
        return
    for src in glob.glob(os.path.join(STAGED, "res", "*", "*")):
        arm = os.path.basename(os.path.dirname(src)); dst = os.path.join(RES, arm, os.path.basename(src))
        if not os.path.exists(dst):
            os.makedirs(os.path.dirname(dst), exist_ok=True); shutil.copy(src, dst)


def fii(Rb, Rc):
    return MS.SW.fii_matrix(Rb, Rc)


_BASE = {}


def base_for(key):
    """key: 'def' | ('var', v) | ('set', s) -> (Rb, ceil) or None if not available yet."""
    if key in _BASE:
        return _BASE[key]
    if key == "def":
        RR = [np.load(f)["R"].astype(np.float32) for f in sorted(glob.glob(os.path.join(STAGED, "dynbase", "*.npz")))]
        val = (np.mean(RR, 0), CEIL_DEF)
    else:
        if key[0] == "var":
            fl = [os.path.join(RES, "S", f"S_{key[1]}_base_s{s}.npz") for s in (100, 101, 102, 103, 104)]
        else:
            fl = sorted(glob.glob(os.path.join(STAGED, "lbase", f"L_set{key[1]:02d}_base_s*.npz")))
        if not fl or not all(os.path.exists(f) for f in fl):
            return None
        RR = [np.load(f)["R"].astype(np.float32) for f in fl]
        loo = [fii(np.mean([RR[j] for j in range(len(RR)) if j != i], 0), RR[i]) for i in range(len(RR))]
        val = (np.mean(RR, 0), float(np.median(loo)))
    _BASE[key] = val
    return val


def basekey(c):
    if c.get("signvar"):
        return ("var", c["signvar"])
    if c.get("lhs") is not None:
        return ("set", int(c["lhs"]["set"]))
    return "def"


def outcome(c):
    """(FII_norm, runaway, rescued) of a finished condition."""
    o = json.load(open(os.path.join(RES, c["arm"], c["id"] + ".json")))
    run = bool((o.get("obs_late_Hz") or 0) > 0.5 or (o.get("FSA") or 0) > 0.05)
    b = base_for(basekey(c))
    if b is None:
        return None, run, False
    R = np.load(os.path.join(RES, c["arm"], c["id"] + ".npz"))["R"].astype(np.float32)
    fn = fii(b[0], R) / b[1]
    return fn, run, bool(fn >= 0.8 and not run)

# ------------------------------------------------------------------------------------------------ scheduler


class Sched:
    def __init__(self, W):
        import multiprocessing as mp_
        self.pool = mp_.get_context("fork").Pool(W, maxtasksperchild=4)
        self.W = W; self.heap = []; self.seq = 0; self.inflight = 0; self.Q = queue.Queue()
        self.done = set(); self.status = {}; self.hooks = []; self.t0 = time.time(); self.nsub = 0

    def submit(self, c, pri=2):
        if c["id"] in self.done or c["id"] in self.status:
            return
        self.status[c["id"]] = "queued"
        if os.path.exists(os.path.join(RES, c["arm"], c["id"] + ".json")):
            self.Q.put((c["id"], "skip", 0)); self.inflight += 1; return
        heapq.heappush(self.heap, (pri, self.seq, c)); self.seq += 1

    def hook(self, ids, fn):
        self.hooks.append((set(ids), fn))

    def _fill(self):
        while self.heap and self.inflight < self.W:
            _, _, c = heapq.heappop(self.heap); self.inflight += 1; self.nsub += 1
            self.pool.apply_async(MS.worker, (c,), callback=self.Q.put,
                                  error_callback=lambda e, cid=c["id"]: self.Q.put((cid, "error:" + repr(e)[:200], 0)))

    def run(self):
        self._fill()
        while self.inflight > 0 or self.heap:
            cid, stt, w = self.Q.get(); self.inflight -= 1
            self.done.add(cid); self.status[cid] = stt
            log(f"[{time.time()-self.t0:7.0f}s] {len(self.done)} {cid} {stt} {w}s q={len(self.heap)} infl={self.inflight}")
            fired = [h for h in self.hooks if h[0] <= self.done]
            self.hooks = [h for h in self.hooks if not h[0] <= self.done]
            for _, fn in fired:
                try:
                    fn()
                except Exception:
                    log("HOOK ERROR\n" + traceback.format_exc())
            self._fill()
        self.pool.close(); self.pool.join()

# ------------------------------------------------------------------------------------------------ condition factories


def C_(arm, name, seed, **kw):
    return dict(arm=arm, id=f"{arm}_{name}_s{seed}", name=name, seed=seed, **kw)


def seq_seed(S, make, pri, tag, seeds0=(0, 1, 2), extra=(3, 4), crit="rescued"):
    """Submit seeds0; when done, add `extra` seeds if seeds disagree on crit."""
    cs = [make(s) for s in seeds0]
    for c in cs:
        S.submit(c, pri)

    def chk():
        oc = [outcome(c) for c in cs]
        vals = [o[2] if crit == "rescued" else o[1] for o in oc]
        if 0 < sum(vals) < len(vals):
            log("TOPUP", tag, vals)
            for s in extra:
                S.submit(make(s), pri)
    S.hook([c["id"] for c in cs], chk)


BIS = {}


class Chain:
    """Log-bisection for the minimal driver ('sustain') set rescuing >=4/5 held-out seeds."""

    def __init__(self, S, var, D, lo=None, hi=None, steps=4):
        self.S, self.var, self.D, self.lo, self.hi, self.steps = S, var, D, lo, hi, steps
        self.res = {}; self.tag = f"{var or 'def'}_D{D}"; BIS[self.tag] = dict(var=var, D=D, tested={}, bracket=None, note="")

    def conds(self, k):
        k = float(f"{k:.3g}")
        return k, [C_("BI", f"{self.var or 'def'}_D{self.D}_k{k:g}", s, insult="anoxia", D=self.D, reperf=f"tblk:sustain:{k:g}",
                      tkind="sustain", tk=k, signvar=self.var) for s in (10, 11, 12, 13, 14)]

    def evaluate(self, ks, then):
        allc = []
        for k in ks:
            k2, cs = self.conds(k); allc += cs
            for c in cs:
                self.S.submit(c, 0)

        def done():
            for k in ks:
                k2, cs = self.conds(k)
                oc = [outcome(c) for c in cs]
                nres = sum(o[2] for o in oc); self.res[k2] = nres >= 4
                BIS[self.tag]["tested"][f"{k2:g}"] = dict(n_rescued=int(nres), FII_norm=[None if o[0] is None else round(o[0], 3) for o in oc])
                log("BISECT", self.tag, k2, nres)
            then()
        self.S.hook([c["id"] for c in allc], done)

    def start(self):
        if self.lo is not None:
            return self.step()
        self.evaluate([0.1, 20.0], self.after_init)

    def after_init(self):
        if self.res[0.1]:
            self.hi = 0.1; self.lo = 0.02; self.evaluate([0.02], self.after_low)
        elif not self.res[20.0]:
            self.evaluate([50.0], self.after_high)
        else:
            self.lo, self.hi = 0.1, 20.0; self.step()

    def after_low(self):
        if self.res[0.02]:
            BIS[self.tag]["bracket"] = [None, 0.02]; BIS[self.tag]["note"] = "rescued at smallest tested set"; self.save()
        else:
            self.lo = 0.02; self.step()

    def after_high(self):
        if self.res[50.0]:
            self.lo, self.hi = 20.0, 50.0; self.step()
        else:
            BIS[self.tag]["bracket"] = [50.0, None]; BIS[self.tag]["note"] = "no rescue up to 50%"; self.save()

    def step(self):
        if self.steps <= 0 or self.hi / self.lo < 1.3:
            BIS[self.tag]["bracket"] = [self.lo, self.hi]; self.save(); return
        self.steps -= 1
        m = float(f"{math.sqrt(self.lo * self.hi):.3g}")
        def upd(m=m):
            if self.res[m]:
                self.hi = m
            else:
                self.lo = m
            self.step()
        self.evaluate([m], upd)

    def save(self):
        json.dump(BIS, open(os.path.join(OUT, "bisection.json"), "w"), indent=1)
        log("BISECT DONE", self.tag, BIS[self.tag]["bracket"])

# ------------------------------------------------------------------------------------------------ groups


def group_A(S):
    S5 = (0, 1, 2, 3, 4)
    for v in ("A", "B", "C"):
        for s in (100, 101, 102, 103, 104):
            S.submit(C_("S", f"{v}_base", s, insult="none", signvar=v), 0)
        for D in (30, 180):
            for s in S5:
                S.submit(C_("SN", f"{v}_anat_D{D}", s, insult="anoxia", D=D, reperf="naive", signvar=v), 0)
                for rp in ("naive", "reset"):
                    S.submit(C_("S", f"{v}_D{D}_{rp}", s, insult="anoxia", D=D, reperf=rp, signvar=v), 1)
        for Kk in (0.0, 9.0, 10.0, 12.0, 15.0):
            for s in S5:
                S.submit(C_("SY", f"{v}_kick_v0_K{Kk:g}", s, kind="K", voff=0.0, Kk=Kk, signvar=v), 1)
        ids = [f"SN_{v}_anat_D{D}_s{s}" for D in (30, 180) for s in S5] + [f"S_{v}_base_s{s}" for s in (100, 101, 102, 103, 104)]

        def start_v(v=v):
            MS.build_drivers(v); log("drivers built", v)
            for D in (30, 180):
                Chain(S, v, D).start()
        S.hook(ids, start_v)
    # default-sign refinement of the minimal set (held-out seeds 10-14; brackets from main sweep)
    Chain(S, None, 30, lo=0.1, hi=0.5, steps=3).start()
    Chain(S, None, 180, lo=10.0, hi=20.0, steps=3).start()
    # closed-loop adaptive reset (i) calibrated set
    for D in (30, 180, 360):
        for T in (17.0, 22.0, 29.0):
            for rp in ("reset", "block40", "CL"):
                seq_seed(S, lambda s, D=D, T=T, rp=rp: C_("CL", f"cal_D{D}_T{T:g}_{rp}", s, insult="anoxia", D=D, T=T, reperf=rp), 1,
                         f"CL_cal_D{D}_T{T:g}_{rp}")
    # (ii) LHS sets where the fixed 20-s block failed at D180
    sets = {p["set"]: p for p in MS.lhs_sets()}
    for ss in FAILING_SETS + STRUCT_SETS:
        lc = MS.lhs_cond(sets[ss])
        seq_seed(S, lambda s, ss=ss, lc=lc: C_("CLL", f"set{ss:02d}_D180_CL", s, insult="anoxia", D=180, reperf="CL", **lc), 1, f"CLL_{ss}")
    # (iii) with hub protection top-1 %
    for D in (180, 360):
        for rp in ("reset", "CL"):
            seq_seed(S, lambda s, D=D, rp=rp: C_("CLH", f"hub1_D{D}_{rp}", s, insult="anoxia", D=D, reperf=rp, protect=["hub", 1]), 1, f"CLH_{D}_{rp}")
    # closed-loop targeted (driver core only, until criterion); held-out seeds 10-12
    for D in (30, 180):
        for k in (0.5, 2.0, 5.0, 20.0):
            seq_seed(S, lambda s, D=D, k=k: C_("CLT", f"D{D}_CL_sustain_{k:g}", s, insult="anoxia", D=D, reperf=f"CL:sustain:{k:g}"), 1,
                     f"CLT_{D}_{k}", seeds0=(10, 11, 12), extra=(13, 14))


def group_B(S):
    sets = MS.lhs_sets()
    for p in sets:
        lc = MS.lhs_cond(p)
        for rp in ("block40", "block60"):
            seq_seed(S, lambda s, p=p, lc=lc, rp=rp: C_("F", f"set{p['set']:02d}_D180_{rp}", s, insult="anoxia", D=180, reperf=rp, **lc), 2,
                     f"F_{p['set']}_{rp}")
    for D in (30, 180, 360):
        for pr in MS.X_PROT:
            seq_seed(S, lambda s, D=D, pr=pr: C_("X", f"k10_D{D}_{pr}", s, insult="anoxia", D=D, reperf=pr, mx=dict(k_exc=10.0), k_exc=10.0), 2,
                     f"X_k10_{D}_{pr}")
    for nz in (0, 1):
        for s in (0, 1, 2):
            S.submit(C_("Y", f"hystwide_noise{nz}", s, kind="hyst", noise=nz, vlo=-30.0, vhi=10.0, Tramp=200.0), 0)


def group_smoke(S):
    S.submit(C_("S", "A_base", 100, insult="none", signvar="A"), 0)
    S.submit(C_("SY", "C_kick_v0_K12", 0, kind="K", voff=0.0, Kk=12.0, signvar="C"), 0)
    S.submit(C_("CL", "cal_D30_T22_CL", 0, insult="anoxia", D=30, T=22.0, reperf="CL"), 0)
    S.submit(C_("CLH", "hub1_D180_CL", 0, insult="anoxia", D=180, reperf="CL", protect=["hub", 1]), 0)
    S.submit(C_("CLT", "D30_CL_sustain_2", 10, insult="anoxia", D=30, reperf="CL:sustain:2"), 0)
    sets = {p["set"]: p for p in MS.lhs_sets()}
    S.submit(C_("F", "set07_D180_block40", 0, insult="anoxia", D=180, reperf="block40", **MS.lhs_cond(sets[7])), 0)
    ch = Chain(S, None, 30, lo=0.1, hi=0.5, steps=1); ch.start()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--group", default="A"); ap.add_argument("--workers", type=int, default=32)
    a = ap.parse_args()
    if a.group == "smoke":
        RES = MS.RES = RES + "_smoke"
    os.makedirs(OUT, exist_ok=True); os.makedirs(RES, exist_ok=True)
    LOGF = open(os.path.join(OUT, f"round3_{a.group}.log"), "a")
    sl = os.path.join(STAGED, "struct_sets.json")
    if os.path.exists(sl):
        STRUCT_SETS = json.load(open(sl))
    t0 = time.time(); MS.setup(); stage()
    log("setup", round(time.time() - t0, 1), "variants", MS.MV_N, "staged files", len(glob.glob(os.path.join(RES, "*", "*"))))
    S = Sched(a.workers)
    {"A": group_A, "B": group_B, "smoke": group_smoke}[a.group](S)
    S.run()
    json.dump(BIS, open(os.path.join(OUT, "bisection.json"), "w"), indent=1)
    for v in ("A", "B", "C"):
        f = os.path.join(RES, f"drivers_{v}.npz")
        if os.path.exists(f):
            shutil.copy(f, OUT)
    log("DONE", round(time.time() - t0))
