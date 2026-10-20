"""w_syn calibration for MCNS / BANC (and FAFB as procedure check), FAFB criterion:
a 150-ms sugar-GRN stimulus (150 Hz) must not reverberate (no spikes anywhere in 1000-1200 ms), probes graded.

python calib.py --workers 32 [--ds FAFB,MCNS,BANC]
Stage 1: coarse log grid (x1.15), 3 seeds.  Stage 2: fine grid (x1.02) between last pass and first fail, 3 seeds.
w_crit = smallest fine-grid w at which any seed reverberates; w_pass = largest w below it with all seeds silent.
Stage 3: probe battery (250 ms, all 7 probes + sugar 50/100/200 Hz, 1 seed) at w = r * w_crit, r in RATIOS, and at
FAFB-matched w = w_crit(DS) * 0.185 / w_crit(FAFB) (computed in-job).  Writes out/calib/*.json + calib_summary.json.
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
import json, time, copy, argparse, traceback
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flysim_cns as fc
fs = fc.fs

OUT = os.path.join(os.environ.get("OUT_DIR", "out"), "calib")
MODELS = {}; BAT = {}
RANGES = {"FAFB": (0.10, 0.40), "MCNS": (0.04, 0.35), "BANC": (0.15, 1.0)}
SEEDS1 = (0, 1)
RATIOS = [0.6, 0.7, 0.8, 0.9, 1.0]
SEEDS = (0, 1, 2)
PROBES = ["sugar", "bitter", "JO_aud", "JO_wg", "LC4", "LC9", "turnL"]
EXPLODE = 60000          # active neurons in 200-400 ms above which the run is stopped (counted as reverberating)


def setup(dsl):
    for ds in dsl:
        m = fs.load_model("princeton5") if ds == "FAFB" else fc.load_cns(ds)
        MODELS[ds] = m
        BAT[ds] = fc.battery(m, ds)


def mw(ds, w):
    m = copy.copy(MODELS[ds]); m.w_syn = float(w); return m


def rev_task(ds, w, seed, probe="sugar"):
    m = mw(ds, w); P, ro, named, _ = BAT[ds]
    sug, prate = P[probe]; t0 = time.time()
    r1 = fs.run(m, 400, stim=[dict(idx=sug, rate=prate, t_off=150.0)], epochs=[(0, 150), (200, 400)], groups=None, seed=seed)
    a = [(r1["counts"][k] > 0).sum() for k in range(2)]
    mn9 = float(r1["counts"][0][named["MN9"]].mean() / 0.15) if len(named["MN9"]) else float("nan")
    out = dict(ds=ds, w=w, seed=seed, probe=probe, act_0_150=int(a[0]), act_200_400=int(a[1]), MN9_Hz=mn9)
    if a[1] > EXPLODE:
        out.update(act_400_1000=None, act_1000_1200=None, exploded=True)
    else:
        r2 = fs.run(m, 600, state=r1["state"], groups=None, seed=seed + 1)
        r3 = fs.run(m, 200, state=r2["state"], groups=None, seed=seed + 2)
        out.update(act_400_1000=int((r2["counts"][0] > 0).sum()), act_1000_1200=int((r3["counts"][0] > 0).sum()), exploded=False)
    out["reverb"] = bool(out["exploded"] or out["act_1000_1200"] > 0)
    out["wall_s"] = round(time.time() - t0, 1)
    return out


def probe_task(ds, w, tag):
    m = mw(ds, w); P, ro, named, _ = BAT[ds]; t0 = time.time()
    res = dict(ds=ds, w=w, tag=tag, probes={})
    plist = list(P.items()) + [("sugar50", (P["sugar"][0], 50.0)), ("sugar100", (P["sugar"][0], 100.0)),
                               ("sugar200", (P["sugar"][0], 200.0))]
    for j, (name, (idx, rate)) in enumerate(plist):
        r = fs.run(m, 250, stim=[dict(idx=idx, rate=rate)], groups=None, seed=50 + j)
        st = r["state"]
        r2 = fs.run(m, 300, state=st, groups=None, seed=80 + j, epochs=[(200, 300)])
        hz = r["counts"][0] / 0.25
        d = dict(rate=rate, n_stim=int(len(idx)), n_active=int((r["counts"][0] > 0).sum()),
                 n_brain_ro_8Hz=int((hz[ro["brain"]] > 8).sum()), n_vnc_ro_8Hz=int((hz[ro["vnc"]] > 8).sum()) if len(ro["vnc"]) else 0,
                 post_active_200_300=int((r2["counts"][0] > 0).sum()),
                 named={k: round(float(hz[v].mean()), 1) if len(v) else None for k, v in named.items()},
                 region_Hz={rg: round(float(hz[(m.neurons.region.fillna("NA") == rg).values].mean()), 3)
                            for rg in ["OL", "CB", "VNC"]} if ds != "FAFB" else {})
        res["probes"][name] = d
    res["wall_s"] = round(time.time() - t0, 1)
    return res


def _worker(task):
    kind = task[0]
    try:
        return kind, (rev_task(*task[1:]) if kind == "rev" else probe_task(*task[1:]))
    except Exception:
        return kind, dict(task=list(map(str, task)), error=traceback.format_exc())


def pmap(tasks, workers):
    import multiprocessing as mp_
    if workers <= 1:
        return [_worker(t) for t in tasks]
    with mp_.get_context("fork").Pool(min(workers, len(tasks))) as pool:
        out = []
        for i, r in enumerate(pool.imap_unordered(_worker, tasks)):
            out.append(r)
            if r[0] == "rev" and "error" not in r[1]:
                x = r[1]; print(f"[{i+1}/{len(tasks)}] rev {x['ds']} {x['probe']} w={x['w']:.4f} s{x['seed']} act {x['act_0_150']}/"
                                f"{x['act_200_400']}/{x['act_1000_1200']} reverb={x['reverb']} {x['wall_s']}s", flush=True)
            else:
                print(f"[{i+1}/{len(tasks)}] {r[0]} {('ERR ' + r[1].get('error','')[-300:]) if 'error' in r[1] else r[1].get('tag')}", flush=True)
        return out


def crit(rows, ds, probes=None):
    rr = [r for r in rows if r["ds"] == ds and (probes is None or r["probe"] in probes)]
    ws = sorted({r["w"] for r in rr})
    rv = {w: any(r["reverb"] for r in rr if r["w"] == w) for w in ws}
    fail = [w for w in ws if rv[w]]
    first_fail = min(fail) if fail else None
    passes = [w for w in ws if (first_fail is None or w < first_fail)]
    return (max(passes) if passes else None), first_fail


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--ds", default="FAFB,MCNS,BANC"); a = ap.parse_args()
    dsl = a.ds.split(","); os.makedirs(OUT, exist_ok=True); t0 = time.time()
    setup(dsl)
    json.dump({ds: BAT[ds][3] for ds in dsl}, open(os.path.join(OUT, "battery_report.json"), "w"), indent=1, default=str)
    print("models", {ds: (MODELS[ds].N, MODELS[ds].info.get("n_edges")) for ds in dsl}, f"{time.time()-t0:.0f}s", flush=True)
    rows = []
    # stage 1
    tasks = []
    for ds in dsl:
        lo, hi = RANGES[ds]; w = lo
        while w <= hi * 1.0001:
            for s in SEEDS1:
                for pn in PROBES:
                    tasks.append(("rev", ds, round(w, 5), s, pn))
            w *= 1.15
    rows += [r for k, r in pmap(tasks, a.workers) if "error" not in r]
    json.dump(rows, open(os.path.join(OUT, "rev_stage1.json"), "w"))
    # stage 2
    tasks = []
    for ds in dsl:
        for crit_name, prs in (("sugar", ["sugar"]), ("any", None)):
            wp, wf = crit(rows, ds, prs)
            if wp is None or wf is None:
                print("WARN no bracket", ds, crit_name, wp, wf, flush=True); continue
            w = wp * 1.02
            while w < wf:
                for s in SEEDS:
                    for pn in (prs or PROBES):
                        t_ = ("rev", ds, round(w, 5), s, pn)
                        if t_ not in tasks:
                            tasks.append(t_)
                w *= 1.02
    if tasks:
        rows += [r for k, r in pmap(tasks, a.workers) if "error" not in r]
    json.dump(rows, open(os.path.join(OUT, "rev_all.json"), "w"))
    summ = {}
    for ds in dsl:
        summ[ds] = {}
        for crit_name, prs in (("sugar", ["sugar"]), ("any", None)):
            wp, wf = crit(rows, ds, prs)
            summ[ds][f"w_pass_{crit_name}"] = wp; summ[ds][f"w_crit_{crit_name}"] = wf
        summ[ds]["w_pass"] = summ[ds]["w_pass_any"]; summ[ds]["w_crit"] = summ[ds]["w_crit_any"]
    if "FAFB" in summ:
        for crit_name in ("sugar", "any"):
            wf = summ["FAFB"].get(f"w_crit_{crit_name}")
            if not wf:
                continue
            k = 0.185 / wf
            for ds in dsl:
                if summ[ds].get(f"w_crit_{crit_name}"):
                    summ[ds][f"ratio_fafb_{crit_name}"] = k
                    summ[ds][f"w_matched_{crit_name}"] = round(summ[ds][f"w_crit_{crit_name}"] * k, 4)
        for ds in dsl:
            summ[ds]["w_matched"] = summ[ds].get("w_matched_any")
    print("SUMMARY", json.dumps(summ), flush=True)
    json.dump(summ, open(os.path.join(OUT, "calib_summary.json"), "w"), indent=1)
    # stage 3: probes
    tasks = []
    for ds in dsl:
        if not summ[ds].get("w_crit"):
            continue
        for r in RATIOS:
            tasks.append(("probe", ds, round(r * summ[ds]["w_crit"], 5), f"r{r}"))
        for kk in ("w_matched_any", "w_matched_sugar", "w_pass_sugar"):
            if summ[ds].get(kk):
                tasks.append(("probe", ds, summ[ds][kk], kk))
        if summ[ds].get("w_pass"):
            tasks.append(("probe", ds, summ[ds]["w_pass"], "w_pass"))
    pr = [r for k, r in pmap(tasks, a.workers)]
    json.dump(pr, open(os.path.join(OUT, "probes.json"), "w"), default=str)
    print("DONE", round(time.time() - t0), flush=True)
