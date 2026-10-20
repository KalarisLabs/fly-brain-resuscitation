"""Metabolic-layer unit tests + demonstration runs (writes results/metab_*.npz, results/metab_tests.json).
  test_a : S=1, k_spk=0 -> spike counts identical to pure Shiu LIF (same seed); default k_spk -> near-identical
  test_b : S=0 -> SD-like depolarization (dep_frac>0.9, K>30 mM) and population silence
  test_c : restoring S after brief anoxia (small J) -> energy, K, firing and sugar->MN9 response recover
  focal  : ischemic core (right optic lobe neuropils S=0) -> K+ spread to neighbouring compartments
python metab_tests.py [which]   which in {a,bc,focal,all}"""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src",):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flysim as fs

OUT = os.path.join(_P.OUT, "sim", "results")
os.makedirs(OUT, exist_ok=True)
which = sys.argv[1] if len(sys.argv) > 1 else "all"
m = fs.load_model("princeton5")
sug = m.idx(fs.SHIU_SUGAR_GRN); mn9 = int(m.idx([fs.SHIU_MN9])[0])
res_path = os.path.join(OUT, "metab_tests.json")
R = json.load(open(res_path)) if os.path.exists(res_path) else {}

def save():
    cur = json.load(open(res_path)) if os.path.exists(res_path) else {}
    cur.update(R); json.dump(cur, open(res_path, "w"), indent=1)

if which in ("a", "all"):
    ra = fs.run(m, 1000, stim=[dict(idx=sug, rate=150)], seed=3, groups=None)
    rb = fs.run(m, 1000, stim=[dict(idx=sug, rate=150)], seed=3, groups=None, metab=dict(k_spk=0.0), supply=1.0)
    rc = fs.run(m, 1000, stim=[dict(idx=sug, rate=150)], seed=3, groups=None, metab={}, supply=1.0)
    x0, x1, x2 = fs.rates(ra), fs.rates(rb), fs.rates(rc)
    r0 = fs.run(m, 500, metab={}, supply=1.0, groups=None)
    R["test_a"] = dict(identical_kspk0=bool(np.array_equal(ra["counts"], rb["counts"])),
                       default_n_neurons_differing=int((x0 != x2).sum()), default_max_abs_dHz=float(np.abs(x0 - x2).max()),
                       default_rate_corr=float(np.corrcoef(x0, x2)[0, 1]), MN9_lif=float(x0[mn9]), MN9_metab_default=float(x2[mn9]),
                       baseline_nostim_spikes=int(r0["counts"].sum()), baseline_K_dev=float(np.abs(r0["slow"]["K"] - 5).max()),
                       baseline_e_min=float(r0["e"].min()),
                       wall_lif=ra["wall_s"], wall_metab=rc["wall_s"], wall_metab_nostim_per_s=r0["wall_s"] / 0.5,
                       passed=bool(np.array_equal(ra["counts"], rb["counts"]) and r0["counts"].sum() == 0))
    save(); print("test_a", R["test_a"], flush=True)

if which in ("bc", "all"):
    # speedup 50: 1 ms of spiking time = 50 ms of metabolic time. Slow-time protocol (s):
    #   0-10 normoxia (sugar probe 1-8.5 s) | 10-40 anoxia S=0 | 40-100 reoxygenation (sugar probe 92.5-100 s)
    sp = 50.0
    sup = fs.piecewise([(0, 1.0), (10.0, 0.0), (40.0, 1.0)])
    probe1 = (20.0, 170.0); probe2 = (1850.0, 2000.0)        # bio ms
    rng = np.random.default_rng(0); rec = np.sort(rng.choice(m.N, 4000, replace=False))
    rec = np.union1d(rec, np.r_[sug, mn9])
    t0 = time.time()
    r = fs.run(m, 2000.0, stim=[dict(idx=sug, rate=150, t_on=probe1[0], t_off=probe1[1]),
                                dict(idx=sug, rate=150, t_on=probe2[0], t_off=probe2[1])],
               metab=dict(speedup=sp), supply=sup, epochs=[probe1, (200, 800), (800, 1850), probe2],
               record_idx=rec, max_raster=3_000_000, seed=11, rec_every_ms=2.0, verbose=True)
    sl = r["slow"]; t = sl["t_s"]
    pre = fs.rates(r, epoch=0); post = fs.rates(r, epoch=3)
    fii = fs.functional_identity_index(pre, post)
    def at(ts, key):
        return float(sl[key][min(np.searchsorted(t, ts), len(t) - 1)])
    K = sl["K"]; Km = lambda ts: float(K[min(np.searchsorted(t, ts), len(t) - 1)].mean())
    dep = np.array(sl["dep_frac"])
    onset = float(t[np.argmax(dep > 0.5)] - 10.0) if (dep > 0.5).any() else None
    rec_t = t[(t > 40) & (dep < 0.01)]
    R["test_bc"] = dict(speedup=sp, wall_s=r["wall_s"], SD_onset_s_after_anoxia=onset,
                        dep_frac_at_35s=at(35, "dep_frac"), K_mean_at_35s=Km(35), block_frac_at_35s=at(35, "block_frac"),
                        e_mean_at_39s=at(39, "e_mean"), syn_frac_at_39s=at(39, "syn_frac"),
                        pop_rate_late_anoxia_Hz_per_neuron=float(r["pop_ms"].sum(1)[600:800].sum() / 0.2 / m.N),
                        repolarized_s_after_reox=(float(rec_t[0] - 40) if rec_t.size else None),
                        e_mean_end=at(99, "e_mean"), dep_frac_end=at(99, "dep_frac"), K_mean_end=Km(99),
                        dead_frac_end=at(99, "dead_frac"), J_mean_at_40s=at(40, "J_mean"),
                        MN9_pre_Hz=float(pre[mn9]), MN9_post_Hz=float(post[mn9]), FII=fii,
                        passed_b=bool(at(35, "dep_frac") > 0.9 and Km(35) > 30),
                        passed_c=bool(at(99, "e_mean") > 0.9 and at(99, "dep_frac") < 0.01 and post[mn9] > 0.5 * pre[mn9]))
    np.savez_compressed(os.path.join(OUT, "metab_anoxia_reox.npz"), t_s=t, K=K, e_np=sl["e_np"], dep_np=sl["dep_frac_np"],
                        **{k: np.asarray(sl[k]) for k in ["e_mean", "e_lt05", "dep_frac", "block_frac", "dead_frac", "vr_mean", "syn_frac", "J_mean", "S_mean"]},
                        pop_ms=r["pop_ms"], group_names=np.array(r["group_names"]), raster_t_ms=r["raster_t_ms"], raster_idx=r["raster_idx"],
                        rec_idx=rec, counts=r["counts"], speedup=sp, neuropils=np.array(m.neuropils))
    save(); print("test_bc", R["test_bc"], time.time() - t0, flush=True)


def protocol_c(name, anoxia_end_s, v_off=None, block=None, metab_extra=None):
    sp = 50.0
    sup = fs.piecewise([(0, 1.0), (10.0, 0.0), (anoxia_end_s, 1.0)])
    probe1 = (20.0, 170.0); probe2 = (1850.0, 2000.0)
    iv = {} if v_off is None else {"v_offset": fs.piecewise(v_off)}
    if block is not None:
        iv["block_all"] = fs.piecewise(block)
    mx = dict(speedup=sp); mx.update(metab_extra or {})
    t0 = time.time()
    r = fs.run(m, 2000.0, stim=[dict(idx=sug, rate=150, t_on=probe1[0], t_off=probe1[1]),
                                dict(idx=sug, rate=150, t_on=probe2[0], t_off=probe2[1])],
               metab=mx, supply=sup, interventions=iv, epochs=[probe1, probe2],
               seed=11, rec_every_ms=2.0, groups=None)
    sl = r["slow"]; t = sl["t_s"]
    pre = fs.rates(r, epoch=0); post = fs.rates(r, epoch=1)
    at = lambda ts, key: float(sl[key][min(np.searchsorted(t, ts), len(t) - 1)])
    dep = np.array(sl["dep_frac"])
    pr = r["pop_ms"].sum(1)
    out = dict(anoxia_s=anoxia_end_s - 10.0, v_offset=v_off, block_all=block, metab_extra=metab_extra, wall_s=r["wall_s"], max_dep_frac=float(dep.max()),
               max_K=float(sl["K"].max()), e_min=float(np.min(sl["e_mean"])), J_mean_max=float(np.max(sl["J_mean"])),
               e_mean_end=at(99, "e_mean"), dep_frac_end=at(99, "dep_frac"), K_mean_end=float(sl["K"][-1].mean()),
               dead_frac_end=at(99, "dead_frac"), pop_Hz_per_neuron_80_92s=float(pr[1600:1850].sum() / 0.25 / m.N),
               MN9_pre_Hz=float(pre[mn9]), MN9_post_Hz=float(post[mn9]), FII=fs.functional_identity_index(pre, post),
               passed=bool(at(99, "e_mean") > 0.9 and at(99, "dep_frac") < 0.01 and post[mn9] > 0.5 * pre[mn9]
                           and fs.functional_identity_index(pre, post)["FII"] > 0.5))
    R[name] = out; save(); print(name, out, time.time() - t0, flush=True)
    np.savez_compressed(os.path.join(OUT, f"metab_{name}.npz"), t_s=t, K=sl["K"], **{k: np.asarray(sl[k]) for k in ["e_mean", "dep_frac", "block_frac", "dead_frac", "vr_mean", "syn_frac", "J_mean"]}, pop_ms=r["pop_ms"], pre=pre, post=post)

if which in ("c2", "all"):
    protocol_c("test_c_short_anoxia", 17.0)                                   # 7 s anoxia (< SD onset)
    protocol_c("test_c_anesthesia", 40.0, v_off=[(0, 0.0), (40.0, -10.0), (60.0, 0.0)])   # 30 s anoxia + -10 mV 40-60 s

if which in ("c3", "all"):
    protocol_c("test_c_reset", 40.0, block=[(0, False), (40.0, True), (60.0, False)])            # 30 s anoxia, total block 40-60 s
    protocol_c("test_c_short_reset", 17.0, block=[(0, False), (17.0, True), (37.0, False)])      # 7 s anoxia, total block 17-37 s

if which in ("focal", "all"):
    sp = 20.0
    core = [i for i, n in enumerate(m.neuropils) if n in ("ME_R", "LO_R", "LOP_R", "AME_R", "LA_R")]
    S = np.ones(len(m.neuropils)); S[core] = 0.0
    sup = fs.piecewise([(0, np.ones(len(m.neuropils))), (2.0, S)])
    r = fs.run(m, 3000.0, metab=dict(speedup=sp), supply=sup, seed=5, groups="super_class", rec_every_ms=5.0)
    sl = r["slow"]
    K = sl["K"]; t = sl["t_s"]
    peak = K.max(0)
    R["focal"] = dict(core=[m.neuropils[i] for i in core], wall_s=r["wall_s"],
                      K_peak_core=float(peak[core].mean()), n_noncore_K_gt_20=int((np.delete(peak, core) > 20).sum()),
                      noncore_top=[(m.neuropils[i], float(peak[i])) for i in np.argsort(-peak) if i not in core][:8],
                      dead_frac_end=float(sl["dead_frac"][-1]), adjacency=m.adj_source)
    np.savez_compressed(os.path.join(OUT, "metab_focal.npz"), t_s=t, K=K, dep_np=sl["dep_frac_np"], e_np=sl["e_np"],
                        core=np.array(core), neuropils=np.array(m.neuropils), pop_ms=r["pop_ms"], group_names=np.array(r["group_names"]))
    save(); print("focal", R["focal"], flush=True)
