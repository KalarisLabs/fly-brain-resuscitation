"""Validation sweep of the calibrated princeton5 model (w_syn 0.185) in ONE process (memory-light).
Also measures post-stimulus reverberation (active neurons 200-400 ms after a 150-ms stimulus ends).
Writes results/validation_calibrated.csv"""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src",):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flysim as fs
import validate_lif as V

OUT = os.path.join(_P.OUT, "sim", "results")
m = fs.load_model("princeton5")          # calibrated w_syn
RO = V.readouts(m); SS = V.stim_sets(m)
rows = []; t0 = time.time()
jobs = [(e, f, t) for e in ["sugar", "LC4_bilateral", "LC9_bilateral"] for f in V.FREQS for t in range(3)]
jobs += [(e, 180, t) for e in ["left_LC9+DNa02", "right_LC9+DNa02"] for t in range(3)]
for i, (e, f, t) in enumerate(jobs):
    r = fs.run(m, 1000.0, stim=[dict(idx=SS[e], rate=float(f))], seed=1000 * t + int(f), groups=None)
    rt = fs.rates(r)
    row = dict(connectome="princeton5", w_syn=m.w_syn, experiment=e, freq_hz=f, trial=t, wall_s=r["wall_s"],
               n_active=int((rt > 0).sum()))
    for k, ix in RO.items():
        row[k] = float(rt[ix].mean())
    rows.append(row); del r
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "validation_calibrated.csv"), index=False)
    print(i + 1, len(jobs), round(time.time() - t0), flush=True)
# reverberation check (both w_syn values, both connectomes)
rv = []
for conn, w in [("princeton5", 0.185), ("princeton5", 0.275)]:
    mm = m if w == 0.185 else None
    for e in ["sugar", "LC4_bilateral", "LC9_bilateral"]:
        r = fs.run(m, 600.0, stim=[dict(idx=SS[e], rate=150.0, t_on=0, t_off=150)], groups=None, seed=2,
                   epochs=[(0, 150), (200, 400), (400, 600)], lif=dict(w_syn=w))
        c = r["counts"]; rv.append(dict(connectome=conn, w_syn=w, experiment=e, active_stim=int((c[0] > 0).sum()),
                                        active_200_400=int((c[1] > 0).sum()), active_400_600=int((c[2] > 0).sum())))
        del r
        pd.DataFrame(rv).to_csv(os.path.join(OUT, "reverberation_check.csv"), index=False)
print("done", round(time.time() - t0))
