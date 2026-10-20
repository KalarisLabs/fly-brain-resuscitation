"""Validation sweeps of the Shiu LIF engine (flysim) on FAFB v783.
Run:  python validate_lif.py [n_proc] [n_trials]
Writes $RESEARCH_ROOT\\sim\\results\\validation_lif.csv (+ mean-rate vectors .npz)."""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src",):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import os, sys, time, itertools
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flysim as fs

OUT = os.path.join(_P.OUT, "sim", "results")
FREQS = [10, 25, 50, 75, 100, 125, 150, 175, 200]
_M = {}

def model(name):
    if name not in _M:
        _M[name] = fs.load_model(name, w_syn=float(os.environ.get("FLYSIM_WSYN", 0.275)))
    return _M[name]

def readouts(m):
    t = m.neurons
    R = {"MN9": m.idx([fs.SHIU_MN9])}
    for typ in ["DNp01", "DNa02", "DNa01", "DNp09", "DNg100", "DNg97", "DNa11", "DNg13", "MDN", "DNp02", "DNp11"]:
        for side in ["left", "right"]:
            ix = m.select(typ, side=side)
            if len(ix): R[f"{typ}_{side[0].upper()}"] = ix
    return R

def stim_sets(m):
    return {
        "sugar": m.idx(fs.SHIU_SUGAR_GRN),
        "LC4_bilateral": m.select("LC4"),
        "LC9_bilateral": m.select("LC9"),                                         # Xenova 'walk'
        "left_LC9+DNa02": np.concatenate([m.select("LC9", side="left"), m.select("DNa02", side="left")]),
        "right_LC9+DNa02": np.concatenate([m.select("LC9", side="right"), m.select("DNa02", side="right")]),
    }

def job(args):
    conn, exp, freq, trial = args
    m = model(conn)
    S = stim_sets(m)[exp]
    r = fs.run(m, 1000.0, stim=[dict(idx=S, rate=float(freq))], seed=1000 * trial + int(freq), groups=None)
    rt = fs.rates(r)
    row = dict(connectome=conn, experiment=exp, freq_hz=freq, trial=trial, wall_s=r["wall_s"],
               n_stim=len(S), n_active=int((rt > 0).sum()), total_spikes=int(r["counts"].sum()))
    for k, ix in readouts(m).items():
        row[k] = float(rt[ix].mean())
    return row, (conn, exp, freq), rt.astype(np.float32)

def all_jobs(ntr):
    jobs = [("princeton5", e, f, t) for e in ["sugar", "LC4_bilateral", "LC9_bilateral"] for f in FREQS for t in range(ntr)]
    jobs += [("princeton5", e, f, t) for e in ["left_LC9+DNa02", "right_LC9+DNa02"] for f in [50, 100, 180] for t in range(ntr)]
    jobs += [("shiu783", "sugar", f, t) for f in FREQS for t in range(ntr)]
    return jobs

def run_shard(i, n, ntr):
    jobs = all_jobs(ntr)[i::n]
    rows, acc = [], {}
    t0 = time.time()
    for j, jb in enumerate(jobs):
        row, key, rt = job(jb)
        rows.append(row); acc.setdefault(key, []).append(rt)
        pd.DataFrame(rows).to_csv(os.path.join(OUT, f"val_shard{i}.csv"), index=False)
        print(f"shard {i}: {j+1}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    np.savez_compressed(os.path.join(OUT, f"val_shard{i}_rates.npz"),
                        **{f"{c}|{e}|{f}|n{len(v)}": np.sum(v, 0) for (c, e, f), v in acc.items()})

def merge(n):
    df = pd.concat([pd.read_csv(os.path.join(OUT, f"val_shard{i}.csv")) for i in range(n)])
    df.to_csv(os.path.join(OUT, "validation_lif.csv"), index=False)
    S, C = {}, {}
    for i in range(n):
        z = np.load(os.path.join(OUT, f"val_shard{i}_rates.npz"))
        for k in z.files:
            c, e, f, nn = k.split("|"); key = f"{c}|{e}|{f}"
            S[key] = S.get(key, 0) + z[k]; C[key] = C.get(key, 0) + int(nn[1:])
    np.savez_compressed(os.path.join(OUT, "validation_mean_rates.npz"), **{k: (S[k] / C[k]).astype(np.float32) for k in S})
    return df

if __name__ == "__main__":
    # python validate_lif.py launch NPROC NTRIALS  |  shard I N NTRIALS  |  merge N
    import subprocess
    os.makedirs(OUT, exist_ok=True)
    cmd = sys.argv[1]
    if cmd == "shard":
        run_shard(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
    elif cmd == "merge":
        merge(int(sys.argv[2]))
    else:
        n, ntr = int(sys.argv[2]), int(sys.argv[3])
        procs = []
        for i in range(n):
            lf = open(os.path.join(OUT, f"val_shard{i}.log"), "w")
            procs.append(subprocess.Popen([sys.executable, os.path.abspath(__file__), "shard", str(i), str(n), str(ntr)],
                                          stdout=lf, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL))
        for pr in procs:
            pr.wait()
        merge(n)
        print("done")
