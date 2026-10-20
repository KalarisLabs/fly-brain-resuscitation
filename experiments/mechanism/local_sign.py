"""Local (no Modal) reduced sign-robustness arm.  Independent processes (sandbox blocks mp pipes).
python local_sign.py PHASE K N     PHASE 1 = S/SN/SY conditions, 2 = ST targeted (needs drivers_A/B), 'drivers' = build A/B rankings
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
import os, sys, time
R0 = os.path.join(_P.OUT, 'followup', 'mechanism', 'raw')
os.environ.setdefault("RES_DIR", os.path.join(R0, 'sign_res'))
os.environ.setdefault("OUT_DIR", os.path.join(R0, 'sign_out'))
os.environ.setdefault("BASE_RES", os.path.join(_P.OUT, 'dynamics', 'raw'))
os.environ.setdefault("FLYSIM_CACHE", _P.CACHE)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mech_sweep as MS

S3 = (0, 1, 2); S10 = (10, 11, 12)


def phase1_ids(Ds=(30,), extra=True):
    ids = []
    for v in ("A", "B"):
        if extra:
            ids += [f"S_{v}_base_s{s}" for s in (100, 101, 102)]
        for D in Ds:
            for rp in ("naive", "reset"):
                ids += [f"S_{v}_D{D}_{rp}_s{s}" for s in S3]
            ids += [f"SN_{v}_anat_D{D}_s{s}" for s in S3]
        if extra:
            for K in ("0", "9", "15"):
                ids += [f"SY_{v}_kick_v0_K{K}_s{s}" for s in S3]
    return ids


def phase2_ids():
    ids = []
    for v in ("A", "B"):
        for p in ("naive", "reset", "tblk_sustain_0.5", "tblk_sustain_2", "tblk_hub_2"):
            ids += [f"ST_{v}_D30_{p}_s{s}" for s in S10]
    return ids


if __name__ == "__main__":
    ph = sys.argv[1]
    os.makedirs(os.environ["RES_DIR"], exist_ok=True); os.makedirs(os.environ["OUT_DIR"], exist_ok=True)
    t0 = time.time(); MS.setup()
    print("setup", round(time.time() - t0, 1), "s; variants changed:", MS.MV_N, flush=True)
    if ph == "drivers":
        for v in ("A", "B"):
            MS.build_drivers(v); print("drivers", v, "built", flush=True)
        sys.exit()
    k, n = int(sys.argv[2]), int(sys.argv[3])
    G = {c["id"]: c for c in MS.grid()}
    ids = {"1": phase1_ids(), "1b": phase1_ids((180,), extra=False), "2": phase2_ids()}[ph]
    miss = [i for i in ids if i not in G]
    assert not miss, miss
    todo = sorted(ids, key=lambda i: -MS.cost(G[i]))[k::n]
    for i in todo:
        r = MS.worker(G[i]); print(f"[{time.time()-t0:7.0f}s]", r, flush=True)
    print("DONE", round(time.time() - t0), flush=True)
