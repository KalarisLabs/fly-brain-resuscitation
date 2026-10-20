"""Runs on Modal (Volume flysim-data at /data). For every dynamics run, read the unstimulated post-protocol
spike counts (cnt_i0/cnt_c0, 500-ms window) and attribute them to neuron groups (groups.npz).
Writes out/spont_attrib.csv (per run x group spikes and active neurons) and out/spont_mean.npz (mean per-neuron
counts for runaway vs non-runaway runs)."""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src",):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import numpy as np, glob, os, json, csv
z = np.load('groups.npz'); names = list(z['names']); N = int(z['N'])
G = np.unpackbits(z['masks'], axis=1)[:, :N].astype(bool)
os.makedirs('out', exist_ok=True)
rows = []; sum_run = np.zeros(N, np.float64); n_run = 0; sum_oth = np.zeros(N, np.float64); n_oth = 0
top_conc = []
for arm in ['A', 'B', 'C', 'D', 'E', 'H', 'base']:
    for f in sorted(glob.glob(_os.environ.get('RES_DIR', '/data/results') + f'/{arm}/*.npz')):
        rid = os.path.basename(f)[:-4]
        try:
            q = np.load(f); i = q['cnt_i0']; c = q['cnt_c0'].astype(np.float64)
        except Exception as e:
            print('ERR', rid, e); continue
        js = json.load(open(f[:-4] + '.json')) if os.path.exists(f[:-4] + '.json') else {}
        v = np.zeros(N); v[i] = c
        tot = v.sum(); r = dict(id=rid, total=tot, n_active=int((v > 0).sum()))
        for k, gm in zip(names, G):
            r['spk_' + k] = float(v[gm].sum()); r['act_' + k] = int((v[gm] > 0).sum())
        srt = np.sort(v)[::-1]
        for q_ in (10, 100, 1000, 10000):
            r[f'top{q_}_frac'] = float(srt[:q_].sum() / tot) if tot > 0 else np.nan
        rows.append(r)
        runaway = (js.get('obs_late_Hz') or 0) > 0.5 or (js.get('FSA') or 0) > 0.05
        if runaway: sum_run += v; n_run += 1
        else: sum_oth += v; n_oth += 1
keys = list(rows[0].keys())
with open('out/spont_attrib.csv', 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=keys); wr.writeheader(); wr.writerows(rows)
np.savez_compressed('out/spont_mean.npz', runaway=(sum_run / max(n_run, 1)).astype(np.float32), other=(sum_oth / max(n_oth, 1)).astype(np.float32), n_run=n_run, n_oth=n_oth)
print('runs', len(rows), 'runaway', n_run, 'other', n_oth)
