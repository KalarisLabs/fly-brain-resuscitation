"""Keystone cell types: drop in strict sensory->sink reach when all removable neurons of one type are removed.
python keystone.py DS [extra_types_file]
Candidates = top-300 types by flow centrality (sum over member neurons of fwd walks<=3 from sources x bwd walks<=3 to sinks,
strict relay) + any names listed in extra_types_file (for cross-dataset replication pass).
Output: tables/<DS>_keystone.csv (appends; resumable).
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
import sys, os, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import perclib as P

DS = sys.argv[1]
extra = []
if len(sys.argv) > 2 and os.path.exists(sys.argv[2]):
    extra = [l.strip() for l in open(sys.argv[2], encoding='utf8') if l.strip()]
OUT = os.path.join(P.ROOT, 'tables', f'{DS}_keystone.csv')
g = P.Graph(DS, 5); g.setup_bits()
nodes = g.nodes
ty = nodes['type'].astype(str).where(nodes['type'].notna(), '')
rem_mask = np.zeros(g.n, bool); rem_mask[g.removable] = True

# flow centrality (strict relay)
relay = ~(g.is_src | g.is_sink)
A = g.AT.T.tocsr()
x = np.zeros(g.n); x[g.src] = 1; fwd = np.zeros(g.n)
for k in range(3):
    x = g.AT @ x; fwd += x * relay; x = x * relay
y = np.zeros(g.n); y[g.sink] = 1; bwd = np.zeros(g.n)
for k in range(3):
    y = A @ y; bwd += y * relay; y = y * relay
flow = fwd * bwd
df = pd.DataFrame({'type': ty, 'flow': flow, 'rem': rem_mask})
df = df[df.rem & (df.type != '') & (df.type != 'nan')]
tstat = df.groupby('type').agg(n=('flow', 'size'), flow=('flow', 'sum')).sort_values('flow', ascending=False)
cands = list(tstat.index[:150]) + [t for t in extra if t in tstat.index and t not in tstat.index[:150]]
tstat.to_csv(os.path.join(P.ROOT, 'tables', f'{DS}_type_flow_centrality.csv'))

alive0 = np.ones(g.n, bool); alive0[g.excluded] = False

# ---- H1 negative control (c): remove every removable neuron that lies on no strict source->sink path ----
fwd_any = np.zeros(g.n, bool); fr = np.zeros(g.n); fr[g.src] = 1
for k in range(200):
    nx = (g.AT @ fr) > 0
    new = nx & ~fwd_any
    if not new.any(): break
    fwd_any |= new; fr = (new & relay).astype(float)
bwd_any = np.zeros(g.n, bool); br = np.zeros(g.n); br[g.sink] = 1
for k in range(200):
    nx = (A @ br) > 0
    new = nx & ~bwd_any
    if not new.any(): break
    bwd_any |= new; br = (new & relay).astype(float)
onpath = fwd_any & bwd_any
off = g.removable[~onpath[g.removable]]
a = alive0.copy(); a[off] = False
mneg = g.evaluate(a, do_inf=True, do_scc=False)
basefull = g.evaluate(alive0, do_inf=True, do_scc=False)
negrow = {'dataset': DS, 'n_removable': int(g.removable.size), 'n_offpath_removed': int(off.size),
          'f_offpath': off.size / g.removable.size}
for k in ['reachN_L6', 'reachN_Linf', 'reachG_L6', 'redund_L6', 'influence_L6']:
    negrow['base_' + k] = basefull[k]; negrow['neg_' + k] = mneg[k]; negrow['delta_' + k] = basefull[k] - mneg[k]
pd.DataFrame([negrow]).to_csv(os.path.join(P.ROOT, 'tables', f'{DS}_negative_control.csv'), index=False)
print('negctl', negrow, flush=True)
base = g.evaluate(alive0, do_inf=True, do_scc=False, maxL=30)
done = set(pd.read_csv(OUT).type) if os.path.exists(OUT) else set()
buf = []
t0 = time.time()
for i, t in enumerate(cands):
    if t in done: continue
    a = alive0.copy(); members = np.where((ty.values == t) & rem_mask)[0]; a[members] = False
    m = g.evaluate(a, do_inf=True, do_scc=False, maxL=30)
    row = {'dataset': DS, 'type': t, 'n_neurons': members.size, 'flow_rank': int(tstat.index.get_loc(t)) + 1,
           'from_extra': t not in tstat.index[:150]}
    for k in ['reachN_L4', 'reachN_L6', 'reachN_Linf', 'reachG_L6', 'reachN_L6_DN', 'reachN_L6_MN', 'redund_L6', 'influence_L6']:
        if k in m:
            row['d_' + k] = base[k] - m[k]
            row['rel_' + k] = (base[k] - m[k]) / base[k] if base[k] else np.nan
    buf.append(row)
    if len(buf) >= 20:
        pd.DataFrame(buf).to_csv(OUT, mode='a', header=not os.path.exists(OUT), index=False); buf = []
        print(DS, i + 1, len(cands), f'{time.time() - t0:.0f}s', flush=True)
if buf:
    pd.DataFrame(buf).to_csv(OUT, mode='a', header=not os.path.exists(OUT), index=False)
pd.Series(base).to_csv(os.path.join(P.ROOT, 'tables', f'{DS}_baseline_metrics.csv'))
print('DONE keystone', DS)
