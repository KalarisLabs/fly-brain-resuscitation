"""Damage-percolation curves for one dataset.  python run_ds.py DS [part]
part = main (default) | controls
Appends one CSV row per (regime, order, f) to tables/<DS>_curves[_controls].csv; skips orders already done.
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
import sys, os, time, gc
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import perclib as P

DS = sys.argv[1]; PART = sys.argv[2] if len(sys.argv) > 2 else 'main'
FGRID_FULL = [0, .01, .02, .03, .05, .075, .1, .125, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9]
# reduced grid (subset of FGRID_FULL) adopted mid-run because of CPU contention; logged as deviation
FGRID = [0, .02, .05, .1, .15, .2, .3, .4, .5, .6, .7, .8, .9]
OUT = os.path.join(P.ROOT, 'tables', f'{DS}_curves.csv' if PART == 'main' else f'{DS}_curves_controls.csv')


def done_keys():
    if not os.path.exists(OUT): return set()
    d = pd.read_csv(OUT, usecols=['graph', 'regime', 'order'])
    return set(zip(d.graph, d.regime, d.order))


def run(g, graph_tag, regime, n_orders, demand=None, seed_base=0, fgrid=FGRID):
    dk = done_keys()
    for o in range(n_orders):
        if (graph_tag, regime, o) in dk: continue
        rng = np.random.default_rng(1000 * seed_base + o)
        t = time.time(); core = ''
        if regime == 'focal':
            order, core = P.neuropil_order(DS, g, rng)
        else:
            order = P.removal_order(g, regime, rng, demand)
        rows = P.run_curve(g, order, fgrid, do_scc=(o < 3))
        df = pd.DataFrame(rows); df.insert(0, 'core', core); df.insert(0, 'order', o)
        df.insert(0, 'regime', regime); df.insert(0, 'graph', graph_tag); df.insert(0, 'dataset', DS)
        df.to_csv(OUT, mode='a', header=not os.path.exists(OUT), index=False)
        print(DS, graph_tag, regime, o, core, f'{time.time() - t:.0f}s', flush=True)


# Reduced design (CPU contention from concurrent simulation; logged deviation). Orders already on disk are kept.
NR = {'FAFB': dict(random=13, degree=3, met=10, met2=5, focal=5, matched=3),
      'BANC': dict(random=12, degree=3, met=8, focal=4, matched=2),
      'MCNS': dict(random=6, degree=2, met=6, focal=3, matched=2),
      'MANC': dict(random=20, degree=5, met=20, focal=10, matched=5)}[DS]
if PART == 'main':
    g = P.Graph(DS, 5); g.setup_bits()
    nodes = g.nodes
    run(g, 'thr5', 'random', NR['random'], seed_base=1)
    run(g, 'thr5', 'degree', NR['degree'], seed_base=2)
    if DS == 'FAFB':
        run(g, 'thr5', 'metabolic_area', NR['met'], demand=nodes.area_nm.values, seed_base=3)
    run(g, 'thr5', 'metabolic_syn', NR.get('met2', NR['met']), demand=nodes.syn_total.values, seed_base=4)
    run(g, 'thr5', 'focal', NR['focal'], seed_base=5)
    run(g, 'thr5', 'degree_matched_random', NR['matched'], seed_base=6)
else:
    # H1 controls: one degree-preserving rewired null; edge-threshold sensitivity
    if DS in ('FAFB', 'MANC'):
        g = P.Graph(DS, 5, rewire_seed=100); g.setup_bits()
        run(g, 'null0', 'random', 3 if DS == 'FAFB' else 2, seed_base=11)
        run(g, 'null0', 'degree', 1, seed_base=21)
        del g; gc.collect()
    thrs = [10] + ([1] if DS == 'MANC' else [])
    for thr in thrs:
        g = P.Graph(DS, thr); g.setup_bits()
        run(g, f'thr{thr}', 'random', 3 if DS == 'FAFB' else 2, seed_base=31)
        run(g, f'thr{thr}', 'degree', 1, seed_base=32)
        del g; gc.collect()
print('DONE', DS, PART)
