"""Fine low-f grid for hub-targeted removal (transition occurs at f<0.05; main grid too coarse there).
python fine_degree.py DS -> tables/<DS>_curves_fine.csv (regime 'degree_finegrid', 2 orders)."""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src",):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import perclib as P
DS = sys.argv[1]
OUT = os.path.join(P.ROOT, 'tables', f'{DS}_curves_fine.csv')
FG = [0, .0025, .005, .0075, .01, .0125, .015, .02, .025, .03, .04, .05, .06, .075, .1, .125, .15, .2, .3]
g = P.Graph(DS, 5); g.setup_bits()
dfs = []
for o in range(2):
    rng = np.random.default_rng(2000 + o)
    order = P.removal_order(g, 'degree', rng)
    rows = P.run_curve(g, order, FG, do_scc=True)
    df = pd.DataFrame(rows); df.insert(0, 'core', ''); df.insert(0, 'order', o)
    df.insert(0, 'regime', 'degree_finegrid'); df.insert(0, 'graph', 'thr5'); df.insert(0, 'dataset', DS)
    dfs.append(df); print(DS, o, flush=True)
pd.concat(dfs).to_csv(OUT, index=False)
print('DONE fine', DS)
