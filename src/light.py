"""Light FAFB princeton5 loader (identical matrix and sign to flysim_dyn.load_model('princeton5')),
for remote runs without the full simulator. Requires research/cache/fafb783_syncount_pre_post_csr.npz and
followup/index/stage/fafb_meta.npz (sign, base_demand, fast-NT mask, types)."""
import os, numpy as np, pandas as pd, scipy.sparse as sp
from types import SimpleNamespace
D = os.environ.get('RESEARCH_ROOT', os.environ.get("RESEARCH_ROOT", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs")))
LIF_DEFAULTS = dict(v_rest=-52.0, v_reset=-52.0, v_th=-45.0, tau_m=20.0, tau_syn=5.0, t_ref=2.2)

def load_light():
    A = sp.load_npz(os.path.join(D, 'cache', 'fafb783_syncount_pre_post_csr.npz')).tocsr(); A.sort_indices()
    m = np.load(os.path.join(D, 'followup', 'index', 'stage', 'fafb_meta.npz'), allow_pickle=True)
    sign = m['sign']; w = A.data.astype(np.float32) * np.repeat(sign, np.diff(A.indptr)).astype(np.float32)
    nt = pd.DataFrame(dict(root_id=m['root_id'], primary_type=m['primary_type'], super_class=m['super_class'], nt=m['nt']))
    return SimpleNamespace(N=A.shape[0], indptr=A.indptr.astype(np.int64), indices=A.indices.astype(np.int32), weight=w,
                           sign=sign, neurons=nt, base_demand=m['base_demand'], fast=m['fast'].astype(bool))
