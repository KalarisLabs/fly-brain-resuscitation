"""Smoke run of the full insult -> reperfusion -> probe pipeline on the 300-neuron toy connectome (seconds, no data).

python experiments/smoke/toy_protocol.py   -> prints a JSON summary; used by CI.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, "src"))
import numpy as np  # noqa: E402

import flysim_mech as fs  # noqa: E402
import toynet  # noqa: E402

m = toynet.make_toy_model(fs, w_syn=2.0)
probe = [dict(idx=toynet.SENSORY, rate=150, t_on=0, t_off=100)]
MET = dict(speedup=50)
ref = fs.run(m, 150, stim=probe, groups=None, epochs=[(0, 150)], seed=0, metab=MET)
r_ref = ref["counts"][0] / 0.15
s1 = fs.run(m, 400, supply=fs.piecewise([(0, 1.0), (2, 0.0)]), metab=MET, groups=None, seed=0)       # 2 s normoxia, anoxia
s2 = fs.run(m, 600, supply=1.0, metab=MET, groups=None, seed=0, state=s1["state"],                   # reperfusion + 20 s block
            interventions=dict(block_all=lambda t: t < s1["slow"]["t_s"][-1] + 20.0))
pr = fs.run(m, 150, stim=probe, groups=None, epochs=[(0, 150)], seed=0, metab=MET, state=s2["state"])
d = fs.functional_identity_index(r_ref, pr["counts"][0] / 0.15)
out = dict(anoxia_end_e=float(s1["slow"]["e_mean"][-1]), anoxia_peak_K=float(np.max(s1["slow"]["K"])),
           post_block_e=float(s2["slow"]["e_mean"][-1]), dead_frac=float(s2["slow"]["dead_frac"][-1]), FII=d["FII"])
print(json.dumps(out, indent=1))
assert np.isfinite(out["post_block_e"])
