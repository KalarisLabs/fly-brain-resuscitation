"""Metabolic / ionic layer unit tests (README section 7, tests a and b) on the toy connectome."""
import numpy as np

from conftest import STIM


def test_a_full_supply_baseline_is_stationary(toy):
    fs, m = toy
    r = fs.run(m, 300, metab=dict(speedup=50), groups=None, epochs=[(0, 300)])
    st = r["state"]
    assert r["counts"].sum() == 0
    assert np.abs(st.K - fs.METAB_DEFAULTS["K0"]).max() == 0
    np.testing.assert_allclose(st.a / st.a0, 1.0, atol=1e-6)
    assert st.alive.all()


def test_a_kspk0_is_bit_identical_to_pure_lif(toy):
    fs, m = toy
    lif = fs.run(m, 300, stim=STIM, groups=None, epochs=[(0, 300)], seed=3)["counts"]
    met = fs.run(m, 300, stim=STIM, groups=None, epochs=[(0, 300)], seed=3,
                 metab=dict(speedup=50, k_spk=0.0))["counts"]
    assert lif.sum() > 0
    assert np.array_equal(lif, met)


def test_b_anoxia_depolarizes_and_raises_extracellular_k(toy):
    fs, m = toy
    r = fs.run(m, 800, supply=0.0, metab=dict(speedup=50), groups=None, epochs=[(0, 800)])
    s = r["slow"]
    assert s["e_mean"][-1] < 0.2
    assert s["dep_frac"][-1] > 0.9
    assert np.max(s["K"]) > 30.0
    assert s["syn_frac"][-1] < 0.05
