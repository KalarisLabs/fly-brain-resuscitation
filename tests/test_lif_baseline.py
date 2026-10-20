"""Fast (Shiu LIF) layer on a 300-neuron toy connectome: stationarity, propagation, determinism, engine identity."""
import numpy as np

from conftest import ENGINES, STIM
import toynet


def test_no_input_is_exactly_stationary(toy):
    fs, m = toy
    r = fs.run(m, 200, groups=None, epochs=[(0, 200)])
    assert r["counts"].sum() == 0
    assert np.all(r["state"].v == np.float32(fs.LIF_DEFAULTS["v_rest"]))
    assert np.all(r["state"].g == 0)


def test_stimulus_propagates_and_network_returns_to_rest(toy):
    fs, m = toy
    r = fs.run(m, 400, stim=STIM, groups=None, epochs=[(0, 100), (150, 250), (300, 400)], seed=1)
    c = r["counts"]
    assert c[0].sum() > 0
    assert (c[0][len(toynet.SENSORY):] > 0).sum() >= 5, "activity should spread beyond the driven neurons"
    assert c[1].sum() == 0 and c[2].sum() == 0, "toy network is subcritical: silent after the stimulus"


def test_same_seed_is_deterministic(toy):
    fs, m = toy
    a = fs.run(m, 150, stim=STIM, groups=None, epochs=[(0, 150)], seed=7)["counts"]
    b = fs.run(m, 150, stim=STIM, groups=None, epochs=[(0, 150)], seed=7)["counts"]
    assert np.array_equal(a, b)


def test_engine_versions_identical_with_default_options():
    """flysim -> flysim_dyn -> flysim_mech only add options whose defaults leave the numerics unchanged."""
    out = {}
    for name, fs in ENGINES.items():
        m = toynet.make_toy_model(fs, w_syn=2.0)
        out[name] = fs.run(m, 200, stim=STIM, groups=None, epochs=[(0, 200)], seed=3,
                           metab=dict(speedup=50))["counts"]
    assert np.array_equal(out["flysim"], out["flysim_dyn"])
    assert np.array_equal(out["flysim_dyn"], out["flysim_mech"])
