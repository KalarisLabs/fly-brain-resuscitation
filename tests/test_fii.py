"""Functional Identity Index: FII = max(r, 0) * sqrt(recall * precision)."""
import math

import numpy as np

from conftest import ENGINES

fii = ENGINES["flysim_mech"].functional_identity_index


def test_identical_profiles_give_one():
    x = np.array([0, 0, 5, 10, 20, 40, 80.0])
    d = fii(x, x)
    assert d["FII"] == 1.0 and d["recall"] == 1.0 and d["precision"] == 1.0


def test_disjoint_responders_give_zero():
    a = np.array([10, 20, 30, 0, 0, 0.0])
    b = np.array([0, 0, 0, 10, 20, 30.0])
    assert fii(a, b)["FII"] == 0.0


def test_novel_responders_reduce_precision():
    a = np.array([10, 20, 30, 40, 0, 0.0])
    b = np.array([10, 20, 30, 40, 25, 25.0])
    d = fii(a, b)
    assert d["recall"] == 1.0 and math.isclose(d["precision"], 4 / 6)
    assert 0 < d["FII"] < 1


def test_too_few_responsive_neurons_is_nan():
    assert math.isnan(fii(np.array([5.0, 0, 0]), np.array([5.0, 0, 0]))["FII"])


def test_fii_bounded_and_identical_across_engines():
    rng = np.random.default_rng(0)
    for _ in range(50):
        a, b = rng.exponential(5, 200) * (rng.random(200) < 0.3), rng.exponential(5, 200) * (rng.random(200) < 0.3)
        vals = [e.functional_identity_index(a, b)["FII"] for e in ENGINES.values()]
        assert 0.0 <= vals[0] <= 1.0
        assert vals[0] == vals[1] == vals[2]
