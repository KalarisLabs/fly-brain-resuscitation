"""Path resolution and the numpy<2 popcount fallback in perclib."""
import importlib
import os

import numpy as np


def test_paths_follow_environment(monkeypatch, tmp_path):
    import fbr_paths
    monkeypatch.setenv("RESEARCH_ROOT", str(tmp_path / "out"))
    monkeypatch.setenv("RAW_ROOT", str(tmp_path / "raw"))
    monkeypatch.delenv("FLYSIM_CACHE", raising=False)
    p = importlib.reload(fbr_paths)
    assert p.OUT == str(tmp_path / "out") and p.RAW == str(tmp_path / "raw")
    assert p.CACHE == os.path.join(str(tmp_path / "out"), "cache")
    monkeypatch.undo()
    importlib.reload(fbr_paths)


def test_popcount_fallback_matches_reference():
    import perclib
    rng = np.random.default_rng(1)
    x = rng.integers(0, 2**63 - 1, size=(17, 3), dtype=np.int64).astype(np.uint64)
    ref = np.vectorize(lambda v: bin(int(v)).count("1"))(x)
    assert np.array_equal(perclib.popcount64(x), ref)
    fb = perclib._T8[np.ascontiguousarray(x).view(np.uint8)].reshape(x.shape + (8,)).sum(-1)
    assert np.array_equal(fb, ref)
