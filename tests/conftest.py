import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import flysim  # noqa: E402
import flysim_dyn  # noqa: E402
import flysim_mech  # noqa: E402
import toynet  # noqa: E402

ENGINES = {"flysim": flysim, "flysim_dyn": flysim_dyn, "flysim_mech": flysim_mech}
STIM = [dict(idx=toynet.SENSORY, rate=150, t_on=0, t_off=100)]


@pytest.fixture(params=sorted(ENGINES))
def engine(request):
    return ENGINES[request.param]


@pytest.fixture
def toy(engine):
    return engine, toynet.make_toy_model(engine, w_syn=2.0)
