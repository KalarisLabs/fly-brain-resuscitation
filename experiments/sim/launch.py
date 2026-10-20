"""Run several scripts as independent processes (the sandbox blocks multiprocessing pipes).
python launch.py "metab_tests.py a" "metab_tests.py bc" ...   -> logs in results/<name>.log"""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src",):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import os, sys, subprocess, time
here = os.path.dirname(os.path.abspath(__file__)); os.makedirs(os.path.join(_P.OUT, "sim", "results"), exist_ok=True)
ps = []
for spec in sys.argv[1:]:
    args = spec.split()
    name = "_".join(a.replace(".py", "") for a in args)
    lf = open(os.path.join(_P.OUT, "sim", "results", f"{name}.log"), "w")
    ps.append((name, subprocess.Popen([sys.executable, "-u"] + [os.path.join(here, args[0])] + args[1:], cwd=here,
                                      stdout=lf, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)))
t0 = time.time()
for name, p in ps:
    p.wait(); print(name, "exit", p.returncode, round(time.time() - t0), "s", flush=True)
