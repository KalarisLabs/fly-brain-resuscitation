"""Reproduce the repository figures from the committed summary tables in results/ (no FlyWire data needed).

python scripts/reproduce_figures.py            # all: stage tables, structure figures, table-derived re-plots
python scripts/reproduce_figures.py structure  # Fig. structure 1-4 with the original plotting code
python scripts/reproduce_figures.py replots    # dynamics / mechanism / replication re-plots from tables

Outputs go to $RESEARCH_ROOT (default <repo>/outputs), mirroring the original research layout:
  outputs/structure/figures/fig{1..4}_*.png|pdf      (experiments/structure/make_figures.py, original code)
  outputs/replots/*.png                               (scripts/replot_from_tables.py, written for this repo)

The original plotting code for the dynamics (fig1-5), mechanism (fig_M1-M4) and simulator-validation figures in
docs/figures/ was run interactively and was not preserved; replot_from_tables.py re-draws the same quantities
from the committed tables so that every headline number can be checked visually.  See docs/REPRODUCE.md.
"""
import os
import shutil
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
sys.path.insert(0, os.path.join(ROOT, "src"))
import fbr_paths as P  # noqa: E402

STAGE = {  # results/<track>  ->  $RESEARCH_ROOT/<layout used by the experiment scripts>
    "structure": os.path.join("structure", "tables"),
    "dynamics": os.path.join("dynamics", "tables"),
    "mechanism": os.path.join("followup", "mechanism", "tables"),
    "replication": os.path.join("followup", "replication", "tables"),
    "index": os.path.join("followup", "index", "tables"),
    "sim": os.path.join("sim", "results"),
}


def stage():
    for src, dst in STAGE.items():
        s, d = os.path.join(ROOT, "results", src), os.path.join(P.OUT, dst)
        os.makedirs(d, exist_ok=True)
        for f in os.listdir(s):
            if os.path.isfile(os.path.join(s, f)):
                shutil.copy2(os.path.join(s, f), os.path.join(d, f))
    print("staged results/ ->", P.OUT)


def run(script):
    env = dict(os.environ, RESEARCH_ROOT=P.OUT, MPLBACKEND="Agg")
    print("+", script, flush=True)
    subprocess.run([sys.executable, os.path.join(ROOT, script)], check=True, cwd=ROOT, env=env)


if __name__ == "__main__":
    what = set(sys.argv[1:]) or {"structure", "replots"}
    stage()
    if "structure" in what:
        os.makedirs(os.path.join(P.OUT, "structure", "figures"), exist_ok=True)
        run(os.path.join("experiments", "structure", "make_figures.py"))
    if "replots" in what:
        run(os.path.join("scripts", "replot_from_tables.py"))
