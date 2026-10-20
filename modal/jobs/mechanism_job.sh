#!/bin/bash
set -e
mkdir -p bundle && tar xzf mech_bundle.tar.gz -C bundle
export FLYSIM_CACHE=$PWD/bundle/cache RES_DIR=/data/results_mechanism BASE_RES=/data/results OUT_DIR=$PWD/out
mkdir -p out
python -c "import numpy,scipy,pandas;print(numpy.__version__,scipy.__version__,pandas.__version__)"
nproc; free -g || true
python bundle/code/mech_sweep.py "$@" 2>&1 | tee out/sweep_stdout.log || true
cd out && tar czf mech_results_${TAG}.tar.gz results $(ls drivers.npz lhs_sets.json progress_*.log sweep_stdout.log errors 2>/dev/null) && cd ..
ls -la out | head; du -sh /data/results_mechanism || true
