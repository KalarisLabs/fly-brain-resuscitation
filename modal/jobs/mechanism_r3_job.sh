#!/bin/bash
set -e
mkdir -p bundle && tar xzf r3_bundle.tar.gz -C bundle
export FLYSIM_CACHE=$PWD/bundle/cache RES_DIR=/data/results_mechanism BASE_RES=/data/none OUT_DIR=$PWD/out
mkdir -p out
nproc; free -g || true
python bundle/code/mech_round3.py "$@" 2>&1 | tee out/round3_stdout.log || true
cd out && tar czf mech_r3_${TAG}.tar.gz results $(ls bisection.json drivers_*.npz round3_*.log errors 2>/dev/null) && cd ..
ls -la out | head; du -sh /data/results_mechanism* || true
