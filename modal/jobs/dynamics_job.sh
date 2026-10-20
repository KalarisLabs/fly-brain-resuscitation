#!/bin/bash
set -e
mkdir -p bundle && tar xzf flysim_bundle.tar.gz -C bundle
export FLYSIM_CACHE=$PWD/bundle/cache RES_DIR=/data/results OUT_DIR=$PWD/out
mkdir -p /data/results out
python -c "import numpy,scipy,pandas;print(numpy.__version__,scipy.__version__,pandas.__version__)"
nproc; free -g || true
python bundle/code/sweep.py "$@" 2>&1 | tee out/sweep_stdout.log
ls /data/results | head; du -sh /data/results || true
