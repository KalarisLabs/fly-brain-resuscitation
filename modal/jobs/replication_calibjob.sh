#!/bin/bash
set -e
mkdir -p bundle && tar xzf rep_bundle.tar.gz -C bundle
export FLYSIM_CACHE=$PWD/bundle/cache CNS_CACHE=$PWD/bundle/cnscache OUT_DIR=$PWD/out
mkdir -p out /data/results_replication
nproc; free -g || true
python bundle/code/calib.py --workers 32 --ds FAFB,MCNS,BANC 2>&1 | tee out/calib_stdout.log
mkdir -p /data/results_replication/calib && cp -r out/calib/* /data/results_replication/calib/ || true
