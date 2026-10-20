#!/bin/bash
set -e
mkdir -p bundle bundle/cache && tar xzf bundleA_MCNS.tar.gz -C bundle
export FLYSIM_CACHE=$PWD/bundle/cache CNS_CACHE=$PWD/bundle/cnscache OUT_DIR=$PWD/out RES_DIR=/data/results_rep_v3
mkdir -p out /data/results_rep_v3
nproc; free -g || true
python bundle/code/rsweep.py --ds MCNS --wsyn 0.0538 --workers 32 2>&1 | tee out/MCNS_main.log
( python bundle/code/rsweep.py --ds MCNS --variant alln --wsyn 0.0538 --workers 16 > out/MCNS_alln.log 2>&1 ) &
python bundle/code/rsweep.py --ds MCNS --variant wlo --wsyn 0.0434 --workers 15 2>&1 | tee out/MCNS_wlo.log
wait
du -sh /data/results_rep_v3 || true
