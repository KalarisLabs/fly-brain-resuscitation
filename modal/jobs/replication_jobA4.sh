#!/bin/bash
set -e
mkdir -p bundle && tar xzf /data/bundles/rep_bundle3.tar.gz -C bundle
export FLYSIM_CACHE=$PWD/bundle/cache CNS_CACHE=$PWD/bundle/cnscache OUT_DIR=$PWD/out RES_DIR=/data/results_replication_v2
mkdir -p out /data/results_replication_v2
nproc; free -g || true
python bundle/code/rsweep.py --ds MCNS --wsyn 0.0538 --workers 32 2>&1 | tee out/MCNS_main.log
python bundle/code/rsweep.py --ds MCNS --variant wlo --wsyn 0.0434 --workers 15 2>&1 | tee out/MCNS_wlo.log
du -sh /data/results_replication_v2 || true
