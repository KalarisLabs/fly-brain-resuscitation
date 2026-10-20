#!/bin/bash
set -e
mkdir -p bundle && tar xzf rep_bundle3.tar.gz -C bundle
export FLYSIM_CACHE=$PWD/bundle/cache CNS_CACHE=$PWD/bundle/cnscache OUT_DIR=$PWD/out RES_DIR=/data/results_replication_v2
mkdir -p out /data/results_replication_v2
nproc; free -g || true
( python bundle/code/rsweep.py --ds FAFBproxy --wsyn 0.185 --workers 10 > out/FAFBproxy.log 2>&1 ) &
P1=$!
python bundle/code/rsweep.py --ds BANC --wsyn 0.2227 --workers 22 2>&1 | tee out/BANC_main.log
python bundle/code/rsweep.py --ds BANC --variant wlo --wsyn 0.1798 --workers 15 2>&1 | tee out/BANC_wlo.log
wait $P1 || true
du -sh /data/results_replication_v2 || true
