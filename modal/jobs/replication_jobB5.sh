#!/bin/bash
set -e
mkdir -p bundle bundle/cache && tar xzf bundleB_BANC_FAFB.tar.gz -C bundle
export FLYSIM_CACHE=$PWD/bundle/cache CNS_CACHE=$PWD/bundle/cnscache OUT_DIR=$PWD/out RES_DIR=/data/results_rep_v3
mkdir -p out /data/results_rep_v3
nproc; free -g || true
( python bundle/code/rsweep.py --ds FAFBproxy --wsyn 0.185 --workers 10 > out/FAFBproxy.log 2>&1 ; python bundle/code/rsweep.py --ds FAFBproxy --variant alln --wsyn 0.185 --workers 10 > out/FAFBproxy_alln.log 2>&1 ) &
python bundle/code/rsweep.py --ds BANC --wsyn 0.2227 --workers 22 2>&1 | tee out/BANC_main.log
( python bundle/code/rsweep.py --ds BANC --variant alln --wsyn 0.2227 --workers 11 > out/BANC_alln.log 2>&1 ) &
python bundle/code/rsweep.py --ds BANC --variant wlo --wsyn 0.1798 --workers 11 2>&1 | tee out/BANC_wlo.log
wait
du -sh /data/results_rep_v3 || true
