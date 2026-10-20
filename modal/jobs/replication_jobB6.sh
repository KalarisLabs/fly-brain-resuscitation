#!/bin/bash
set -e
mkdir -p bundle/cache bundle/cnscache
for b in b*.tar.gz; do tar xzf $b -C bundle; done
export FLYSIM_CACHE=$PWD/bundle/cache CNS_CACHE=$PWD/bundle/cnscache OUT_DIR=$PWD/out RES_DIR=/data/results_rep_v4
mkdir -p out /data/results_rep_v4
nproc; free -g || true
( while true; do echo "$(date +%s) $(free -m | awk 'NR==2{print $3, $7}')" >> out/mem.log; sleep 10; done ) &
MON=$!
python bundle/code/rsweep.py --ds BANC --wsyn 0.2227 --workers 32 2>&1 | tee out/BANC_main.log
( python bundle/code/rsweep.py --ds BANC --variant alln --wsyn 0.2227 --workers 16 > out/BANC_alln.log 2>&1 ) &
P2=$!
python bundle/code/rsweep.py --ds BANC --variant wlo --wsyn 0.1798 --workers 16 2>&1 | tee out/BANC_wlo.log
wait $P2 || true
kill $MON || true
sort -k2 -n out/mem.log | tail -1
