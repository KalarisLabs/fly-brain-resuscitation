#!/bin/bash
set -e
mkdir -p bundle/cache bundle/cnscache
for b in b*.tar.gz; do tar xzf $b -C bundle; done
export FLYSIM_CACHE=$PWD/bundle/cache CNS_CACHE=$PWD/bundle/cnscache OUT_DIR=$PWD/out RES_DIR=/data/results_rep_v4
mkdir -p out /data/results_rep_v4
nproc; free -g || true
( while true; do echo "$(date +%s) $(free -m | awk 'NR==2{print $3, $7}')" >> out/mem.log; sleep 10; done ) &
MON=$!
( python bundle/code/rsweep.py --ds FAFBproxy --wsyn 0.185 --workers 8 > out/FAFBproxy.log 2>&1 ; python bundle/code/rsweep.py --ds FAFBproxy --variant alln --wsyn 0.185 --workers 8 > out/FAFBproxy_alln.log 2>&1 ) &
P1=$!
python bundle/code/rsweep.py --ds MCNS --wsyn 0.0538 --workers 24 2>&1 | tee out/MCNS_main.log
( python bundle/code/rsweep.py --ds MCNS --variant alln --wsyn 0.0538 --workers 12 > out/MCNS_alln.log 2>&1 ) &
P2=$!
python bundle/code/rsweep.py --ds MCNS --variant wlo --wsyn 0.0434 --workers 12 2>&1 | tee out/MCNS_wlo.log
wait $P1 $P2 || true
kill $MON || true
sort -k2 -n out/mem.log | tail -1
