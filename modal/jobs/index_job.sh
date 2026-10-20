#!/bin/bash
set -e
mkdir -p R && tar xf index_bundle.tar -C R
export RESEARCH_ROOT=$PWD/R/research RAW_ROOT=$PWD/R/raw OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# (perclib numpy<2 popcount fallback is now built into src/perclib.py; runtime patch removed)
nproc; free -g || true
cd R/research/followup/index; mkdir -p data
python index_fafb.py precompute > data/precompute.log 2>&1; cat data/precompute.log
NS=${NSHARD:-12}
for k in $(seq 0 $((NS-1))); do python index_fafb.py $k $NS > data/shard_$k.log 2>&1 & done
( while sleep 120; do free -g | head -2 | tail -1; ls data/run_indices_*.csv 2>/dev/null | xargs -r wc -l | tail -1; done ) &
MON=$!
wait $(jobs -p | grep -v $MON) 2>/dev/null || true
while pgrep -f 'index_fafb.py [0-9]' >/dev/null; do sleep 10; done
kill $MON 2>/dev/null || true
echo "FAFB shards done"; tail -n 2 data/shard_*.log | tail -n 40
ls data
printf "%s\n" FAFB:random FAFB:metabolic_syn FAFB:degree FAFB:focal BANC:random BANC:metabolic_syn BANC:degree BANC:focal MCNS:random MCNS:metabolic_syn MCNS:degree MCNS:focal MANC:random MANC:metabolic_syn MANC:degree MANC:focal | xargs -P 8 -I{} bash -c 'ds=${1%%:*}; rg=${1##*:}; python cross.py $ds $rg 2 > data/cross_${ds}_${rg}.log 2>&1 || echo FAIL $1' _ {}
echo "cross done"; tail -n 3 data/cross_*.log | tail -n 80
cd -; mkdir -p out; cp R/research/followup/index/data/*.csv R/research/followup/index/data/*.json R/research/followup/index/data/*.log out/ 2>/dev/null || true
mkdir -p /data/results_index && cp out/* /data/results_index/ || true
ls -la out
