"""Build all derived caches from the raw downloads, in dependency order.

python scripts/prepare_caches.py --fafb          # simulator caches (bit-exact vs. the original; ~2 min, ~1.5 GB RAM)
python scripts/prepare_caches.py --geometry      # + K+ geometric adjacency (streams the 2.7 GB synapse table)
python scripts/prepare_caches.py --skeletons     # + skeleton stats (streams the 13.9 GB zip; 6 shards) then rebuilds neurons
python scripts/prepare_caches.py --structure     # structure-track graph caches for FAFB, BANC, MANC, MCNS
python scripts/prepare_caches.py --cns           # replication caches for MCNS and BANC (needs --structure first)
python scripts/prepare_caches.py --all
"""
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
B = os.path.join(ROOT, "src", "build")


def run(*args):
    print("+", " ".join(args), flush=True)
    subprocess.run([sys.executable, *args], check=True, cwd=ROOT)


a = set(sys.argv[1:]) or {"--fafb"}
if "--all" in a:
    a |= {"--fafb", "--geometry", "--skeletons", "--structure", "--cns"}
if "--fafb" in a:
    run(os.path.join(B, "build_fafb_cache.py"))
if "--geometry" in a:
    run(os.path.join(B, "geom_synapses.py"))
if "--skeletons" in a:
    run(os.path.join(B, "geom_skeletons.py"), "launch", "6")
    run(os.path.join(B, "build_fafb_cache.py"))          # merge skeleton columns into the neuron table
if "--structure" in a:
    for ds in ["FAFB", "BANC", "MANC", "MCNS"]:
        run(os.path.join(B, "build_graphs.py"), ds)
if "--cns" in a:
    for ds in ["MCNS", "BANC"]:
        run(os.path.join(B, "build_cns.py"), ds)
