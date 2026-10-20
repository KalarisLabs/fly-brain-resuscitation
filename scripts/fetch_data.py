"""Check / download the FlyWire Codex files listed in data/README.md.

python scripts/fetch_data.py --check            # which required files are present under $RAW_ROOT
python scripts/fetch_data.py --list fafb        # list Codex data products for a dataset (JSON)
python scripts/fetch_data.py --get [--all]      # download missing files (needs CODEX_API_TOKEN; --all adds large optional files)

Dataset ids and product names follow the Codex download API (see the Codex FAQ). Verify them with --list before use,
because Codex renames or versions products over time. Downloads are written under the local filenames the build
scripts expect.
"""
import json
import os
import sys
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import fbr_paths as P  # noqa: E402

API = "https://codex.flywire.ai/api"
# local filename -> (dataset id, data_product, required)
FILES = {
    "connections_princeton (2).csv.gz": ("fafb", "connections_princeton", True),
    "neurons (2).csv.gz": ("fafb", "neurons", True),
    "classification.csv.gz": ("fafb", "classification", True),
    "consolidated_cell_types.csv.gz": ("fafb", "consolidated_cell_types", True),
    "cell_stats.csv.gz": ("fafb", "cell_stats", True),
    "names.csv.gz": ("fafb", "names", True),
    "connections_princeton (3).csv.gz": ("banc", "connections_princeton", False),
    "neurons (3).csv.gz": ("banc", "neurons", False),
    "connections_princeton (4).csv.gz": ("manc", "connections_princeton", False),
    "neurons (4).csv.gz": ("manc", "neurons", False),
    "connections_princeton.csv.gz": ("mcns", "connections_princeton", False),
    "neurons.csv.gz": ("mcns", "neurons", False),
    "fafb_v783_princeton_synapse_table.csv.gz": ("fafb", "fafb_v783_princeton_synapse_table", False),
    "sk_lod1_783_healed.zip": ("fafb", "sk_lod1_783_healed", False),
}
LARGE = {"fafb_v783_princeton_synapse_table.csv.gz", "sk_lod1_783_healed.zip"}


def check():
    miss = 0
    for f, (ds, prod, req) in FILES.items():
        p = P.raw(f); ok = os.path.exists(p)
        miss += (not ok) and req
        print(f"{'OK ' if ok else '-- '} {'required' if req else 'optional'} {ds:5s} {prod:36s} {p}")
    print(f"RAW_ROOT = {P.RAW}; missing required: {miss}")
    return miss


def list_products(ds):
    with urllib.request.urlopen(f"{API}/download?dataset={urllib.parse.quote(ds)}", timeout=60) as r:
        print(json.dumps(json.loads(r.read().decode()), indent=1)[:20000])


def get(include_large=False):
    tok = os.environ.get("CODEX_API_TOKEN")
    if not tok:
        sys.exit("Set CODEX_API_TOKEN (Codex account -> API token) or download manually: see data/README.md")
    os.makedirs(P.RAW, exist_ok=True)
    for f, (ds, prod, req) in FILES.items():
        if os.path.exists(P.raw(f)) or (f in LARGE and not include_large):
            continue
        q = urllib.parse.urlencode(dict(data_product=prod, dataset=ds, api_token=tok))
        print("downloading", ds, prod, "->", f, flush=True)
        try:
            with urllib.request.urlopen(f"{API}/download_resource?{q}", timeout=600) as r, open(P.raw(f) + ".part", "wb") as o:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    o.write(b)
            os.replace(P.raw(f) + ".part", P.raw(f))
        except Exception as e:  # dataset not publicly downloadable, renamed product, auth ...
            print(f"   FAILED ({type(e).__name__}: {e}); download {f} manually (data/README.md)")


if __name__ == "__main__":
    a = sys.argv[1:]
    if "--list" in a:
        list_products(a[a.index("--list") + 1])
    elif "--get" in a:
        get("--all" in a)
    else:
        sys.exit(1 if check() else 0)
