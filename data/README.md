# Data (not committed)

Download the files below into `data/raw/`, or set `RAW_ROOT` to wherever they already are. Keep the exact filenames: the build scripts use the names that Codex browser downloads produce, including the ` (2)` suffixes.

| local filename | dataset | Codex data product | size (gz) | used by |
|---|---|---|---|---|
| `connections_princeton (2).csv.gz` | FAFB v783 | connections_princeton (≥5 syn/pair) | 68.5 MB | simulator, structure |
| `neurons (2).csv.gz` | FAFB v783 | neurons | 1.7 MB | simulator, structure |
| `classification.csv.gz` | FAFB v783 | classification | 0.9 MB | simulator, structure |
| `consolidated_cell_types.csv.gz` | FAFB v783 | consolidated_cell_types | 0.9 MB | simulator, structure |
| `cell_stats.csv.gz` | FAFB v783 | cell_stats | 2.5 MB | simulator (metabolic demand), structure |
| `names.csv.gz` | FAFB v783 | names | 1.2 MB | simulator |
| `fafb_v783_princeton_synapse_table.csv.gz` | FAFB v783 | synapse table | 2.7 GB | geometry (K⁺ adjacency), optional |
| `sk_lod1_783_healed.zip` | FAFB v783 | skeletons (SWC) | 13.9 GB | skeleton stats, optional |
| `connections_princeton (3).csv.gz`, `neurons (3).csv.gz` | BANC | connections / neurons | 29.2 / 2.9 MB | structure, replication |
| `connections_princeton (4).csv.gz`, `neurons (4).csv.gz` | MANC v1.2.1 | connections / neurons | 26.7 / 1.0 MB | structure |
| `connections_princeton.csv.gz`, `neurons.csv.gz` | MCNS | connections / neurons | 45.1 / 2.3 MB | structure, replication |

Optional, for the `shiu783` comparison matrix: `Connectivity_783.parquet` and `Completeness_783.csv` from github.com/philshiu/Drosophila_brain_model.

How to get them:
1. **Browser:** open https://codex.flywire.ai/, select the dataset, then go to *Info → Download Data*.
2. **Programmatic:** `python scripts/fetch_data.py --list fafb` lists the available products. With `CODEX_API_TOKEN` set, `python scripts/fetch_data.py --get` downloads them (Codex FAQ: `/api/download?dataset=…`, `/api/download_resource?data_product=…&dataset=…&api_token=…`). Codex currently offers public downloads for FAFB and BANC. For MANC and MCNS, use the Janelia homepages linked in the Codex FAQ, or their Codex pages if your account has access.

All files are CC-BY; see [../DATA_LICENSE.md](../DATA_LICENSE.md). Large or rare caches can be rebuilt with `scripts/prepare_caches.py`.
