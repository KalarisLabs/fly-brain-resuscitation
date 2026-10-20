# Reproducing tables and figures

Two levels of reproduction are supported.

## Level 1: figures from committed tables (minutes, no FlyWire data)

```bash
python scripts/reproduce_figures.py             # stages results/ into $RESEARCH_ROOT (default ./outputs), then:
#   outputs/structure/figures/fig1..fig4_*.png|pdf   original code: experiments/structure/make_figures.py
#   outputs/replots/R1..R4_*.png                     table-derived re-plots: scripts/replot_from_tables.py
pytest -q                                         # 26 unit tests on a 300-neuron toy connectome
python experiments/smoke/toy_protocol.py          # insult -> reperfusion + block -> probe, toy net
```

| figure in docs/figures | reproduced by | from table(s) |
|---|---|---|
| structure/fig1_reach_vs_damage | `experiments/structure/make_figures.py` (original) | `results/structure/curves_normalised_all.csv.gz`, `fits_all.csv` |
| structure/fig2_targeted_vs_random_ratio | same | `results/structure/ratios_all.csv` |
| structure/fig3_metrics_and_controls | same | `results/structure/*_curves_controls.csv`, `negative_control_all.csv` |
| structure/fig4_keystone_types | same | `results/structure/*_keystone.csv` |
| dynamics/fig1_brain_death_map | `scripts/replot_from_tables.py` -> R1 (re-plot) | `results/dynamics/B_brain_death_map.csv` |
| dynamics/fig3_resuscitation_ranking | re-plot R2 | `results/dynamics/C_resuscitation_ranking.csv` |
| mechanism/fig_M1_robustness | re-plot R3 | `results/mechanism/L_fractions.csv` (+ `L_set_summary.csv`) |
| replication (paper cross-connectome panel) | re-plot R4 | `results/replication/v4_condition_summary.csv` |
| dynamics/fig2, fig4, fig5; mechanism/fig_M2-M4; sim/* | not reproduced (see note) | `results/dynamics/condition_summary.csv`, `results/mechanism/{X,N,T,Y}_*.csv`, `results/sim/*` |

Note: the original plotting code for the dynamics, mechanism and simulator-validation figures was run interactively and was not preserved. The PNGs in `docs/figures/` are the published versions. The re-plots draw the same stored values with new code and are not pixel-identical.

## Level 2: full pipeline from FlyWire data (hours; Modal or a 16+ core machine)

```bash
python scripts/fetch_data.py --check              # what to download from https://codex.flywire.ai/ and where to put it
python scripts/prepare_caches.py --all            # simulator caches, geometry, structure graphs, CNS caches
```

| track | command(s) | writes (under $RESEARCH_ROOT) |
|---|---|---|
| Simulator validation | `python experiments/sim/validate_lif.py 8 30`; `python experiments/sim/validate_calibrated.py`; `python experiments/sim/metab_tests.py a` (also `bc`, `focal`) | `sim/results/validation_*.csv`, `metab_tests.json` |
| Structure (4 connectomes) | `python experiments/structure/run_ds.py <DS>` (+ `controls`), `fine_degree.py <DS>`, `keystone.py <DS>`, `keystone_replication.py prep|final`, `analyze.py` | `structure/tables/*` |
| Dynamics (FAFB) | `python experiments/dynamics/sweep.py --shard K --nshard N --workers W` (Modal: `modal/jobs/dynamics_job.sh`), then `python experiments/dynamics/analyze.py <results_dir>` | `dynamics/tables/*` |
| Mechanism (FAFB) | `python experiments/mechanism/mech_sweep.py --shard A|B --workers W`; round 3: `mech_round3.py --group A|B`; local sign arm: `local_sign.py 1 K N`; analysis: `analyze_load.py`, `analyze_sign.py` | `followup/mechanism/tables/*` |
| Replication (MCNS, BANC, FAFB-proxy) | `python experiments/replication/calib.py --workers 32`; `rsweep.py --ds MCNS --wsyn <w> --workers 32`; `analyze_rep.py` | `followup/replication/tables/*` |
| Recoverability index | `experiments/index/index_fafb.py`, `run_model.py`, `model_index.py`, `cross.py <DS> <regime>` | `followup/index/tables/*` |

All paths resolve through `src/fbr_paths.py` (`RAW_ROOT`, `RESEARCH_ROOT`, `FLYSIM_CACHE`). On hosts that block `multiprocessing` pipes, use `experiments/sim/launch.py` (independent processes).
