> Copied from the research tree at repository build time (2026-10-19). All results are **model predictions** of a phenomenological simulator, not evidence about real brains or real resuscitation. Paths refer to the original layout under `$RESEARCH_ROOT`.

# flysim — FAFB v783 whole-brain spiking + metabolic simulator

`flysim.py` (pure numpy/scipy; no brian2/numba — neither is installable on this host) re-implements the
Shiu et al. 2024 (*Nature* 634:210) leaky-integrate-and-fire model of the whole FlyWire brain and adds an
optional slow **metabolic / ionic layer** (ATP, Na/K-pump failure, extracellular K+ per neuropil, spreading
depolarization, depolarization block, injury → permanent death, reversible synaptic failure).

```python
import sys; sys.path.insert(0, r"$RESEARCH_ROOT\sim")
import flysim as fs
m = fs.load_model("princeton5")                 # or "shiu783"; glut_inhibitory=True (Shiu convention)
sug = m.idx(fs.SHIU_SUGAR_GRN); mn9 = m.idx([fs.SHIU_MN9])
r = fs.run(m, 1000, stim=[dict(idx=sug, rate=150)])          # 1 s, pure Shiu LIF
fs.rates(r, mn9)                                               # Hz
```

## 1. Connectome (`load_model`)
| option | matrix | edges | synapses | sign |
|---|---|---|---|---|
| `princeton5` (default) | `connections_princeton (2).csv.gz` aggregated per (pre,post) over neuropils; file is the **≥5-synapse/pair thresholded** Princeton v783 release (single neuropil rows go down to 1) | 3,732,460 | 50,666,648 | presynaptic neuron NT from `neurons (2).csv.gz` (`nt_type`); 18,032 neurons with no neuron-level call get the synapse-weighted majority `nt_type` of their outgoing rows (100 % agreement where both exist); 1,626 neurons with no NT → `unknown_sign` (+1) |
| `shiu783` | Shiu's own `Connectivity_783.parquet` (GitHub philshiu/Drosophila_brain_model) | 15,091,983 (no threshold) | 54,492,922 | Shiu's `Excitatory` column |

Sign convention (switches): GABA −1; GLUT −1 if `glut_inhibitory=True` (Shiu) else +1; ACh +1; DA/SER/OCT
`modulatory_sign` (+1, Shiu); unknown `unknown_sign` (+1). With Shiu convention 61.5 % of synapses (62.3 % of
edges, 70.5 % of neurons) are excitatory; with Glu excitatory 76.5 % of synapses.
NB: Shiu's matrix uses a different synapse prediction set: only 2.70 M of its pairs have ≥5 synapses
(34.2 M synapses) and only 392,587 pair counts equal the Princeton counts, so `princeton5` is a related but not
identical network; `princeton5` is ~2× more excitable (sugar→MN9 curve, figure).

`Model` fields: `N`, out-CSR `indptr/indices/weight` (rows = presynaptic, weight = sign×syn_count),
`neurons` (DataFrame in index order: idx, root_id, nt, super_class, class, sub_class, primary_type, side,
length_nm, area_nm, size_nm, dominant_neuropil, name, group …), `neuropils` (79 compartments incl. `UNASGD`),
`Mn` (N×79 row-normalized synapse fractions), `adj` (79×79 K+ diffusion adjacency), `adj_source`.
Helpers: `m.idx(root_ids)`, `m.select(primary_type, side=None, **column_filters)`.

Caches ($RESEARCH_ROOT\cache): `fafb783_syncount_pre_post_csr.npz`, `fafb783_neurons.csv.gz`,
`fafb783_neuron_neuropil_syn.npz`, `fafb783_neuropils.csv`, `fafb783_connectome_stats.json`,
`shiu783_syncount_pre_post_csr.npz`, `shiu783_presyn_sign.npy`; geometry (streamed from
`fafb_v783_princeton_synapse_table.csv.gz`, 80,215,790 synapses): `fafb783_neuropil_geometry.csv`
(centroid/SD/bbox µm), `fafb783_neuropil_contact.npy` (4-µm voxel face contacts), `fafb783_neuropil_centroid_dist_um.npy`,
`fafb783_neuropil_adjacency_geometric.npy` (contact area / centroid distance, row-max normalized — used automatically
when present), `fafb783_neuron_synapse_centroid.npz`; skeletons (streamed from `sk_lod1_783_healed.zip`):
`fafb783_skeleton_stats.csv.gz` (cable_um, area_um2, soma xyz).

## 2. Fast layer (Shiu LIF, verified against `model.py` on GitHub)
dv/dt = (v_0 − v + g)/τ_m, dg/dt = −g/τ_syn, exact (brian2 `linear`) integration, dt 0.1 ms;
v_0 = v_reset = −52 mV, v_th = −45 mV, τ_m 20 ms, τ_syn 5 ms, refractory 2.2 ms (v **and g** frozen), delay 1.8 ms,
spike → g_post += w_syn·sign·n_syn, w_syn 0.275 mV; Poisson input: independent per stimulated neuron, adds
f_poi·w_syn = 68.75 mV **to v** (not g), stimulated neurons have refractory 0 (as in Shiu). Brian2 schedule order
kept (update → threshold → delayed synaptic delivery + Poisson → reset). Event-driven: only CSR rows of spiking
neurons are gathered; only an *active set* of neurons away from rest is integrated (`lif=dict(active_tol=1e-5)`
mV snaps back to exact rest; `0` disables; identical results observed). Spontaneous activity is zero (as in Shiu).
Runtime on this host (16 slow cores, single process): 13.5 s wall per biological second for sugar-GRN 150 Hz on the
calibrated princeton5 model (16.8 s at w_syn 0.275; ~10 k active neurons); metabolic layer adds ~17.6 s/s when
quiescent and more under activity (72.8 s for the same 1-s sugar run); grows with activity (whole-brain hyperexcitation is much slower); multiprocessing via
pipes is blocked by the sandbox → run independent processes (see `validate_lif.py launch`).

## 3. Slow metabolic / ionic layer (`metab=dict(...)`, defaults `fs.METAB_DEFAULTS`)
Updated every `dt_slow_ms` (5 ms default; 1 ms allowed). `speedup` multiplies all slow rates (time
compression: 1 ms of spiking = `speedup` ms of metabolic time) — needed for minute-scale protocols.
* **ATP** a_i: da/dt = [m·b_i·Q_s·S_i·(1−a)/(K_adp+1−a) − d_i·a]/τ_atp, e_i = a_i/a0.
  Saturating, finite-reserve production (m = `reserve` = 4 × resting demand), K_adp set so a = a0 = 0.95 at
  rest, S = 1. Demand d_i = b_i·Q_d·(1 + r_i/r_half + k_K·(K_i−K0)/K0); b_i = clip((area_i/median)^β, 0.25, 4),
  β = 0.5, area from `cell_stats` (`demand_mode` 'area' | 'length' | 'surface_to_volume' | 'uniform' | any
  neuron-table column, e.g. skeleton `area_um2` after merging); r_i = 100-ms rate estimate; r_half = 100 Hz.
  τ_atp = 20 s → median neuron at rest loses ATP with e-fold time 20 s in anoxia. Temperature:
  Q_d = q10_demand^((T−22)/10) (2.5), Q_s (1.5), pump/clearance Q10 1.8.
  Steady state under partial supply: electrical failure begins at S ≈ 0.15–0.2 (cf. classical ~20 % CBF threshold).
* **Pump failure** dV_atp = 40·((e_ok−e)/e_ok)^1.5 for e < e_ok = 0.75 (dead-zone keeps the LIF identical at S = 1).
* **Extracellular K+** K_c (mM) in 79 neuropil pools, K0 = 5:
  dK_c/dt = k_spk·r̄_c + k_rel·dep_c·R_c·(1−K/K_max) + k_leak·leak(pump_c) − k_clr·pump_c·(K_c−K0) + k_diff·L·K
  (L = Laplacian of the geometric adjacency). Neuron K_i = Σ_c Mn_ic K_c; dV_K = g_K·RT/F·ln(K/K0) mixed over
  compartments (g_K 0.6). Total depolarization soft-capped: v_rest_eff = −52 + 50·tanh((dV_atp+dV_K)/50).
  Regenerative release (spreading depolarization): dep_i = σ((v_rest_eff − V_sd,i)/2), V_sd = −30 ± 2 mV,
  k_rel 60 mM/s; release store R_c depletes (τ 8 s) and recovers with pumps (τ 30 s) so SD in energized tissue
  terminates; clearance k_clr = 1/s × pump (50 % glial, glial ATP follows S).
* **Depolarization block**: h → 1 − σ((v_rest_eff + 32)/2) with τ_h 10 ms; v_th,eff = −45 + 100·(1−h).
* **Injury / death**: dJ/dt = dep·(1−e) − k_rep·e·(1−dep)·J (Ca2+-load proxy, s); hazard
  h_max·J⁴/(J⁴+J50⁴), h_max = 1/30 s⁻¹, J50 = 60 s → death is permanent (mask; threshold +∞, output 0).
* **Synaptic failure** (reversible): presynaptic release fraction = clip((e−0.3)/0.3, 0, 1).
All constants are phenomenological and adjustable; they are **not fitted to fly data** (bibliography track should
calibrate τ_atp, J50 etc.; Drosophila are anoxia tolerant).

## 4. `run(model, duration_ms, ...)`
* `stim=[{'idx', 'rate' (Hz or f(t_ms)), 't_on', 't_off'}]` — `fs.xenova_pulse()` gives the Xenova envelope.
* `metab=None | dict` — enables slow layer. `supply` = scalar | (79,) | f(t_s) (`fs.piecewise`, `fs.ramp`).
  `temperature` = °C or f(t_s). Times in schedules are **slow-layer seconds** (= bio_ms·1e-3·speedup).
* `interventions={...}`: `block_all` (bool or f(t_s): total spike + synaptic block, TTX / isoelectric-coma analog), `k_rel_scale` (SD blocker / NMDA-antagonist analog), `clr_scale` (K+ buffering),
  `injury_scale` (neuroprotection), `atp_floor`, `k_adapt` is a *metab* key (spike-frequency adaptation, mV/Hz, default 0), `v_offset` (mV, scalar or per-neuron), `K_clamp` (mM),
  `silence` (indices), `revive` (counterfactual: resurrect dead neurons). Scalars or f(t_s).
* `state=` continue from `result['state']` (multi-stage protocols: insult → reperfusion → probe).
* `epochs=[(t0,t1), ...]` (ms) per-neuron counts; `groups` column for population rate per ms; `record_idx`
  (indices or 'all') raster; `rec_every_ms` slow-trace sampling; `seed`.
Returns `Result`: `counts` (epochs×N), `pop_ms`, `group_names`, `raster_t_ms/raster_idx`, `state`, `wall_s`,
`mean_active`, and with metab `slow` (t_s, e_mean, e_lt05, K (T×79), e_np, dep_frac_np, dep_frac, block_frac,
dead_frac, vr_mean, syn_frac, J_mean, S_mean, T), `alive`, `e`, `J`.

## 5. Metrics
`fs.rates(res, idx, epoch)`, `fs.run_trials(...)`, `fs.population_rate(res, group, bin_ms)`,
`fs.functional_identity_index(rate_ref, rate_test)` → FII = max(r,0)·√(recall·precision), where r = Pearson
correlation of log(rate+1) over neurons responsive (≥1 Hz) in either condition, recall = fraction of
reference-responsive neurons still responsive, precision = fraction of test-responsive neurons that were
responsive before.

## 6. Scripts / results
* `validate_lif.py launch 7 5` → `results/validation_lif.csv`, `results/validation_mean_rates.npz`
  (mean rate vectors per condition — usable as FII references).
* `metab_tests.py a|bc|focal` → `results/metab_tests.json`, `results/metab_anoxia_reox.npz`, `results/metab_focal.npz`.
* `geom_synapses.py`, `geom_skeletons.py launch 6` → geometry caches (streamed; never load the big files).
* Figures in `figures/`.
Porting: the per-step work is the body of the loop in `run` (dense ops + CSR gather + bincount) and
`slow_update` (dense ops + two sparse mat-mats); both translate directly to torch (`index_add_`, sparse mm).


## 7. Calibration, validation and unit-test results (all numbers from files in `results/`)
**Criticality / w_syn.** With Shiu's literal w_syn = 0.275 mV the Princeton ≥5 matrix is *supercritical*: after a
150-ms sugar-GRN stimulus the active population keeps growing (538 → 1,805 → 2,721 neurons at 0–150 / 200–400 /
1000–1200 ms) whereas Shiu's own matrix returns to silence (326 → 33 → 0). `load_model('princeton5')` therefore
defaults to **w_syn = 0.185** (344 → 0 → 0 active; MN9 107 Hz vs 97 Hz on shiu783 at 150 Hz).
`validate_lif.py` runs with literal 0.275 (`results/validation_lif.csv`, 5 trials × 1 s; figure
`figures/fig_lif_validation_v3.png`); `validate_calibrated.py` repeats sugar/LC4/LC9 at 0.185
(`results/validation_calibrated.csv`, 3 trials; 79/87 jobs completed — turn presets not rerun at 0.185).
| readout (mean Hz) | princeton5 w0.275 | princeton5 w0.185 | shiu783 |
|---|---|---|---|
| MN9, sugar 50 / 100 / 200 Hz | 65.6 / 162.4 / 195.6 | 3.3 / 117.7 / 148.3 | 15.8 / 70.6 / 111.8 |
| DNp01 L/R, LC4 100 Hz | 220 / 236 | 213 / 229 | – |
| DNg100 R, LC9 100 Hz | 124 | 116 | – |
| DNa02 L/R, left LC9+DNa02 180 Hz | 197 / 0 | n.r. | – |

**Metabolic unit tests** (`results/metab_tests.json`, calibrated model, speedup 50):
* (a) S = 1: no-stimulus baseline exactly stationary (0 spikes, K dev 0, e = 1). With activity-K term off
  (`k_spk=0`) spike counts are **bit-identical** to the pure LIF; with default k_spk = 0.05, 312 neurons differ, max
  14 Hz, rate correlation 0.999, MN9 141 vs 140 Hz. PASSED.
* (b) S = 0 from t = 10 s: SD onset 9.25 s after anoxia; at 35 s 99.5 % depolarized/blocked, mean K 39.6 mM,
  synaptic release 0.2 %, firing 0.47 Hz/neuron. PASSED.
* (c) reoxygenation after 30 s anoxia: energy and ions recover (repolarized 4.25 s after reoxygenation) but the
  whole-brain repolarization burst **ignites a self-sustained high-activity attractor** (Shiu network has no
  adaptation/noise-balanced inhibition): firing climbs to tens of Hz/neuron, K re-accumulates (12.8 mM at 100 s),
  MN9 probe 0 Hz, FII 0 → FAIL as a naive protocol. Same after only 7 s anoxia (below SD onset): FAIL.
  −10 mV hyperpolarization 40–60 s does not terminate it: 0.75 Hz/neuron persists during hyperpolarization
  (50–59 s) and 4.36 Hz/neuron after release (80–92.5 s; `pop_Hz_per_neuron_80_92s` = 4.36 in metab_tests.json).
  **Reset protocol** (`block_all` for 20 s after reperfusion) restores the silent baseline:
  30-s anoxia → e 1.00, K 5.00 mM, 0.7 % neurons dead, MN9 113 Hz (pre 127), FII 0.635;
  7-s anoxia → no deaths, FII 0.686. Noise ceiling (same protocol, no anoxia): FII 0.686 (single 150-ms probe;
  repeat probes / longer epochs for tighter FII). PASSED with reset.
* Focal ischemia (right optic lobe AME/LA/LO/LOP/ME_R, S = 0 from 2 s, geometric adjacency): core K peaks
  42.5 mM; 26 of 74 non-core compartments exceed 20 mM, including the contralateral optic lobe (LA_L 28 mM);
  peak K outside the core is not distance-ordered → spread is partly network-driven (spike K) rather than
  diffusive; K shows ~5-s oscillations (release-store depletion/recovery cycles); 1.9 % dead at 60 s.
Figure: `figures/fig_metab_anoxia_v2.png`. Earlier runs at w_syn 0.275 kept as `*_wsyn0.275*`.

## 8. Known limitations
* Metabolic/ionic constants are phenomenological, not fitted to Drosophila data; SD timing (9 s) scales with
  tau_atp; J50/h_max set the death kinetics (Drosophila are far more anoxia tolerant than mammals).
* `speedup` compresses slow vs fast time; minute-scale protocols at speedup 1 are prohibitively slow on this host
  (whole-brain hyperexcitable phases produce ~10^7–10^8 spikes per biological second).
* No spontaneous activity, adaptation or plasticity in the Shiu model: the silent state and the runaway state are
  both attractors; FII depends on returning to the silent attractor (see reset protocol).
* `UNASGD` (unassigned-synapse) compartment has no geometric neighbours (isolated K pool).
* Memory: one model ≈ 0.5–1 GB in RAM; free results between runs (several parallel processes were OOM-killed).
