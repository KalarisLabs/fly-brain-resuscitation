# Methods

This file collects the methods sections of the track reports in [docs/reports/](reports/). Section text is reproduced from those reports; the reports remain authoritative.

## A. Simulator (`src/flysim.py`, `src/flysim_dyn.py`, `src/flysim_mech.py`)

`flysim_dyn` adds options to `flysim`, and `flysim_mech` adds options to `flysim_dyn`. Every added option defaults to the original numerics; `tests/test_lif_baseline.py::test_engine_versions_identical_with_default_options` checks this.

### 1. Connectome (`load_model`)
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

### 2. Fast layer (Shiu LIF, verified against `model.py` on GitHub)
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

### 3. Slow metabolic / ionic layer (`metab=dict(...)`, defaults `fs.METAB_DEFAULTS`)
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

### 4. `run(model, duration_ms, ...)`
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

### 5. Metrics
`fs.rates(res, idx, epoch)`, `fs.run_trials(...)`, `fs.population_rate(res, group, bin_ms)`,
`fs.functional_identity_index(rate_ref, rate_test)` → FII = max(r,0)·√(recall·precision), where r = Pearson
correlation of log(rate+1) over neurons responsive (≥1 Hz) in either condition, recall = fraction of
reference-responsive neurons still responsive, precision = fraction of test-responsive neurons that were
responsive before.


## B. Structure track (percolation, `src/perclib.py`)

### 1. Data and pairing QC

| dataset | neurons | conn. ids found in annotation | annotated ids present in conn. | pair-level synapse floor in file | pairs (all) | pairs ≥5 | sources (sensory) | sinks | of which DN / MN | removable interneurons | modality groups |
|---|---|---|---|---|---|---|---|---|---|---|---|
| FAFB | 139255 | 100.00% | 99.52% | ≥5 | 3732460 | 3732460 | 17550 | 1415 | 1305 / 110 | 120290 | 22 |
| BANC | 158262 | 100.00% | 97.28% | ≥3 | 3037361 | 1528585 | 17086 | 2347 | 1316 / 1031 | 138356 | 50 |
| MCNS | 166700 | 100.00% | 99.48% | ≥5 | 6242118 | 6242118 | 17937 | 2241 | 1320 / 921 | 146522 | 41 |
| MANC | 23665 | 100.00% | 99.90% | ≥1 | 5305638 | 1372588 | 6460 | 838 | 0 / 838 | 16352 | 21 |

- Files: FAFB v783 = `connections_princeton (2)` + `neurons (2)` + classification/consolidated_cell_types/cell_stats; BANC v888 = `(3)`; MANC v1.2.1 = `(4)`; MCNS v1.0 = `connections_princeton.csv.gz` + `neurons.csv.gz` (the `(1)` copies are byte-identical and were ignored). The per-neuropil rows were summed per (pre, post) pair. The pair counts match the codex totals: FAFB 3,732,460; BANC 3,037,361; MANC 5,305,638.
- **MAOL v1.1 is not available locally and was not analysed.**
- Harmonised edge set: pairs with ≥5 synapses. FAFB and MCNS files are already floored at 5, so a ≥1 sensitivity analysis was possible only for MANC; BANC is floored at 3. A ≥10 sensitivity analysis was run for all four datasets.
- Roles: **sources** = all super-classes containing "sensory" (FAFB sensory + sensory_ascending; MCNS cb/ol/vnc_sensory etc.). **Sinks** = descending + motor (FAFB); descending + motor + efferent (BANC, MCNS); motor + efferent (MANC; MANC DNs are interneurons there). Glia, trachea and not_a_neuron (BANC, 473 total; MANC glia 15) were excluded. **Removable** = everything else. Damage fraction f is the fraction of removable neurons removed; sources and sinks are never removed.

### 2. Methods (as run)

- **Primary metric `Reach_N(L)`** = fraction of (source neuron, sink neuron) pairs connected by a directed path of ≤L hops through surviving neurons, as the preregistration specifies. Paths are *strict*: they may relay only through removable interneurons, so sensory→sensory and DN→DN chains cannot act as an undamageable backbone. Sources are a fixed sample of 128 sensory neurons stratified by modality group (class × side). L = 4, 6 (primary) and "∞" (iterated to convergence, capped at 30 hops in reduced runs; the uncapped first runs used 60).
- **Secondary metrics:** `Reach_G` (modality-group × sink pairs, strict); `Reach_G,loose` (relays through any neuron); **path redundancy** = Σ_sinks log(1 + number of walks of length 1–6 from the whole source set) on the binarised strict graph; **weighted influence** = Σ_sinks of summed products of input-normalised weights over walks ≤6 hops (prereg "weighted variant"); **giant SCC** fraction (scipy csgraph, computed for the first 3 orders of each regime). Sink subsets DN vs MN are also reported in `curves_normalised_all.csv.gz`.
- Engine: pure numpy bitset propagation (64 sources per uint64 word) with frontier-delta updates, plus scipy.sparse; no networkx on whole graphs. One dataset was processed per process with int32 indices.
- **Damage regimes:** (a) random; (b) metabolic-vulnerability-weighted, removal probability ∝ demand (Efraimidis–Spirakis weighted order). Demand is membrane area (`cell_stats.area_nm`, FAFB only, `metabolic_area`) or total synapse count (all datasets, `metabolic_syn`). (c) Hub-targeted static total degree in the ≥5 graph (`degree`), plus a fine low-f grid (`degree_finegrid`, f ≤ 0.3, 19 points). (d) Focal neuropil lesion: grows from a random core neuropil along the neuropil co-innervation (cosine) graph; neurons are removed by the rank of their dominant neuropil. (e) H2 control: degree-bin-matched random (40 log-degree bins).
- f grid: 24 points from 0 to 0.9 for the first orders, then a 13-point subset (all points are on the 24-point grid).
- **Fits:** pooled points from all orders. Models: 4-parameter logistic, linear, and 3-parameter exponential, compared by Gaussian AIC. ΔAIC = AIC_logistic − min(AIC_linear, AIC_exp). f_c = logistic midpoint; w = 10–90% width = 2·ln9·s. H1 categories follow the preregistered thresholds exactly. D50 = first f at which the normalised curve falls to ≤0.5, per order; it is censored when this never happens by f = 0.9.
- **Ratios targeted/random:** bootstrap over removal orders with logistic refits. There were 2000 resamples for the primary metric on the ≥5 graph and 100 for controls and secondary metrics; the CI is the 2.5–97.5 percentile. Note that degree-targeted orders differ only by tie-breaking, so ratio CIs mostly reflect random-arm variability and order sampling, not connectome uncertainty.
- **Controls:** degree-preserving directed rewired null (target shuffling, duplicates and self-loops removed; FAFB and MANC); edge cut-offs ≥10 (all) and ≥1 (MANC). Negative control: remove all removable neurons lying on no strict source→sink path.
- **Keystones:** candidate types = top-150 by flow centrality, i.e. the sum over a type's neurons of (walks ≤3 from sources) × (walks ≤3 to sinks), strict. Then each other dataset's top-50 names present in a given dataset were also evaluated there (replication pass). For each type, all its removable neurons were removed and relative drops in Reach_N(L6, ∞), Reach_G, redundancy and influence were recorded. Ranking: relative ΔReach_N(L6), then relative Δinfluence. Names were matched exactly, case-insensitive.


## C. Dynamics track (insult / brain death / resuscitation)

### 1. Methods (as run)
* **Protocol** (slow-layer model seconds, speedup 50 = 50 model-s per biological s of spiking): 2 s normoxia → insult (D s) → 60 s reperfusion protocol → 20 s free observation (S = 1, no interventions) → evaluation at speedup 1.
* **Probe battery** (7 stimuli, each run for 250 ms from a copy of the post-protocol state): sugar GRNs 150 Hz, bitter GRNs 150 Hz, JO auditory 150 Hz, JO wind/gravity 150 Hz, LC4 bilateral 100 Hz, LC9 bilateral 100 Hz, left turn (LC9_L + DNa02_L) 180 Hz. Readouts: all 1,415 descending + motor neurons (+MN9). Response = probe rate − rate in a 500-ms unstimulated window from the same state.
* **FII** (preregistered matrix version): Pearson r of vec(R_condition) vs vec(R_baseline) over entries > 0 in either. Baseline = mean of 10 intact seeds. **Noise ceiling** = median leave-one-out single-seed r = **0.9940** (range 0.9888–0.9954); FII_norm = FII / ceiling. Controls: readout-shuffled baseline r = -0.220 (floor); 0.5× gain-scaled baseline r = 1.000 (FII is gain-insensitive, as required for H5-R3). Neuron-level FII (`fs.functional_identity_index`) and Spearman FII are also in `all_runs.csv`.
* **Other metrics:** dead fraction; Reach = fraction of (probe, readout) pairs joined by a ≤6-hop path through surviving neurons (baseline 0.9866); Act = number of neurons evoked ≥8 Hz above spontaneous by any probe (baseline 3466); FSA = fraction of neurons spiking in the unstimulated 500-ms window; SD onset = first time >50 % of neurons depolarized; K⁺ peak; reperfusion surge = peak population rate.
* **Outcome classes:** *runaway* (late observation rate >0.5 Hz/neuron or FSA >0.05); *brain death* (FSA <0.01 **and** no readout response >8 Hz, i.e. silent and unresponsive after the protocol has ended); *recovered* (FII_norm ≥0.8 and Reach ≥0.8×baseline); *structural loss* (dead ≥10 %, not recovered); *impaired* (otherwise).
* **Calibrated parameter set** used for arms B–E/H: τ_ATP 30 s, J50 600 s, q10_injury 2.5 (= q10_demand), all other `METAB_DEFAULTS`, speedup 50. Seeds: 5 for B and C, 3 for A, D, E, H; bootstrap 95 % CIs (2,000 resamples) over seeds.


## D. Mechanism follow-up

### 1. Common methods
* **Protocol** (as in the dynamics track): 2 s normoxia → anoxia D s → 60 s reperfusion protocol → 20 s free observation → probe battery. Times are model seconds; speedup 50 means 50 model-s per biological second of spiking. The battery is 7 probes × 1,415 descending/motor readouts.
* **Calibrated set:** princeton5, w_syn 0.185, τ_ATP 30 s, J50 600 s, q10_injury 2.5, 22 °C.
* **FII_norm.**
  * Arms X/T/M: FII is computed against the 10-seed dynamics-track baseline, with ceiling 0.9940.
  * **Arm L:** each parameter set has its own intact baseline (insult none, 4 seeds 100–103, which do not overlap the run seeds). The ceiling is the median leave-one-out r. Per-set ceilings ranged 0.985–0.999 (median 0.993).
* **Classes:**
  * *runaway*: late observation rate > 0.5 Hz/neuron or FSA > 0.05, the dynamics-track rule.
  * *rescued*: FII_norm ≥ 0.8. For arm T, rescue also requires not runaway.
* **Statistics:** bootstrap 95 % CIs (2,000–4,000 resamples) over seeds, or over parameter sets for the L fractions. Paired-by-seed bootstrap for protocol contrasts. PRCC on ranks with Holm correction over all 250 (set-level outcome × parameter) tests.


## E. Replication models (MCNS v1.0, BANC v888)

### 1. Models — `flysim_cns.py` (flysim.py / flysim_dyn.py untouched)
| | MCNS | BANC | FAFB princeton5 |
|---|---|---|---|
| neurons | 166,700 | 157,789 (473 non-neurons dropped) | 139,255 |
| edges / synapses | 6,242,085 / 89,859,938 | 3,036,600 / 23,553,107 | 3,732,460 / 50,666,648 |
| pair floor in release | ≥5 | ≥3 | ≥5 |
| excitatory neurons / synapses | 69.8 % / 62.5 % | 66.4 % / 56.7 % | 70.5 % / 61.5 % |
| synapses from unknown-NT / HIST neurons | 1.81 % / 0.36 % | 0.01 % / 2.77 % | – |
| K+ compartments | 83 (from 134 neuropils) | 82 (from 111) | 79 |
| neurons by region | OL 101,097 · CB 44,523 · VNC 21,023 · NA 57 | OL 88,799 · CB 42,544 · VNC 24,177 · NA 2,269 | brain only |

Conventions and substitutions:
* **Sign:** presynaptic predicted NT (BANC falls back to the verified NT only where no prediction exists).
  - GABA −1 and Glu −1 (Shiu); HIST −1; ACh +1; DA/SER/OCT/TYR +1; unknown +1.
  - Xenova instead uses HIST 0 and unknown 0, on the unthresholded MCNS graph at w_syn 0.275.
  - The edge-level NT columns are empty, so there is no edge-majority fallback.
* **Demand proxy** = total synapse count. Neither release has size data.
* **K+ adjacency** = connection-based Jaccard 9-NN (no geometry available; see §4).
* Neuropils with <0.1 % of the synapse mass are merged into per-region MISC pools.

### 2. w_syn calibration (pure LIF; criterion = no activity at 1.0–1.2 s after a 150-ms probe)
First reverberating w_syn per probe (coarse ×1.15 scan with 2 seeds, then fine ×1.02 scan with 3 seeds):

| probe | FAFB | MCNS | BANC |
|---|---|---|---|
| sugar | 0.2713 | 0.2014 | 0.6568 |
| bitter | 0.3518 | 0.1407 | 0.6979 |
| JO auditory | 0.3059 | 0.2140 | 0.2569 |
| JO wind/gravity | 0.2221 | 0.1064 | 0.3990 |
| LC4 | 0.2265 | 0.0620 | 0.3470 |
| LC9 | 0.2660 | 0.1861 | 0.3017 |
| turnL | 0.2134 | 0.1618 | 0.3470 |

* The FAFB procedure check passes: Shiu's 0.275 is supercritical on princeton5, and 0.185 = 0.867 × the any-probe threshold.
* **Matched w_syn: MCNS 0.0538, BANC 0.2227.** The sensitivity values 0.7 × threshold (MCNS 0.0434, BANC 0.1798) were prepared but not run.

Responses at the calibrated w_syn (Hz, 1 seed):
| readout | FAFB 0.185 | MCNS 0.0538 | BANC 0.2227 |
|---|---|---|---|
| MN9 ← sugar 50/100/150/200 Hz | 56/92/140/136 | 0/0/0/0 | 0/0/0/0 |
| DNp01 L/R ← LC4 | 208/216 | 204/172 | 0/0 |
| DNg100 L/R ← LC9 | 108/116 | 48/44 | 12/8 |
| DNp09 L/R ← LC9 | 28/16 | 68/64 | 0/72 |
| DNa02 L/R ← turnL | 176/0 | 144/0 | 184/0 |

### 3. Probe battery
| probe | MCNS (n) | BANC (n) |
|---|---|---|
| sugar = left LB3a–d | 38 | 62 |
| bitter = LB1a–e | 57 | 72 |
| JO auditory | 114 | 511 |
| JO wind/gravity | 475 | 475 |
| LC4 / LC9 | 126 / 219 | 114 / 164 |
| turnL (LC9_L + DNa02_L) | 105 | 84 |
| brain readouts (DN + brain MNs; FAFB 1,415) | 1,439 | 1,435 |
| VNC motor/efferent readouts | 802 | 840 |

* All named readouts are present (MN9, DNp01, DNg100, DNp09, DNa02, DNa01, DNg97, MDN).
* Output synapses per sugar GRN: FAFB 405, MCNS 214, **BANC 8**. Sugar→MN9 is silent in both CNS models; LC4→DNp01 is silent in BANC.

### 4. Bug: dense connection-based K+ adjacency is numerically unstable
* **v1 MCNS sweep:** every run, including the 0-s sham, ended in runaway with K+ at 80 mM.
* **Local reproduction:** a period-2 K+ oscillation (5.0→7.4→5.3→23.7→7.8→29 mM) appears before any spike.
* **Cause:** k_diff·dt·λmax exceeds the Euler limit of 2. Shared-neuron adjacency gives 2.87 (FAFB), 4.82 (MCNS), 3.65 (BANC). FAFB geometric gives 0.89, so the FAFB dynamics track is unaffected.
* **Fix (jknn9):** keep the 9 nearest neighbours by Jaccard overlap, union-symmetrise, row-max normalise, and rescale to the FAFB geometric median row sum (2.62).
  - Stability factor 0.69 / 0.72 / 0.78 (FAFB / MCNS / BANC).
  - On FAFB it recovers geometric contacts with precision 0.59 and recall 0.68 (chance 0.135).
  - The MCNS sham is now exactly stationary.
  - `check_stability()` guards future loads.

