> Copied from the research tree at repository build time (2026-10-19). All results are **model predictions** of a phenomenological simulator, not evidence about real brains or real resuscitation. Paths refer to the original layout under `$RESEARCH_ROOT`.

# Dynamics track — insult, brain-death and resuscitation simulations (FAFB v783)

**Status:** sweep complete, 1,209 condition-seed runs (274 conditions), 0 errors. Model: `flysim_dyn` = Shiu-2024 LIF on Princeton v783 ≥5-synapse matrix (w_syn 0.185) + phenomenological metabolic/ionic layer. **All results describe a model with metabolic constants that are not fitted to Drosophila data.** They are model-mechanistic statements, not evidence about real brain death or human resuscitation.

## 0. Provenance
* Code: `dynamics/flysim_dyn.py` (copy of `sim/flysim.py` plus: `demand_scale`, per-neuron `atp_floor`/`injury_scale`, `block_mask`, `coma_arrest`, `q10_injury`, `inh_gain`), `dynamics/sweep.py` (grid, protocols, probe battery, Reach), `dynamics/analyze.py`, `dynamics/flyjob.sh`.
* Compute: Modal, 4 × (32 CPU, 64 GiB) containers + 1 smoke job. Jobs: smoke 79482286-adfb-4c8d-b57c-64d67726799e; shards a4ddd6e8-9ca3-4803-a889-fa5fdd58ccf9, 31c05682-cf05-4a51-b769-97b562f60172, 4d88396b-d4af-4b98-935c-ea62f206acfb, 937a47ef-8836-46d7-8764-a2a4d3655763 (wall 513–922 s). Total 24.3 CPU-h; max run 511 s, max RSS 533 MB. All handles closed and confirmed in `host.compute.ledger()`.
* Raw results: Modal Volume `flysim-data:/results/<arm>/<id>.json|.npz` (full neuron-level sparse spike counts and alive masks) and local copy `dynamics/raw/<arm>/` (json + readout matrix + traces). Tables: `dynamics/tables/` (`all_runs.csv`, `condition_summary.csv`, `B_brain_death_map.csv`, `C_resuscitation_ranking.csv`, `analysis_results.json`).

## 1. Methods (as run)
* **Protocol** (slow-layer model seconds, speedup 50 = 50 model-s per biological s of spiking): 2 s normoxia → insult (D s) → 60 s reperfusion protocol → 20 s free observation (S = 1, no interventions) → evaluation at speedup 1.
* **Probe battery** (7 stimuli, each run for 250 ms from a copy of the post-protocol state): sugar GRNs 150 Hz, bitter GRNs 150 Hz, JO auditory 150 Hz, JO wind/gravity 150 Hz, LC4 bilateral 100 Hz, LC9 bilateral 100 Hz, left turn (LC9_L + DNa02_L) 180 Hz. Readouts: all 1,415 descending + motor neurons (+MN9). Response = probe rate − rate in a 500-ms unstimulated window from the same state.
* **FII** (preregistered matrix version): Pearson r of vec(R_condition) vs vec(R_baseline) over entries > 0 in either. Baseline = mean of 10 intact seeds. **Noise ceiling** = median leave-one-out single-seed r = **0.9940** (range 0.9888–0.9954); FII_norm = FII / ceiling. Controls: readout-shuffled baseline r = -0.220 (floor); 0.5× gain-scaled baseline r = 1.000 (FII is gain-insensitive, as required for H5-R3). Neuron-level FII (`fs.functional_identity_index`) and Spearman FII are also in `all_runs.csv`.
* **Other metrics:** dead fraction; Reach = fraction of (probe, readout) pairs joined by a ≤6-hop path through surviving neurons (baseline 0.9866); Act = number of neurons evoked ≥8 Hz above spontaneous by any probe (baseline 3466); FSA = fraction of neurons spiking in the unstimulated 500-ms window; SD onset = first time >50 % of neurons depolarized; K⁺ peak; reperfusion surge = peak population rate.
* **Outcome classes:** *runaway* (late observation rate >0.5 Hz/neuron or FSA >0.05); *brain death* (FSA <0.01 **and** no readout response >8 Hz, i.e. silent and unresponsive after the protocol has ended); *recovered* (FII_norm ≥0.8 and Reach ≥0.8×baseline); *structural loss* (dead ≥10 %, not recovered); *impaired* (otherwise).
* **Calibrated parameter set** used for arms B–E/H: τ_ATP 30 s, J50 600 s, q10_injury 2.5 (= q10_demand), all other `METAB_DEFAULTS`, speedup 50. Seeds: 5 for B and C, 3 for A, D, E, H; bootstrap 95 % CIs (2,000 resamples) over seeds.

## 2. Arm A — calibration
| condition | SD onset (model s) | dead | FII_norm |
|---|---|---|---|
| tau20_sp25_D60_reset | 9.0 | 0.000 | 1.001 |
| tau20_sp50_D60_reset | 9.5 | 0.000 | 1.001 |
| tau20_sp100_D60_reset | 10.0 | 0.000 | 0.998 |
| tau30_sp25_D60_reset | 11.5 | 0.000 | 1.001 |
| tau30_sp50_D60_reset | 12.5 | 0.000 | 1.001 |
| tau30_sp100_D60_reset | 13.5 | 0.000 | 0.998 |
| tau40_sp25_D60_reset | 14.0 | 0.000 | 1.001 |
| tau40_sp50_D60_reset | 15.5 | 0.000 | 1.001 |
| tau40_sp100_D60_reset | 16.5 | 0.000 | 0.998 |
| J50_60_D120_reset | 12.5 | 0.907 | 0.131 |
| J50_300_D120_reset | 12.5 | 0.017 | 0.951 |
| sham_naive | – | 0.000 | 1.001 |
| sham_reset | – | 0.000 | 1.001 |
| sp25_D30_naive | 11.5 | 0.000 | -0.104 |
| sp100_D30_naive | 13.5 | 0.000 | 0.065 |

* **SD onset 12.5 s at 22 °C with τ_ATP 30 s** (target 10–20 s; literature 14.2 s). The onset scales with τ_ATP (9.0–10.0 s at τ 20 s, 14.0–16.5 s at τ 40 s). It moves by only 2–3 s across a 4-fold change in speedup (25→100), so time compression does not drive SD timing. SD onset does not vary across seeds (it is quantised by the 0.5-s trace sampling).
* **Brief anoxia is structurally survivable** with J50 600 s: 120 s → 0.1 % dead (B map). J50 60 s, the simulator default, killed 90.6 % after 120 s, so it was rejected. J50 300 s → 1.7 %.
* **The runaway after naive reperfusion is not a compression artefact.** It appears at speedup 25, 50 and 100 (FII_norm −0.10 / −0.06 / 0.065 at D = 30 s). Sham protocols give FII_norm 1.001, so neither the protocol nor the 20-s block changes identity on its own.
* **Fly-like metabolic arrest does not by itself give hours-scale tolerance.** 'Coma arrest' cut demand to 5 % in depolarization-blocked neurons. With J50 3,600 s and 960 s anoxia, it halved death (1.5 % vs 3.2 %). With J50 600 s, every neuron still died at 960 s, and at 3,600 s all neurons died for J50 ≤3,600 s with or without arrest. Arrest barely slows injury accumulation: the J slope fell from 0.996 to 0.886 per s, because energy collapses to ~0 either way. In the injury rule, J grows as dep·(1−e). Extrapolated LT50 for the median neuron ≈ 1.51 h at J50 14,400 s (arrest 0.05: 1.52 h). Reaching the fly LT50 of ~8 h would need J50 ≈ 7.7×10⁴ s (~21 h; LT50 scales ~linearly with J50 here), i.e. intrinsic injury resistance, not energy economy alone (*provisional; extrapolation*). The long-coma runs used speedup 100 after the first 60 s of anoxia.

## 3. Arm B — brain-death map (plan step 1)
*[Figure: Brain-death map — see [docs/figures/dynamics/](../figures/dynamics/)]*

*Fig. 1 — FII_norm (a, b) and dead fraction (c) by anoxia duration × temperature; letters give the modal outcome class (n = 5 seeds per cell).*

**Recoverability boundary (block reset; D_max = duration at which mean FII_norm crosses 0.8, log-interpolated, bootstrap CI):**
| T (°C) | SD onset (s) | D_max (model s) | dead at 240 s | outcome at 480 s | brain death from |
|---|---|---|---|---|---|
| 17 | 15.0 | 359 [342, 381] | 0.5 % | structural loss (16.8 % dead, FII_norm 0.68) | not reached (960 s: 98.7 % dead, survivors in runaway) |
| 22 | 12.5 | 271 [261, 278] | 3.8 % | structural loss (66 % dead) | 960 s (100 % dead) |
| 29 | 8.5 | 149 [145, 153] | 40 % | brain death (99.7 % dead) | 480 s |
| 33 | 6.5 | 124 [112, 129] | 81 % | brain death | 480 s |

* **The brain-death boundary has two layers.** (1) A *functional* boundary at ~7 s of anoxia. Naive reoxygenation at 22 °C ignited a self-sustained runaway attractor in all 30 runs from 7 to 480 s of anoxia (6 durations × 5 seeds); at 960 s everything was dead. In the arm-C naive runs FSA was 0.13–0.49, late firing 105–119 Hz/neuron and FII_norm ≈ 0. (2) A *structural* boundary near 4–5 min at 22 °C. Even with the block reset, identity is lost by neuron death beyond D_max, and FII falls roughly linearly with dead fraction.
* **'Simulated brain death' (silent and unresponsive after the protocol) occurred in 92 runs, and only when ≥97.9 % of neurons were dead.** The model never produced a silent-but-structurally-intact unresponsive state. Every surviving network either responded to probes or ran away. In this model, therefore, irreversible functional death requires near-total structural death. The recoverable window is set by how long it takes the injury integral to kill neurons.
* K⁺: 30–60 s anoxia peaks at 45.2 mM (literature 44–58 mM ✓). For ≥120 s, K⁺ reaches the model ceiling K_max = 80 mM (literature terminal plateau 54 ± 14 mM; the ceiling is a model limitation). Repolarization takes 4.5–8.5 s after reperfusion.
* **Focal ischemia** (right optic lobe, 5 neuropils): SD in the core at 53–56 s, K⁺ peak 45 mM. With the reset, 960 s gave 32.5 % of the whole brain dead but FII_norm 0.997–1.000, because the probe battery barely uses the right optic lobe (see §9, item 5). Right-hemisphere core: 960 s → 43.3 % dead, FII_norm 0.773 (structural loss). Naive reperfusion → runaway in every focal condition.
* **Partial ischemia:** S = 0.25 is tolerated indefinitely (960 s: no SD, 0 dead). S = 0.10 gives SD at 17 s, 0.6 % dead at 240 s and 97.9 % at 960 s (brain death in 4/5 seeds).
* **Hyperthermia alone** (33/37/40 °C for 300 s, S = 1) had no effect (FII_norm 1.00, 0 dead, K⁺ 5.0 mM). The reserve (4 × Q10_supply) outpaces demand Q10 at full supply; heat only acts through faster SD onset and injury during anoxia.

## 4. Arm C — resuscitation protocols (plan step 2)
*[Figure: Example traces — see [docs/figures/dynamics/](../figures/dynamics/)]*

*Fig. 2 — 180-s anoxia, seed 0: naive reperfusion (red) vs a 20-s total block at reperfusion (blue).*

*[Figure: Resuscitation ranking — see [docs/figures/dynamics/](../figures/dynamics/)]*

*Fig. 3 — FII_norm per protocol for 30-, 180- and 360-s anoxia (n = 5 seeds; 95 % bootstrap CI).*

Ranking (mean FII_norm over the three insults; `tables/C_resuscitation_ranking.csv`):
| rank | protocol | 30 s | 180 s | 360 s | dead @360 s | outcome 30/180/360 |
|---|---|---|---|---|---|---|
| 1 | revive_reset | 1.001 | 1.001 | 1.001 | 0.000 | recovered/recovered/recovered |
| 2 | clr3_reset | 1.001 | 0.963 | 0.588 | 0.244 | recovered/recovered/structural_loss |
| 3 | noise_reset | 1.001 | 0.960 | 0.585 | 0.236 | recovered/recovered/structural_loss |
| 4 | adapt_reset | 1.001 | 0.959 | 0.577 | 0.248 | recovered/recovered/structural_loss |
| 5 | block20_del5 | 1.001 | 0.959 | 0.577 | 0.248 | recovered/recovered/structural_loss |
| 6 | staged_fwd | 1.001 | 0.959 | 0.577 | 0.248 | recovered/recovered/structural_loss |
| 7 | block40 | 1.001 | 0.959 | 0.577 | 0.248 | recovered/recovered/structural_loss |
| 8 | block20 | 1.001 | 0.953 | 0.578 | 0.248 | recovered/recovered/structural_loss |
| 9 | staged_rev | 1.001 | 0.953 | 0.578 | 0.248 | recovered/recovered/structural_loss |
| 10 | sdblock_reset | 1.001 | 0.953 | 0.577 | 0.248 | recovered/recovered/structural_loss |
| 11 | block20_del20 | 1.001 | 0.957 | 0.566 | 0.250 | recovered/recovered/structural_loss |
| 12 | block10 | 1.001 | 0.959 | 0.428 | 0.248 | recovered/recovered/runaway |
| 13 | adapt | 0.734 | 0.731 | 0.548 | 0.248 | runaway/runaway/runaway |
| 14 | voff-10 | 0.155 | 0.361 | 0.464 | 0.246 | runaway/runaway/runaway |
| 15 | voff-20 | 0.137 | 0.320 | 0.435 | 0.243 | runaway/runaway/runaway |
| 16 | sdblock | 0.409 | 0.010 | 0.319 | 0.252 | runaway/runaway/runaway |
| 17 | clr3 | 0.718 | -0.170 | 0.020 | 0.244 | runaway/runaway/runaway |
| 18 | ramp10 | 0.402 | 0.161 | -0.157 | 0.267 | runaway/runaway/runaway |
| 19 | gaba4 | 0.183 | 0.349 | -0.173 | 0.239 | runaway/runaway/runaway |
| 20 | noise | 0.260 | 0.015 | 0.008 | 0.240 | runaway/runaway/runaway |
| 21 | naive | -0.058 | 0.084 | 0.085 | 0.252 | runaway/runaway/runaway |
| 22 | gaba2 | 0.191 | 0.100 | -0.195 | 0.239 | runaway/runaway/runaway |
| 23 | block5 | -0.060 | 0.097 | 0.018 | 0.252 | runaway/runaway/runaway |
| 24 | block2 | -0.088 | 0.098 | 0.010 | 0.252 | runaway/runaway/runaway |
| 25 | naive_noinj | -0.079 | -0.078 | 0.105 | 0.250 | runaway/runaway/runaway |
| 26 | pulsed6x5 | -0.104 | -0.032 | -0.123 | 0.281 | runaway/runaway/runaway |
| 27 | pulsed3x10 | -0.057 | -0.093 | -0.136 | 0.273 | runaway/runaway/runaway |
| 28 | ramp60 | -0.104 | -0.065 | -0.212 | 0.318 | runaway/runaway/runaway |
| 29 | ramp30_noinj | -0.305 | 0.064 | -0.175 | 0.281 | runaway/runaway/runaway |
| 30 | ramp30 | -0.301 | 0.063 | -0.246 | 0.292 | runaway/runaway/runaway |
| 31 | cool12 | -0.161 | -0.127 | -0.227 | 0.260 | runaway/runaway/runaway |
| 32 | cool17 | -0.229 | -0.305 | -0.161 | 0.256 | runaway/runaway/runaway |

* **Every protocol without a ≥10-s total block failed (runaway at all three insults).** This includes −10/−20 mV hyperpolarization, 2–4× GABAergic gain, SD blocker (k_rel 0.2), 3× K⁺ clearance, ramped (10/30/60 s) or pulsed (3×10 s, 6×5 s) reperfusion, cooling to 12/17 °C during reperfusion, low-level noise and spike-frequency adaptation (k_adapt 0.2; best of the failures, FII_norm 0.73 at 30/180 s).
* **Minimum effective reset:** 2-s and 5-s blocks failed; a 10-s block worked at 30/180 s but was inconsistent at 360 s (FII_norm 0.428, runaway); 20–40 s worked everywhere.
* **Therapeutic window:** starting the 20-s block 5 or 20 s *after* reperfusion, once runaway is established, still worked (FII_norm 0.841–0.846, the same as immediate block). The runaway attractor is erasable, not imprinting.
* **Remaining loss after the reset is purely structural:** at 360 s, 24.8 % dead → FII_norm 0.578. The counterfactual revive of all dead neurons + reset gave 1.001 at every insult. Add-ons to the block (3× clearance, noise, adaptation, SD blocker, staged release forward or reverse) changed FII_norm by ≤0.01 (all CIs overlap).
* **Staged reawakening** (sensory→central→descending vs reverse, 10 s per stage after a 20-s block): no difference (0.846 vs 0.844).

## 5. Arm D — tolerance and preconditioning (reset protocol)
| condition | SD onset (s) | dead 180 s | FII_norm 180 s | dead 360 s | FII_norm 360 s |
|---|---|---|---|---|---|
| hypomet0.2 | 33.0 | 0.001 | 0.995 | 0.088 | 0.704 [0.573, 0.784] |
| hypomet0.5 | 19.5 | 0.005 | 0.983 | 0.198 | 0.631 [0.576, 0.662] |
| hypomet0.8 | 14.5 | 0.008 | 0.970 | 0.235 | 0.568 [0.452, 0.657] |
| precool12 | 17.0 | 0.000 | 1.001 | 0.006 | 0.977 [0.966, 0.997] |
| precool17 | 15.0 | 0.001 | 0.985 | 0.041 | 0.908 [0.890, 0.938] |
| precond15 | 12.5 | 0.009 | 0.961 | 0.249 | 0.591 [0.442, 0.681] |
| arrest0.05 | 12.5 | 0.001 | 0.983 | 0.059 | 0.711 [0.608, 0.785] |
| unprotected (C block20) | 12.5 | 0.009 | 0.953 | 0.248 | 0.578 |

* **Cooling during the insult was the strongest tolerance intervention:** 12 °C gave 0.6 % dead and FII_norm 0.977 after 360 s. Hypometabolism helped in a graded way (demand 20 %: SD delayed to 33 s; 8.8 % dead, FII_norm 0.704). Coma-coupled arrest at 5 % gave 5.9 % dead and FII_norm 0.711.
* **Preconditioning (15 s anoxia + reset, 60 s earlier) gave no protection** (24.9 % vs 24.8 % dead), as expected for a model with no gene-expression or plasticity mechanisms. Any preconditioning benefit in real tissue must come from mechanisms absent here.

## 6. Arm E — ageing and neurocognitive reserve
*[Figure: Ageing and hub protection — see [docs/figures/dynamics/](../figures/dynamics/)]*

*Fig. 5 — (a) FII_norm vs fraction of synapses lost (n = 3 seeds per point). (b) Targeted neuroprotection (injury_scale 0 on the protected set) after 360-s anoxia.*

* The age parameter a ∈ [0, 1] combines synapse loss of 0.3·a, reserve ×(1 − 0.5a) and clearance ×(1 − 0.5a). FII_norm fell roughly linearly with synapse loss. Demand-weighted loss was more damaging than random loss at every level: 30 % lost → 0.668 vs 0.745; 15 % → 0.803 vs 0.856. Demand-weighted ageing crosses the 0.8 'recovered' line at ~15 % synapse loss; random loss crosses it at ~25 %.
* A 60-s anoxia + reset added no further identity loss at any age; ATP/clearance decline did not matter at this insult size. **SD susceptibility rose only at the oldest level:** partial supply S = 0.3 for 120 s triggered SD (onset 27 s) at a = 1.0 but at no younger age. Full anoxia SD latency changed little (12.5 → 11.5 s).
* Reserve boost, clearance boost and mild hypometabolism (0.8) at a = 1.0 did not restore identity (all 0.668). The age deficit is structural (wiring) and cannot be reached by metabolic support in this model.

## 7. Arm H — targeted (hub) neuroprotection (lead request)
| protected set | n neurons | FII_norm (reset) | dead | FII gain per 1,000 protected | FII_norm (naive) |
|---|---|---|---|---|---|
| top-degree 1% | 1393 | 0.703 [0.668, 0.743] | 0.245 | 0.090 | -0.218 |
| random 1% | 1393 | 0.568 [0.530, 0.591] | 0.246 | -0.007 | 0.236 |
| top-degree 2% | 2785 | 0.744 [0.714, 0.791] | 0.243 | 0.059 | -0.254 |
| random 2% | 2785 | 0.590 [0.559, 0.617] | 0.243 | 0.004 | 0.172 |
| top-degree 5% | 6963 | 0.847 [0.789, 0.905] | 0.235 | 0.039 | -0.253 |
| random 5% | 6963 | 0.574 [0.531, 0.596] | 0.236 | -0.001 | 0.117 |
| top-degree 10% | 13926 | 0.906 [0.845, 0.944] | 0.222 | 0.024 | 0.342 |
| random 10% | 13926 | 0.588 [0.563, 0.637] | 0.224 | 0.001 | -0.118 |
| whole brain | 139255 | 1.001 [1.001, 1.001] | 0.000 | 0.0030 | -0.041 |
| none (C block20 / naive) | 0 | 0.578 | 0.248 | – | 0.085 |

* **Protecting the top 1 % of neurons by degree (1,393 neurons) recovered more than half the identity that whole-brain protection recovers:** FII_norm 0.578 → 0.703 vs 1.001 for whole brain. Per protected neuron this is ~30× more efficient (0.090 vs 0.003 FII per 1,000). The top 10 % gave 0.906 ('recovered'), although 22.2 % of neurons still died. Random sets of the same size gave no gain (0.568–0.590), even though random 10 % prevented a similar number of deaths. Hub protection works by saving the *right* neurons, not more of them.
* **Protection is useless without the block reset:** all naive arms stayed in runaway (FII_norm −0.25 to 0.34).

## 8. Hypothesis tests (preregistration H1, H3, H4, H5; H2/H6 belong to the Structure track)
*[Figure: H4 and H1 — see [docs/figures/dynamics/](../figures/dynamics/)]*

*Fig. 4 — (a) D_max vs temperature with Q10 scaling predictions. (b) FII_norm vs dead fraction for non-runaway reset runs (n = 372).*

* **H1 (sharp percolation-like threshold): FALSIFIED.** Pooled B/C/H reset runs without runaway (n = 372): logistic ΔAIC vs the best alternative (linear) = +0.11 (> −2). 10–90 % width = 1.56 (> 0.25). FII declines approximately linearly with dead fraction. Recovered runs reached up to 0.325 dead (focal optic lobe), while global anoxia lost recovery at ~10 % dead. Under this battery, *where* neurons die matters more than how many.
* **H3 (staged/gradual reperfusion beats abrupt): FALSIFIED (direction reversed, rival H3-R2).** At 360 s, every staged protocol had lower FII_norm than naive: ramp30 Δ = -0.331 [-0.583, -0.130]; ramp60 -0.297 [-0.540, -0.092]; pulsed3x10 -0.221 [-0.460, -0.011]. Staged protocols also had *more* deaths (ramp60 +6.6 pp, ramp30 +4.0 pp; CIs exclude 0), because slower re-energization prolongs depolarization. The peak surge was not lower (ramp30 283 vs naive 265 Hz at 360 s). The only positive cell was ramp10 at 30 s (+0.46, both arms runaway). The negative control was not null (ramp30_noinj − naive_noinj: -0.23 / +0.14 / -0.28), so FII differences between runaway states are attractor noise, not excitotoxicity. The model has no spike-driven injury term, so H3's mechanism is absent by construction. The decisive variable is the block reset (Δ FII_norm vs naive = +1.06 / +0.87 / +0.49, all CIs exclude 0).
* **H4 (hypothermia widens the window by ≈Q10): SUPPORTED.** D_max 17 °C / 29 °C ratio CI = [2.29, 2.57] (excludes 1). Q10_window = 2.09 [1.99, 2.19] at model Q10 2.5; acceptance band [1.67, 3.75]. Sensitivity: model Q10 2.0 → 1.71 [1.60, 1.83] (band [1.33, 3.0]); Q10 3.0 → 2.44 [2.26, 2.64] (band [2.0, 4.5]). All three fall inside the 'approximately Q10' band, but sit consistently at Q10^0.8 (2.5^0.8 = 2.08, 2^0.8 = 1.74, 3^0.8 = 2.41). That is exactly what the Hill-4 hazard predicts (cumulative hazard ∝ rate⁴·D⁵ ⇒ D_max ∝ rate^−4/5). Per the prereg circularity warning, the widening is set by the damage rule, not by network nonlinearity.
* **H5 (identity lost before activity/structure): NOT SUPPORTED as preregistered; mixed.** With the block reset, D50 values coincide: FII_norm 374 s [363, 385], Act 368 s [365, 373], Reach 390 s [388, 391]. So D50_FII ≥ D50_Act, and the 'resuscitated-but-different' region excluding runaway covers 0.0 % of runs ⇒ falsified for reset recovery. With naive reperfusion, D50_FII = 42 s [21, 59] ≪ Reach 388 s, but the CI overlaps Act (291 s [21, 384]). The active + connected + FII_norm ≤0.5 region covers 22.1 % of runs (22.8 % of conditions), all of them runaway. 'Resuscitated but different' exists in this model only as the pathological runaway attractor; once it is cleared, identity tracks structure.
* Multiplicity: H1, H3 and H5 fail on effect-size/structure criteria, not on p-values, so Holm correction does not change any verdict. H4's CI excludes 1 by a wide margin.

## 9. Unexpected observations
1. Reoxygenation ignites the runaway even after 7 s of anoxia, *before* SD at 22 °C, but not at 17 °C (7 s naive: FII_norm 0.988, recovered).
2. Runaway is weaker at high temperature (naive surge 32–54 Hz at 29–33 °C vs 180–263 Hz at 17–22 °C), so naive FII_norm stays at ~0.7 at 29–33 °C, *higher* than at 22 °C. Faster pumps/clearance (Q10 1.8) damp the K⁺ feedback.
3. The runaway is a self-sustained, K⁺-coupled oscillation (Fig. 2: 10–20 mM K⁺ swings, recurrent partial SD waves). A block started up to 20 s after reperfusion still erases it.
4. Silent, unresponsive 'brain death' occurred only with ≥97.9 % of neurons dead. The model has no functional-only death state.
5. Losing one third of the brain (right optic lobe, 32.5 % of neurons) left probe-battery identity intact (FII_norm ≈1.0). FII depends on the battery.
6. Hyperthermia up to 40 °C without hypoxia had no effect.
7. Cooling *during reperfusion* was harmful (runaway, FII_norm −0.30 to −0.13), whereas cooling *before/during the insult* was the best tolerance measure.
8. Targeted hub protection was ~30× more efficient per neuron than whole-brain protection, matching the Structure track's hub-vulnerability result.
9. Time compression (speedup 25–100) does not alter outcome classes.

## 10. Deviations from the task/preregistration
* Single stimulus rate per probe (7 stimuli × 1 rate) instead of 3 rates; 3–5 seeds per condition (prereg: ≥10 seeds per R matrix); baseline 10 seeds.
* FSA under 'weak background drive' was replaced by unstimulated FSA plus 'Act' (evoked-neuron count). The Shiu model is silent at rest, so FSA_baseline = 0 and a ratio is undefined. H5 uses Act instead of FSA.
* Reach is computed on the ≥5-synapse matrix only (no ≥1/≥10 cut-off sensitivity); ageing synapse loss is not reflected in Reach.
* Model choices: J50 raised 60 → 600 s and τ_ATP 20 → 30 s; added q10_injury 2.5; K⁺ hits the 80-mM ceiling for ≥120 s anoxia.
* H1 was fitted on FII vs dead fraction from insult runs, not on targeted/random removal curves (Structure track). H3 protocols are not matched for integrated supply.
* LT50 is extrapolated, not simulated (≤1 h anoxia simulated). 20/26 °C were not run (17/22/29/33 °C were). Hypoglycemia is represented by partial supply only. Keystone rescue was done as revive-all (counterfactual), not per type.
