# Follow-up track 1 — mechanism and robustness of the reperfusion runaway attractor (FAFB v783)

**Status:** main sweep complete, 2,061 condition-seed runs, 0 errors (Modal, 2 × 32-CPU containers, wall 34–37 min). A planned 605-run follow-up job (longer blocks for all LHS sets, k_exc = 10 H3 grid, finer targeted-set brackets, wide hysteresis ramp) was **not run**: it was killed after 3 min because the user's Modal workspace exceeded its spend limit (see §7). All statements describe the `flysim_mech` model (Shiu-2024 LIF on the Princeton ≥5-synapse v783 matrix plus a phenomenological metabolic/ionic layer). The metabolic constants are not fitted to *Drosophila* data. The results are model-mechanistic statements, not evidence about real brain death or human resuscitation.

## 0. Provenance
* **Code** (`followup/mechanism/code/`):
  * `flysim_mech.py` is `dynamics/flysim_dyn.py` plus one optional term (below).
  * `mech_sweep.py` is the driver: arms L/X/N/T/M/Y. It reuses `dynamics/sweep.py` for protocols, probe battery, Reach and `finish`.
  * `flyjob_mech.sh` is the Modal entry script. `drivers_extra.npz` and `done_ids.txt` exist only for the follow-up, which did not run.
  * **Identity check:** with k_exc = 0 the variant reproduces the stored dynamics-track run `C_D30_naive_s0` exactly (FII 0.3617, late rate 122.998 Hz/neuron).
* **Modal jobs:** smoke 900ff527-e669-4ffa-bcd0-5eed11e9297e (84 s); shard A 65780d67-c068-43a9-8cf0-1976a09783a9 (2,229 s); shard B c6fe6ce5-53f5-4e4e-afe0-dc930dc25632 (2,042 s); follow-up 12a56f59-2ec3-4e3c-8934-253f564c9204 (killed after 186 s by the spend limit; results discarded). All handles were closed; `host.compute.ledger()` shows 0 live sandboxes.
* **Raw data:**
  * Modal Volume `flysim-data:/results_mechanism/<arm>/<id>.json|.npz`: full neuron-level spike counts and alive masks. Anatomy runs also carry full epoch × neuron count matrices (`*.heavy.npz`).
  * Local light copies: `followup/mechanism/raw/results/<arm>/` (json + response matrix + traces; anatomy per-neuron arrays), plus `raw/drivers.npz` and `raw/lhs_sets.json`.
  * Result tarballs saved as artifacts: mech_results_A.tar.gz (ff6ccd34-8fe2-465f-abbb-3bbff4fb1adf) and mech_results_B.tar.gz (fe14f8c8-719f-4418-8c9e-86c17c7f4267).
* **Tables** (`followup/mechanism/tables/`):
  * `all_runs_mech.csv`
  * `L_set_summary.csv`, `L_fractions.csv`, `L_prcc.csv`
  * `X_condition_summary.csv`, `X_H3_contrasts.csv`, `X_all_runs_incl_dynamicsC_k0.csv`
  * `T_targeted_summary.csv`
  * `N_driver_scc.csv`, `N_{sc,ty,np}_D{30,180}.csv`
  * `M_dissection_summary.csv`
  * `Y_kick_map.csv`, `Y_kick_threshold.csv`, `Y_gain_kick_map.csv`, `Y_quench_spike_postreset.csv`, `Y_hysteresis.csv`

## 1. Common methods
* **Protocol** (as in the dynamics track): 2 s normoxia → anoxia D s → 60 s reperfusion protocol → 20 s free observation → probe battery. Times are model seconds; speedup 50 means 50 model-s per biological second of spiking. The battery is 7 probes × 1,415 descending/motor readouts.
* **Calibrated set:** princeton5, w_syn 0.185, τ_ATP 30 s, J50 600 s, q10_injury 2.5, 22 °C.
* **FII_norm.**
  * Arms X/T/M: FII is computed against the 10-seed dynamics-track baseline, with ceiling 0.9940.
  * **Arm L:** each parameter set has its own intact baseline (insult none, 4 seeds 100–103, which do not overlap the run seeds). The ceiling is the median leave-one-out r. Per-set ceilings ranged 0.985–0.999 (median 0.993).
* **Classes:**
  * *runaway*: late observation rate > 0.5 Hz/neuron or FSA > 0.05, the dynamics-track rule.
  * *rescued*: FII_norm ≥ 0.8. For arm T, rescue also requires not runaway.
* **Statistics:** bootstrap 95 % CIs (2,000–4,000 resamples) over seeds, or over parameter sets for the L fractions. Paired-by-seed bootstrap for protocol contrasts. PRCC on ranks with Holm correction over all 250 (set-level outcome × parameter) tests.

## 2. Robustness of the runaway attractor and the block reset (arm L)
![Robustness]({{artifact:3ddf2995-f748-4e7c-a839-351fe1f8af0f}})

*Fig. M1. (a) Set-mean FII_norm, 48 Latin-hypercube parameter sets × 3 seeds; the diamond is the calibrated set; dashed line = 0.8. (b) PRCC of set-level outcomes with each parameter; * = Holm p < 0.05. (c) Outcome of the 20-s block after 180-s anoxia in the (τ_ATP, reserve) plane.*

The design was 48 LHS sets plus the calibrated set, varying 10 uncertain constants:

| parameter | range |
|---|---|
| w_syn | 0.15–0.25 |
| τ_ATP | 15–60 s, log scale |
| J50 | 300–1,200 s, log scale |
| g_K | 0.3–0.9 |
| k_rel | 20–120, log scale |
| k_clr | 0.5–2, log scale |
| reserve | 2–8, log scale |
| demand β | 0.25–1 |
| K_max | 50–100 mM |
| glutamate | inhibitory or excitatory |

Each set ran 2 anoxia durations × 2 protocols × 3 seeds, plus its baseline: 784 runs.

| anoxia | protocol | sets with runaway (majority of seeds) | sets rescued (FII_norm ≥ 0.8) | mean FII_norm [95 % CI] |
|---|---|---|---|---|
| 30 s | abrupt | **48/48 = 1.00 [1.00, 1.00]** | 4/48 = 0.08 [0.02, 0.17] | 0.24 [0.14, 0.34] |
| 30 s | abrupt + 20-s block | 0/48 | **48/48 = 1.00 [1.00, 1.00]** | 1.001 [1.000, 1.001] |
| 180 s | abrupt | **48/48 = 1.00 [1.00, 1.00]** | 2/48 = 0.04 [0.00, 0.10] | 0.17 [0.08, 0.27] |
| 180 s | abrupt + 20-s block | 11/48 = 0.23 [0.13, 0.35] | **34/48 = 0.71 [0.58, 0.83]** | 0.75 [0.63, 0.85] |

* **(a) Runaway after abrupt reperfusion is a property of the whole explored parameter box.** It occurred in 288/288 runs. Its intensity varies: late rate 2–250 Hz/neuron, median 92–95. In 41/48 sets (0.85 [0.75, 0.96]) it is the high-rate (≥10 Hz/neuron) K⁺-coupled state. In 7 sets it is a low-rate (2–5 Hz/neuron) reverberation; these sets have w_syn ≤ 0.19 and glutamate inhibitory, and 4 of them keep FII_norm 0.81–0.86 despite being classed runaway.
* **(b) The 20-s block always rescues after 30 s of anoxia, but after 180 s only in 71 % of sets.** The 14 failing sets fail in two ways:
  * *Re-ignition* (11 sets, FII_norm ≈ 0; late rate 4–216 Hz/neuron). These sets have slow ATP recovery (median τ_ATP 44 s vs 26 s in rescued sets) and low reserve (2.6 vs 4.7). Energy has not recovered when the block ends, so the network falls back into the attractor. PRCC with re-ignition: τ_ATP +0.55 (Holm p = 0.02); reserve −0.48 (unadjusted p = 0.002, Holm 0.16). Two of the 11 (sets 7 and 14) re-ignite only into the low-rate state.
  * *Structural loss* (3 sets, J50 ≈ 325 s, ~10 % dead).
* **Dominant determinants.** Dead fraction after reset at 180 s: J50 PRCC −1.00 and τ_ATP −0.75 (both Holm p < 0.001). Abrupt-reperfusion runaway intensity: glutamate sign PRCC −0.71 to −0.73 (Holm p < 0.001; excitatory glutamate gives median late 156 vs 50 Hz/neuron at 180 s) and τ_ATP +0.55. w_syn, g_K, k_rel, k_clr and K_max had no Holm-significant effect on any outcome.
* **Interpretation.** The existence of the runaway attractor is robust. The adequacy of a *fixed* 20-s reset is not: a block shorter than energy recovery re-ignites. Longer blocks were planned for exactly this test, but that follow-up did not run (§7). The dynamics track, at the calibrated set, found 40-s blocks no better than 20-s ones, but it did not test the failing sets.

## 3. Excitotoxic reperfusion term and H3 (arm X)
![Excitotoxicity and H3]({{artifact:02319029-3e70-4e56-8499-72b4b8db57ea}})

*Fig. M2. (a) Dead fraction after 180-s anoxia vs k_exc (n = 5 seeds per point). (b) Paired Δ dead, graded − abrupt, no block. (c) Paired Δ FII_norm, graded + block − abrupt + matched block. Dots are means with paired-bootstrap 95 % CIs; k_exc = 0 cells for protocols without a block come from the dynamics-track arm C (identical code path).*

* **Term:** dJ_i/dt += q·k_exc · max(r_i − r0, 0)/r0 · Kex_i/(Kex_i + K0), with r0 = 50 Hz.
  * Kex_i is the neuron-mixed extracellular K⁺ excess, a Ca²⁺-entry/glutamate proxy. It vanishes at rest and during probe responses.
  * It accrues only while neurons fire above 50 Hz in high K⁺, i.e. during the reperfusion surge and the runaway. Repair and the Hill-4 hazard are unchanged.
* **Calibration** (180-s anoxia, 5 seeds; abrupt vs abrupt + 20-s block):

  | k_exc | dead, abrupt | dead, + block |
  |---|---|---|
  | 0 | 0.9 % | 0.9 % |
  | 0.3 | 0.9 % | 0.9 % |
  | 1 | 1.1 % | 0.9 % |
  | 3 | 3.5 % | 0.9 % |
  | 10 | 39.5 % | 1.0 % |

  * Mean excitotoxic J under abrupt reperfusion was 22 / 72 / 218 / 717 at k_exc 0.3 / 1 / 3 / 10, against ≤2.3 with the block.
  * k_exc = 3 is the smallest tested value at which abrupt reperfusion clearly adds injury (+2.6 pp). H3 was tested at k_exc 1 and 3 with 12 protocols × 3 durations × 5 seeds; k_exc = 10 has calibration cells only.
* **Without a block** (every run of every graded protocol stayed in runaway, so FII differences are attractor noise):
  * At k_exc 0 and 1, graded reperfusion added death at every duration: +0.1 to +0.6 pp at 180 s and +1.5 to +6.6 pp at 360 s.
  * At **k_exc = 3** the direction partly flipped. Pulsed 3×10 s reduced death by −1.6 pp [−1.8, −1.4] at 180 s and −0.8 pp [−1.2, −0.4] at 360 s. Ramp 60 s reduced it by −0.8 pp [−0.9, −0.6] at 180 s. Ramps of 10–60 s still added +1.7 to +3.1 pp at 360 s.
* **With a block** (the clinically relevant comparison):
  * No graded protocol beat abrupt reperfusion plus a block, at any k_exc or duration. Δ FII_norm ranged −0.045 to +0.006 for the late-block combinations, and every combination added +0.1 to +7.6 pp dead.
  * *Ramp 30 s + block at 0 s* was much worse: Δ FII_norm −0.15 to −0.90, with 100 % re-ignition, because the block ends while supply is only at 67 %.
  * The block's advantage over abrupt reperfusion grew with k_exc: Δ dead −0.2 pp at k_exc 1, −2.6 pp at 3 and −38.6 pp at 10 (180 s).
* **Verdict.** With the excitotoxic term, the H3 reversal persists for identity (FII) and for every comparison that includes a reset. A weak H3-direction effect on death appears without a reset, for pulsed and slow-ramp reperfusion at k_exc = 3; it is small (≤1.6 pp) and the network stays in runaway. In this model the decisive variable remains the post-reperfusion activity block. The block is also the best protection against excitotoxic injury, because it removes the surge that injures. *Provisional:* k_exc = 10 was tested only for abrupt vs block; the full k_exc = 10 H3 grid did not run (§7).

## 4. Attractor anatomy and targeted reset (arms N, T)
![Anatomy and targeted reset]({{artifact:56ded7bf-fc03-4491-a98a-62e77f9e176e}})

*Fig. M3. (a) Population rate after reoxygenation (5 seeds thin, mean thick). (b, c) FII_norm on held-out seeds 10–14 (n = 5, bootstrap 95 % CI) vs size of the set blocked for 20 s from reperfusion onset; black = whole-brain block, dotted = no block. (d) Fraction of each set in its largest strongly connected component (SCC) of the excitatory subgraph.*

**Ignition** (abrupt reperfusion, seeds 0–4; 5-ms epochs over the first 15 model-s, then 50–100-ms epochs):
* After **30 s** of anoxia, ignition is near-synchronous: the population jumps from ~0 to ~100 Hz/neuron within ≈1.5 model-s, starting ≈4 model-s after reoxygenation, then relaxes to ~7 Hz before the self-sustained state builds.
* After **180 s**, ignition starts later (≈6.5 model-s) and climbs over ~10 model-s.
* Early firers are central and descending neurons in the antennal lobe (local neurons lLN1_bc, lLN2P/X12; DA1 lPNs) and AVLP/ATL (median first spike 80–85 ms of biological time at 30 s). Optic neurons (lobula plate) and sensory neurons fire last (105–110 ms).
* In the sustained state, 73 % (30 s) and 86 % (180 s) of neurons fire ≥20 Hz (mean 97 and 121 Hz/neuron).
* Only 13–20 % of neurons are *intrinsically* suprathreshold (effective resting potential above effective threshold) at the end of reperfusion, so most of the firing is recruited.

**Driver sets.** These were defined before any T run, from the N runs (pooled 30 + 180 s, seeds 0–4), excitatory neurons only:
* *sustain:* rate × excitatory synapses onto neurons active ≥20 Hz.
* *ignite:* earliest mean first spike.
* The top 0.5 % sustain set (696 neurons) is 74 % central brain and enriched in visual-centrifugal (22×), ascending and central neurons. It is concentrated in AL_L/AL_R (279 neurons; lLN1_bc, DA1/VA1v/VM5d/DL2v PNs) and includes hΔB.
* **96 % of this set forms a single excitatory SCC (31 excitatory edges per node).** The figures for comparison sets are:

  | set (same size) | in largest SCC | excitatory edges per node |
  |---|---|---|
  | sustain drivers | 96 % | 31 |
  | top-degree hubs | 40 % | 9.4 |
  | random | 0.4 % | 0.1 |

  The ignite set reaches 81–95 %.
* The Modal-built rankings were verified locally (identical top 700).

**Targeted 20-s block on held-out seeds:**

| anoxia | sustain drivers | earliest igniters | top-degree hubs | random | whole brain |
|---|---|---|---|---|---|
| 30 s: smallest rescuing set | **0.5 % (696 neurons), 5/5 seeds, FII_norm 0.947 [0.943, 0.951]** | **0.5 %, 5/5, 0.945** | 2 % (2,785), 5/5 | none up to 20 % (best 0.33) | 1.000 |
| 180 s: smallest rescuing set | **20 % (27,851), 5/5, 0.979** (10 %: 0/5) | **20 %, 5/5** (10 %: 0/5) | none up to 20 % (0/5) | none up to 20 % (0/5) | 0.979 |

* After 30-s anoxia, blocking ≤0.1 % of sustain drivers did not rescue but cut the late runaway from 112 to ~17 Hz/neuron. Rescue appears between 0.1 % and 0.5 %; a finer bracket did not run.
* **Driver block is ≥4× more efficient than hub block and infinitely more efficient than random block at 30 s.** At 180 s the minimal set is 20 % of the brain (bracket 10–20 %, not refined; §7). This is consistent with the K⁺ mechanism in §5: after long anoxia K⁺ is near its ceiling everywhere, so a larger fraction of neurons must be kept from firing while the pumps recover.
* *Model artefact noted:* 671 neurons without neuropil synapse assignments receive no K⁺ coupling. 109 of them enter the ignite ranking, because at t = 0 they fire as the global supply returns.

## 5. Mechanism dissection and bistability (arms M, Y)
![Bistability and dissection]({{artifact:858aee6e-225e-4029-b5dc-94e1b59b6c43}})

*Fig. M4. (a) Intact network: state 30–40 s after a 2-s K⁺ clamp kick, as a function of tonic offset δ (5 seeds at δ = 0, 3 otherwise); blue line = analytic intrinsic threshold K_th = K0·exp((V_th − V_rest − δ)/(g_K·RT/F)). (b) Slow ramp of δ, −6 → +10 → −6 mV over 2 × 160 s (3 seeds; thick = seed 0). (c) Manipulation applied 20–60 s after abrupt reperfusion (n = 5). (d) Quench of a kicked runaway by a total block of duration B (n = 3).*

**Dissection** (arm M; rate in the last 10 s of a 40-s manipulation, Hz/neuron, 30 s / 180 s anoxia; unmanipulated runaway 105 / 121):

| manipulation | rate during manipulation | late rate after release |
|---|---|---|
| all synapses off | 109 / 134 | 138 / 120 |
| excitatory synapses off | 14 / 120 | 38 / 91 |
| inhibitory synapses off | 133 / 135 | 111 / 115 |
| **K⁺ clamped at 5 mM** | **2.4 / 2.3** | 4.3 / 3.9 |
| **no spike-driven K⁺ release (k_spk = 0)** | **2.4 / 2.3** | 4.3 / 4.1 |
| no SD K⁺ release (k_rel = 0) | 56 / 153 | 129 / 25 |

* **The high-rate attractor is ionic, not synaptic.** It survives removal of all synaptic transmission. It collapses within seconds when activity-dependent extracellular K⁺ accumulation is removed.
* Synaptic excitation contributes after short anoxia: removing excitatory synapses cuts 30-s runaway from ~56 to 14 Hz.
* Releasing a K⁺-clamp leaves a residual low-rate (≈4 Hz/neuron) synaptic reverberation with FII_norm 0.71–0.72. This is a third state, distinct from both quiet and runaway.

**Bistability in the intact network** (arm Y, no insult):
* **Coexistence.** With identical parameters (δ = 0), the unperturbed network stays silent (5/5 seeds; K⁺ 5.0 mM), while a 2-s K⁺ kick to ≥9 mM drives it into a self-sustained state (5/5 seeds; 65–129 Hz/neuron at 9–12 mM; K⁺ ≈ 15 mM). The same holds after a recovered 30-s anoxia + reset: no kick → quiet (3/3); 10- or 15-mM kick → runaway (3/3 each).
* **Basin boundary.** The minimal kick that reaches the high-rate state (≥10 Hz/neuron):

  | offset δ | +6 mV | +4 mV | +2 mV | 0 mV | −2 mV | −4 mV | −6 mV |
  |---|---|---|---|---|---|---|---|
  | minimal kick | 6 mM | 7 mM | 8 mM | 9 mM | 12 mM | 15 mM | none up to 30 mM |
  | analytic K_th | 5.3 | 6.1 | 6.9 | 7.9 | 9.0 | 10.2 | 11.7 mM |

  * At δ ≥ 0 the simulated boundary sits one grid step above the analytic value. It diverges at hyperpolarized offsets, where kicks reach only a ~1–2 Hz low-rate state.
  * Synaptic gain did not move the boundary: at 0.5, 0.75, 1.25 and 1.5× excitatory gain, 8 mM failed and 10 mM succeeded.
  * Spike kicks (100 Hz Poisson to 1–50 % of neurons for 2 model-s) produced only the low-rate state (5.6–19 Hz/neuron, K⁺ 5.6–8.2 mM).
* **Hysteresis.**
  * Ramping δ up with no noise, the silent network ignites at δ = +7.15 mV in 3/3 seeds. This equals V_th − V_rest = 7 mV, where resting neurons become intrinsically suprathreshold. With weak noise (2 % of neurons at 2 Hz) it ignites at +2.2 to +3.7 mV (>20 Hz criterion).
  * Ramping back down, the active state never extinguishes within the range: it is still at 252–276 Hz/neuron (no noise) or 137–156 Hz/neuron (noise) at δ = −6 mV.
  * Along the down-ramp, energy falls (e 0.94 → 0.73) and mean K⁺ rises (9.5 → 20 mM). The hyperpolarized branch is therefore sustained by an ionic–metabolic loop: firing → K⁺ accumulation and ATP demand → pump failure → depolarization.
  * The loop is at least 13 mV wide (+7 → < −6 mV). Its lower edge was not found; the planned −30 mV ramp did not run.
* **Escape.** A total block quenches a kicked runaway only if it lasts ≥10 model-s at δ = 0 (1, 2 and 5 s fail, 3/3 seeds each), or ≥20 s at δ = +4 mV (10 s fails). This matches the dynamics-track finding that 2- and 5-s blocks fail and 10–20-s blocks work: the block must outlast K⁺ clearance (k_clr = 1 s⁻¹) and ATP recovery before release.

**Dynamical-systems interpretation (model-internal).**
* The slow variables (compartment K⁺, ATP, release store) and fast spiking form a positive feedback loop: spiking → K⁺ release (k_spk) → Nernst depolarization (g_K·RT/F·ln K/K0) → spiking.
* This creates two stable fixed points, quiet (K⁺ = 5 mM) and active (K⁺ ≈ 10–20 mM, oscillating between firing and depolarization block), separated by a K⁺ threshold set by the single-neuron gap V_th − V_rest − δ. The system shows a hysteretic, discontinuous transition: a saddle-node-like ignition at δ ≈ +7 mV and a lower turning point below −6 mV.
* Anoxia does not create the attractor. It provides the perturbation: the extracellular K⁺ load (45–80 mM), which on reoxygenation passes through the window between the K_th (~8 mM) threshold and depolarization block (~−32 mV), places the network inside the active basin.
* The block reset works by holding the system silent until K⁺ and ATP relax back into the quiet basin. Targeted blocks work when they remove enough recurrent excitatory K⁺ sources to keep the compartments below K_th during that relaxation.

## 6. Hypothesis and claim verdicts
* **Runaway attractor is robust to parameter uncertainty: SUPPORTED.** It occurred in 48/48 LHS sets and 288/288 runs; it was high-rate in 41/48.
* **20-s block reset generalizes: SUPPORTED at 30 s (48/48), CONDITIONAL at 180 s (34/48).** Failures are re-ignition (slow ATP recovery, low reserve) or structural loss (low J50).
* **H3 (graded beats abrupt reperfusion) with an excitotoxic term: still FALSIFIED for FII and for every block combination.** There is a weak, non-decisive H3-direction effect on death without a block for pulsed / 60-s ramp at k_exc = 3 (−0.8 to −1.6 pp). The H3 reversal from the dynamics track is not an artefact of the missing excitotoxicity term, within k_exc ≤ 3. k_exc = 10 is provisional (abrupt vs block only).
* **Attractor has an identifiable driver core: SUPPORTED.** The top 0.5 % sustain/ignite set is a single strongly connected excitatory core (81–96 %, vs 40 % hubs and 0.4 % random), centred on antennal-lobe local and projection neurons.
* **Targeted reset can replace whole-brain block: SUPPORTED after 30 s** (0.5 % driver block rescues 5/5 held-out seeds; hubs need 2 %, random fails at ≤20 %). **Largely NOT SUPPORTED after 180 s** (20 % needed; hubs and random fail at ≤20 %).
* **Bistability and hysteresis in the intact network: SUPPORTED.** Quiet and active states coexist at the same parameters. The basin boundary follows the analytic K⁺ threshold at δ ≥ 0. Hysteresis is ≥13 mV wide. The escape time is ≥10–20 s. The mechanism is ionic (K⁺) and survives removal of all synapses.

## 7. Deviations
* **Follow-up runs not executed** (Modal spend limit reached at ~11:30Z; 605 runs). Affected items:
  * 40/60-s blocks for all 49 LHS sets at 180 s.
  * The full k_exc = 10 H3 grid.
  * Targeted-set brackets: 12.5–17.5 % at 180 s; 0.2–0.3 % at 30 s; hubs at 30–50 % and 1.5 %; rankings from 180-s runs only.
  * A wide hysteresis ramp to −30 mV.

  The corresponding statements are labelled as brackets or provisional.
* **LHS:** 48 sets + calibrated set × 3 seeds, as requested; D = 30 and 180 s at 22 °C only. The planned 'demand beta' parameter was varied; q10s and temperatures were not.
* **Per-set FII baselines** use 4 seeds and a leave-one-out ceiling (dynamics track: 10 seeds). Arm X k_exc = 0 cells for naive, ramp10/30/60, pulsed3x10, block20 and block20_del20 were reused from dynamics arm C, under the same seeds and an identical code path (verified on one run), not rerun.
* **Excitotoxic term:** one functional form (rate × K⁺-excess proxy, r0 = 50 Hz) and k_exc ∈ {0.3, 1, 3, 10}. No literature calibration exists for this proxy; the 'adds injury' criterion is the model's own.
* **Anatomy:** first-spike resolution is 5 ms of biological time (0.25 model-s). Driver rankings were built from 10 runs. Inhibitory neurons are excluded from driver sets (silencing them would raise excitation). Only the 20-s block from reperfusion onset was tested for targeted sets.
* **Arm Y:**
  * δ is a tonic depolarizing/hyperpolarizing offset standing in for "K⁺ offset". Exact K⁺ baseline shifts were not run.
  * Kicks are 2-s K⁺ clamps or Poisson spike drives.
  * The 'runaway' rule (>0.5 Hz/neuron) also captures the low-rate state, so the high-rate state is reported separately (≥10 Hz/neuron).
  * Y runs have no probe battery (no FII).
* **Statistics:** n = 3 seeds per Y cell (5 at δ = 0) and per LHS cell. Several between-seed CIs are degenerate (zero width) because outcomes are nearly deterministic given the protocol.

## 8. Limitations
* All conclusions are internal to a LIF + phenomenological K⁺/ATP model with unfitted constants. The ionic attractor depends on g_K·RT/F coupling, compartment-mixed K⁺ and the model's depolarization-block rule, none of which is validated against fly tissue.
* Real neurons have K⁺-dependent conductances, glial buffering, Na⁺/Ca²⁺ dynamics and neuromodulators that could abolish or create such bistability.
* FII depends on a 7-probe battery read out at descending/motor neurons. The low-rate state can keep FII_norm ≈ 0.7–0.86 while being 'runaway'.
* The excitotoxic proxy uses extracellular K⁺, not Ca²⁺ or glutamate. Results at k_exc = 10 are provisional.
* Driver sets were identified and tested in the same model at the calibrated set, on held-out seeds only. Their generality across parameter sets and connectomes (MCNS/BANC) is untested. Antennal-lobe enrichment may reflect local recurrent excitation that is specific to this sign assignment (glutamate inhibitory, unknown NT excitatory).
* Speedup 50 compresses slow dynamics. The dynamics track showed outcome classes are invariant for speedup 25–100, but escape times in model-s and biological-s differ by that factor.
