# Supplementary Information

**A whole-brain *Drosophila* connectome model predicts a duration- and sign-dependent potassium trap after anoxia and a reset-based route to functional recovery**

Srihari Muralikrishnan (Kalaris Labs). Preprint v1.0, 2026-10-20. Contact: hello@kalarislabs.com

All tables below were extracted programmatically from the analysis result files (`research/structure/STRUCTURE_RESULTS.md`, `research/dynamics/DYNAMICS_RESULTS.md`, `research/followup/mechanism/MECHANISM_RESULTS.md`, `research/followup/replication/REPLICATION_RESULTS.md`, `research/sim/README.md`, `research/lit/parameters.csv`, `research/followup/replication/tables/v4_condition_summary.csv`, `research/followup/mechanism/raw/sign_res/`, `research/followup/index/tables/`) by `build_paper.py`. The claim-by-claim synthesis is `research/SYNTHESIS.md`. They are not retyped. Every value describes the model, not real tissue.

---

## S1. Preregistration (registered 2026-10-16, before any insult simulation)

The full text is in `research/lit/preregistration.md`. Summary of hypotheses and falsification criteria:

| H | Key observable | Supported if | Falsified if |
|---|---|---|---|
| H1 sharp threshold | FII(f), Reach(f) | logistic preferred by ΔAIC ≤ −10, 10–90% width ≤ 0.10, irrecoverable beyond f_c | width > 0.25, or ΔAIC > −2, or recovery at f ≥ f_c + 0.10 |
| H2 targeted vulnerability | ρ(demand, death time); f_c,targeted/f_c,random | ρ ≤ −0.3 and ratio ≤ 0.7 (CI excludes 1) | ρ > −0.1 and ratio CI includes 1 |
| H3 graded reperfusion | ΔFII staged − abrupt | ≥ +0.05 (CI excludes 0) in ≥2/3 severities, null negative control | ΔFII ≤ 0, or benefit survives without excitotoxic term |
| H4 hypothermia | Q10_window | within [0.67, 1.5] × model Q10 | D_max ratio CI includes 1, or outside band |
| H5 identity first | D50 ordering; regime coverage | D50_FII lowest by ≥0.05; regime ≥5% of grid | D50_FII not lowest and coverage < 1% |
| H6 replication | outcome concordance | ≥3/4 connectomes match; f_c max/min ≤ 1.5 | ≤2/4 match, or f_c max/min > 2 |

Operational "brain death" (preregistered): <1% of neurons spontaneously active for ≥1 s under baseline drive, no readout response (< 2 SD of noise), persisting after the protocol ends. As run (deviation S6): spontaneous fraction < 0.01 and no readout response > 8 Hz after the protocol. This is a model analogue of the clinical construct [Greer et al. 2023], not a clinical determination.

---

## S2. Simulator constants

### Table S2a. Fast LIF layer and calibration (from `research/sim/README.md`)

| constant | value | source |
|---|---|---|
| V_rest = V_reset | −52 mV | Shiu et al. 2024 |
| V_th | −45 mV | Shiu et al. 2024 |
| τ_m / τ_syn | 20 ms / 5 ms | Shiu et al. 2024 |
| refractory / delay | 2.2 ms / 1.8 ms | Shiu et al. 2024 |
| w_syn | 0.275 mV published; **0.185 mV used** (Princeton ≥5 matrix supercritical at 0.275) | calibrated here |
| Poisson kick | 68.75 mV to v | Shiu et al. 2024 |
| dt | 0.1 ms (fast), 5 ms (slow layer) | — |

{{T:sim_readout}}

### Table S2b. Slow metabolic/ionic layer (phenomenological, **not fitted to *Drosophila* data**)

| constant | default | calibrated set (dynamics, mechanism) | LHS range (mechanism arm L) |
|---|---|---|---|
| τ_ATP | 20 s | **30 s** | 15–60 s (log) |
| ATP reserve m | 4 × resting demand | 4 | 2–8 (log) |
| demand exponent β (area^β) | 0.5 | 0.5 | 0.25–1 |
| r_half (rate demand) | 100 Hz | 100 Hz | — |
| pump failure e_ok | 0.75 | 0.75 | — |
| K₀ / K_max | 5 / 80 mM | 5 / 80 mM | K_max 50–100 mM |
| g_K (Nernst coupling) | 0.6 | 0.6 | 0.3–0.9 |
| k_rel (SD release) | 60 mM/s | 60 | 20–120 (log) |
| k_clr (clearance × pump) | 1 s⁻¹ | 1 | 0.5–2 (log) |
| V_sd | −30 ± 2 mV | same | — |
| depolarization block | h → 1 − σ((v_rest,eff + 32)/2), τ_h 10 ms | same | — |
| injury hazard | h_max·J⁴/(J⁴ + J50⁴), h_max 1/30 s⁻¹ | — | — |
| J50 | 60 s | **600 s** | 300–1,200 s (log) |
| Q10 demand / supply / pump / injury | 2.5 / 1.5 / 1.8 / — | 2.5 / 1.5 / 1.8 / **2.5** | not varied |
| speedup (slow-time compression) | — | **50** | not varied (25–100 tested in arm A) |
| glutamate sign | inhibitory | inhibitory | inhibitory or excitatory |
| w_syn | 0.185 | 0.185 | 0.15–0.25 |

### Table S2c. Calibration runs (dynamics arm A)

{{T:dyn_armA}}

---

## S3. Literature parameter table (Table S1; `research/lit/parameters.csv`, 81 rows, 6 flagged UNVERIFIED)

Values are as extracted from the cited sources. "UNVERIFIED" rows were not used as calibration targets.

{{PARAMS}}

---

## S4. Structural percolation (all four connectomes)

### Table S4a. Data and pairing QC
{{T:st_qc}}

### Table S4b. Baseline reachability (intact, strict)
{{T:st_base}}

### Table S4c. Damage curves, primary metric Reach_N(L≤6)
{{T:st_main}}

### Table S4d. D50 (f at 50% of baseline) per structural metric
{{T:st_d50}}

**Controls.** Degree-preserving rewired null: FAFB random f_c 0.893 (real 0.800), hub-targeted 0.266 (real 0.035); MANC random 1.107 (extrapolated; real 1.027), hub-targeted 0.613 (real 0.388). Edge cut-off ≥10: FAFB random f_c 0.607, w 0.960; hub-targeted f_c 0.013, w 0.087, ΔAIC 0.2 (falsified). MANC ≥1: random f_c 0.787, w 0.124 (indeterminate); hub-targeted f_c 0.770. The H1 category for the hub-targeted regime (FAFB) and the random regime (MANC) is therefore not robust to the edge cut-off. Negative control (off-path neurons removed): Δ = 0 for every metric in every dataset (off-path fractions FAFB 2.1%, BANC 27.0%, MCNS 2.4%, MANC 3.2%).

**Keystones.** Max relative ΔReach_N(L6) for single-type removal: FAFB 0.88% (CB1231; 154 types evaluated), BANC 1.03% (IN13B010; 218), MCNS 0.86% (AN04A001; 194), MANC 0.00% (161). Top-50 membership overlap did not exceed chance once the universe was restricted to types evaluated in both datasets (p = 0.129–0.634). Graded importance replicated (Spearman ρ of influence loss: FAFB–BANC 0.730, FAFB–MCNS 0.710, BANC–MCNS 0.911, BANC–MANC 0.652, MCNS–MANC 0.776).

![Figure S1](figures/FigS1_structure_metrics_controls.png)

**Figure S1.** (a) FAFB random loss, all structural metrics; (b, c) degree-preserving null and edge cut-off controls for FAFB and MANC.

![Figure S2](figures/FigS2_keystone_types.png)

**Figure S2.** (a) Relative reach lost when all neurons of one cell type are removed (each dot a type; bar = median); (b–d) cross-dataset agreement of relative weighted-influence loss for shared types (Spearman ρ).

---

## S5. Dynamics: brain-death map, protocols, tolerance, neuroprotection

### Table S5a. Recoverability boundary (20-s block reset)
{{T:dyn_boundary}}

### Table S5b. All 32 reperfusion protocols, ranked by mean FII_norm (30/180/360 s anoxia, 5 seeds)
{{T:dyn_rank}}

Protocol key: `naive` = abrupt step to S = 1; `blockN` = total spike + synaptic block for N s from reperfusion; `block20_delN` = block starting N s after reperfusion; `rampN` = linear supply ramp over N s; `pulsedAxB` = A cycles of B-s reperfusion; `voff-N` = −N mV offset; `gabaN` = N× GABAergic gain; `sdblock` = k_rel × 0.2; `clr3` = 3× clearance; `noise` = low-level background Poisson; `adapt` = spike-frequency adaptation k_adapt 0.2; `coolN` = cooling to N °C during reperfusion; `staged_fwd/rev` = staged release sensory → central → descending (or reverse), 10 s per stage after a 20-s block; `revive_reset` = counterfactual revival of dead neurons + reset; `_noinj` = injury off; `_reset` suffix = followed by 20-s block.

### Table S5c. Tolerance and preconditioning (arm D; reset protocol)
{{T:dyn_armD}}

### Table S5d. Targeted (hub) neuroprotection (arm H; 360-s anoxia)
{{T:dyn_armH}}

---

## S6. Mechanism follow-up

### Table S6a. Latin-hypercube ranges (48 sets + calibrated set)
{{T:m_lhs}}

### Table S6b. Runaway and rescue across parameter sets
{{T:m_lhs_out}}

### Table S6c. Excitotoxic term calibration (180-s anoxia)
{{T:m_kexc}}

### Table S6d. Driver-set connectivity
{{T:m_scc}}

### Table S6e. Targeted 20-s block on held-out seeds
{{T:m_target}}

### Table S6f. Dissection (rate in last 10 s of a 40-s manipulation, Hz/neuron, 30 s / 180 s anoxia)
{{T:m_dissect}}

### Table S6g. Basin boundary: minimal K⁺ kick reaching the high-rate state vs tonic offset δ
{{T:m_basin}}

![Figure S3](figures/FigS3_excitotoxic_H3.png)

**Figure S3.** (a) Dead fraction after 180-s anoxia vs excitotoxic gain k_exc (n = 5 seeds per point). (b) Paired Δ dead, graded − abrupt, no block. (c) Paired Δ FII_norm, graded + block − abrupt + matched block. Dots = means with paired-bootstrap 95% CIs.

### S6h. Transmitter-sign robustness (antennal-lobe local neurons inhibitory)

Variant A sets every neuron with class ALLN or primary type matching lLN/vLN/v2LN/ilLN/il2LN to inhibitory. In FAFB the mask covers 431 neurons; 284 were already inhibitory and 147 changed sign (90 predicted cholinergic, 57 predicted monoaminergic). The same rule changed 147 neurons in the FAFB-proxy model and 219 in MCNS.

**Table S6h-1. FAFB full model, variant A, 30-s anoxia, naive reperfusion (per seed; local runs).** FII is raw (vs the default-sign intact baseline; no variant-A baseline was run); FII_norm ≈ FII/0.9940. Runaway rule: late rate > 0.5 Hz/neuron or FSA > 0.05. Default-sign comparison (dynamics arm C, 5 seeds): FII_norm −0.058, runaway 5/5; anatomy runs N_anat_D30 s0–s4: mean observation rate 91.3–104.3 Hz/neuron, 64.0–69.5 % of neurons ≥ 20 Hz.

{{SIGNA}}

**Table S6h-2. MCNS and FAFB-proxy, default signs vs antennal-lobe local neurons inhibitory (`_alln`), 5 seeds per condition.**

{{SIGNREP}}

Variants B (all low-confidence calls inhibitory) and C (neurons without a call silenced) were prepared in code but not run. FAFB variant A was not run at 180 s, with the reset, for bistability kicks, or for targeted blocks.

---

## S7. Cross-connectome dynamic replication (MCNS v1.0; BANC v888 not analysed)

Source: `research/followup/replication/REPLICATION_RESULTS.md` (211 analysed runs: MCNS 146, FAFB-proxy 65). MCNS: 166,700 neurons, 6,242,085 edges (≥5 synapses), w_syn 0.0538 mV (same distance below the reverberation threshold as FAFB). Proxy substitutions: demand = total synapse count; K⁺ coupling = connection-based Jaccard 9-nearest-neighbour neuropil graph. A dense shared-neuron K⁺ adjacency was numerically unstable (spurious oscillation and runaway in sham runs) and was replaced; the FAFB dynamics track used the stable geometric adjacency and is unaffected. Responsive probes only (sugar → MN9 is silent in MCNS); noise ceilings MCNS 0.9972, FAFB-proxy 0.9922.

### Table S7a. Condition summary
{{T:rep_cond}}

### Table S7b. Hub-protection efficiency in MCNS (360-s anoxia + 20-s block)
{{T:rep_hub}}

### Table S7c. Replication verdicts
{{T:rep_verdict}}

**BANC.** A BANC sweep (146 result files on the Modal Volume `flysim-data:/results_rep_v4/BANC*`) completed but was not retrieved or analysed because further cloud work was stopped (job recorded as failed; cost unknown, ≤ USD 1.76 by the container cap). MANC was not simulated dynamically.

---

## S8. Recoverability index (cross-validated; FAFB)

Source: `research/followup/index/tables/model_results.json` and `univariate_*.csv`; 1,118 of 1,209 dynamics runs processed (91 not processed); condition-grouped 5-fold cross-validation (4–20 repeats depending on the model) with cluster-bootstrap 95% CIs. Indices were computed on default-sign graphs. Intact FAFB: R_eff (mean-field) 19.62; confident-transmitter R_eff 7.38; signed 17.05 (`index/data/fafb_intact_index.json`). `ecore100` = fraction of neurons in the excitatory bootstrap core (maximal set in which every member receives enough excitatory synapses from other members to sustain firing at an assumed 100 Hz); `recip_n` = surviving reciprocal excitatory pairs relative to intact; `reach_sim` = simulator Reach.

### Table S8a. Post-reset FII_norm, non-runaway runs (top single indices)
{{IDX_FII}}

### Table S8b. Runaway after naive reperfusion among viable networks (top single indices)
{{IDX_RUN}}

Combined models: post-reset FII random forest CV R² 0.986 [0.974, 0.992]; recovered-after-reset logistic CV AUC 0.999; runaway among viable networks random forest CV AUC 0.61 [0.38, 0.76] (chance). Leave-arm-out transfer of `ecore100` for post-reset FII: arm D R² 0.89, arm H 0.84, arm E (ageing) −0.50. The planned BANC/MCNS/MANC recoverability maps were started (FAFB random-removal grid to f = 0.6 only) and not completed; no standalone INDEX_RESULTS report was written.

---

## S9. Deviations from the preregistration and task plan (all tracks)

1. **Seeds and stimulus rates.** 3–5 seeds per condition (preregistered ≥10); one rate per probe (preregistered 3). Baseline used 10 seeds.
2. **Spontaneous-activity metric.** The Shiu model is silent at rest, so FSA_baseline = 0 and the preregistered ratio is undefined. H5 used an evoked-neuron count (Act).
3. **Brain-death definition as run.** Spontaneous fraction < 0.01 and no readout response > 8 Hz (preregistered: < 2 SD of noise).
4. **Calibration choices made after inspecting runs.** J50 60 → 600 s, τ_ATP 20 → 30 s, q10_injury 2.5 added, w_syn 0.275 → 0.185. All are reported with the rejected alternatives.
5. **H1 dynamical test** was fitted on FII vs dead fraction from insult runs, not on removal curves.
6. **H2(a)** (death-time correlation with demand) was not computed.
7. **H3 protocols** were not matched for integrated supply.
8. **Temperatures** 20 and 26 °C were not run (17/22/29/33 °C were).
9. **Reach in dynamics** used the ≥5 matrix only, with no cut-off sensitivity. Ageing synapse loss is not reflected in Reach.
10. **LT50** was extrapolated, not simulated (≤1 h anoxia simulated).
11. **Structural removal orders** were below the preregistered ≥20 for several regimes (FAFB random 13, metabolic_area 10, metabolic_syn 5, degree 3 + 2 fine-grid, focal 5, matched 3; BANC 15/8/3/4/2; MCNS 6/6/2/3/2; MANC full) because of CPU contention. The f grid was reduced from 24 to 13 points mid-run, L∞ was capped at 30 hops for later orders, and GSCC was computed for orders 0–2 only. The primary metric samples 128 source neurons.
12. **Structural controls reduced.** One rewired null (FAFB, MANC only); ≥1 cut-off only for MANC (FAFB/MCNS files are floored at 5).
13. **A CSV column-shift bug** (rows without GSCC) affected 143 FAFB, 78 BANC and 65 MCNS rows. These were repaired deterministically.
14. **MAOL v1.1** was not available locally and was not analysed (planned deviation).
15. **Mechanism follow-up of 605 runs not executed** (Modal spend limit): 40/60-s blocks for all LHS sets, the full k_exc = 10 H3 grid, finer targeted-set brackets (0.2–0.3% at 30 s; 12.5–17.5% at 180 s), and a hysteresis ramp to −30 mV.
16. **Per-set FII baselines in arm L** used 4 seeds. Arm X k_exc = 0 cells reused dynamics arm C runs (identical code path).
17. **Dynamic replication.** Earlier MCNS sweeps (v1, adjacency bug; v2, infrastructure failure) were discarded. The v4 sweep gave 211 analysed runs (MCNS + FAFB-proxy). The BANC v4 sweep (146 files) was not harvested or analysed. MCNS and FAFB-proxy use proxies for demand and K⁺ coupling. MANC was not simulated dynamically.
18. **Recoverability index.** 91 of 1,209 runs were not processed. Cross-connectome maps were not completed. Indices use default signs. No standalone report was written; results are taken from `model_results.json`.
19. **Sign robustness.** Only the antennal-lobe local-neuron variant was run: FAFB full model 3 seeds at 30 s naive only (no variant-A intact baseline, so FII is relative to the default-sign baseline); MCNS and FAFB-proxy 5 seeds at 30/180 s naive and reset. Variants B and C were not run.
20. **Closed-loop release** was tested only in MCNS and the FAFB-proxy (180 s, default signs), not in the FAFB full model or across LHS sets.

---

## S10. Negative and null results (reported in full)

1. No sharp percolation threshold under random structural loss in any of four connectomes (H1 falsified; FAFB w 0.671).
2. No sharp threshold for FII vs dead fraction (ΔAIC +0.11 vs linear; width 1.56).
3. Graded, ramped or pulsed reperfusion was worse than abrupt reperfusion for identity and survival (H3 reversed). Adding an excitotoxic term did not change this for any comparison that included a reset.
4. Hyperpolarization (−10/−20 mV), GABA gain ×2–4, SD-release block, 3× clearance, noise, adaptation and cooling during reperfusion all failed to prevent the runaway.
5. Blocks of 2 and 5 s failed, and 10 s was unreliable at 360 s.
6. Preconditioning gave no protection (24.9% vs 24.8% dead).
7. Hyperthermia without hypoxia (33/37/40 °C) had no effect.
8. Metabolic support (reserve boost, clearance boost, mild hypometabolism) did not rescue the aged network.
9. "Resuscitated but different" did not exist outside the runaway state (0.0% coverage with the reset).
10. No single cell type is a structural keystone (≤1.03% reach loss).
11. Fly-like metabolic arrest alone did not produce hours-scale tolerance (extrapolated LT50 ≈ 1.51 h at J50 14,400 s vs fly 8.3 h).
12. No connectome index predicted naive-reperfusion runaway among viable networks (cross-validated best single-index AUC 0.53; combined 0.61 [0.38, 0.76]).
13. Targeted driver or hub blocks failed after 180 s of anoxia up to 10% of the brain.
14. A closed-loop release at K⁺ < 6 mM and ATP > 0.9 (reached after 11 s) re-ignited the runaway in MCNS (runaway 5/5) and the FAFB-proxy.
15. With antennal-lobe local neurons inhibitory, the 30-s high-rate runaway did not form (FAFB 0/3 runaway; MCNS and FAFB-proxy kept a weak 3.0 / 0.82 Hz/neuron oscillation classed runaway by the preregistered rule). The 30-s driver-core result is therefore not established.

---

## S11. Provenance and compute

- Dynamics: Modal jobs a4ddd6e8-9ca3-4803-a889-fa5fdd58ccf9, 31c05682-cf05-4a51-b769-97b562f60172, 4d88396b-d4af-4b98-935c-ea62f206acfb, 937a47ef-8836-46d7-8764-a2a4d3655763 (smoke 79482286-adfb-4c8d-b57c-64d67726799e); 1,209 runs, 0 errors, 24.3 CPU-h.
- Mechanism: smoke 900ff527-e669-4ffa-bcd0-5eed11e9297e; shard A 65780d67-c068-43a9-8cf0-1976a09783a9; shard B c6fe6ce5-53f5-4e4e-afe0-dc930dc25632; 2,061 runs, 0 errors. Follow-up 12a56f59-2ec3-4e3c-8934-253f564c9204 was killed by the spend limit and its results discarded.
- Replication: calibration 6fd134b9-a8a3-4e02-8438-143947ce30ea; v1 MCNS 6e86553b-044e-4155-becd-21845a5eb524 (invalid); v2 04f516ec-6e0a-4889-86d3-d5d9a903a4f6 (infrastructure failure).
- Replication (v4, analysed): smoke 14de3ae1-d0bf-412c-8841-c1a279754447 (MCNS 360-s hub arm, 45 runs); job A a9ef9894-0c76-4de2-89c9-6879f59167ff (166 runs). BANC job ef7d8554-1bbc-4a61-a027-11bb646177dc recorded as failed; 146 files on the Volume, not analysed.
- Sign variant (FAFB full model): local runs, `research/followup/mechanism/raw/sign_res/` (S_A_D30_naive_s0–2, SN_A_anat_D30_s0–2).
- Recoverability index: local; `research/followup/index/`.
- Structure: local 16-core host; pure NumPy bitset engine.
- Identity check: the mechanism-track simulator with k_exc = 0 reproduced dynamics run `C_D30_naive_s0` exactly (FII 0.3617, late rate 122.998 Hz/neuron).
- Run totals in the paper: 1,209 dynamics + 2,061 mechanism (784 of them in the LHS arm) + 211 replication + 6 FAFB sign-variant = 3,487 analysed simulated condition–seed runs.
