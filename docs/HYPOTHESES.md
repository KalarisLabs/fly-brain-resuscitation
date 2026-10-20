# Hypotheses (preregistration) and verdicts

Part 1 reproduces the preregistration from `research/lit/preregistration.md` unchanged; it was frozen before any simulation was run. Part 2 collects the verdicts from the track reports. Every verdict is internal to the model.

---

# Part 1 — Preregistration (verbatim)

## Preregistration — In-silico anoxia, "brain death" and resuscitation on Drosophila connectomes

**Registered:** 2026-10-16, before any insult/resuscitation simulation has been run or inspected by the literature track.
**Status of every statement below:** candidate hypothesis. None of this is a finding.
**Method guidance followed:** Kalaris `hypothesis-generation` skill v2.2 (observation → question → evidence boundary → rivals → discriminating predictions → operationalisation → falsification/controls → analysis plan).
**Evidence boundary:** sources and parameter values in `bibliography.md` (82 Crossref-verified records) and `parameters.csv` (81 rows; 6 flagged UNVERIFIED), search date 2026-10-16.

---

## 0. Scope and safety gate

| Item | Record |
|---|---|
| Owner / intended use | Project owner; exploratory computational research using connectomes as a proxy brain. |
| System | Leaky integrate-and-fire (LIF) network built from connectome graphs (Shiu et al. 2024 parameterisation) plus an added metabolic layer (ATP pool, Na+/K+-pump-dependent excitability, shared extracellular K+, damage accumulation, death). |
| Data | Public connectomes stored in `$RAW_ROOT`: FAFB v783 (139,255 neurons), BANC v888, MANC v1.2.1, MCNS v1.0. MAOL v1.1 not found locally, so it is excluded from H6 unless obtained (to be logged as a deviation). |
| Ethics / biosafety | Purely computational; no animals, humans, pathogens or patient data. No clinical claims. |
| Claim type | **Model-mechanistic / descriptive of the model.** All results describe the behaviour of a simplified model; they are not evidence about clinical brain death or the reversibility of human death. |
| Key limitation, stated in advance | The fly CNS enters SD as a *protective, reversible* shutdown (adult LT50 for anoxia ≈ 8 h); mammalian anoxic depolarisation begins injury within 2–5 min. Simulated "irreversibility" depends on a damage rule that we write ourselves, so any threshold reflects model assumptions plus connectome structure. |
| "Brain death" in the model | Operational state only: (i) <1% of neurons spontaneously active for ≥1 s of simulated time under baseline drive, **and** (ii) no readout responds to any stimulus (response < 2 SD of noise), **and** (iii) this state persists after the reperfusion protocol ends. This is an analogue of the *construct* in the AAN 2023 guideline (10.1212/WNL.0000000000207740), not a clinical determination. |

---

## 1. Frozen observations motivating the hypotheses (from the literature)

1. Insect anoxia triggers SD within seconds (D. melanogaster: 14.2 ± 0.3 s to anoxic SD at room temperature; 10.1242/jeb.179598). [K+]o rises from ~6–13 mM to ~44–58 mM (10.1242/jeb.074468; 10.1371/journal.pone.0028994; 10.1242/jeb.179598).
2. Repeated anoxia, or Na+/K+-ATPase block, can end in a terminal [K+]o plateau (54.1 ± 13.9 mM) with no recovery, followed by death (10.1371/journal.pone.0028994). So the fly CNS has an irreversible end state.
3. Recovery slows nonlinearly with anoxia duration: first movement within ≤15 min after <2 h of anoxia, ~60 min after 4 h (10.1016/S0022-1910(96)00084-4 as cited in 10.1242/jeb.177147).
4. The reperfusion pattern matters. Long anoxia with short normoxia at high cycle frequency degrades fly K+ homeostasis (10.1242/jeb.074468). Fly mortality peaks at reoxygenation (10.1371/journal.pone.0005422). In rats, gradual reperfusion and postconditioning reduce infarct size (10.1002/jnr.21703).
5. Cooling protects. Fly [K+]o drift during repetitive anoxia is ~0 at 17 °C vs 333 mM cumulative at 29 °C (10.1242/jeb.074468). SD propagation has Q10 ≈ 2.1–3.0 (10.1152/ajpregu.00040.2025). Clinical hypothermia improves outcome (10.1056/NEJMoa012689).
6. Recovery at the cellular level can dissociate from network function. BrainEx restored cellular and synaptic activity 4 h post-mortem without global ECoG activity (10.1038/s41586-019-1099-1).
7. Signalling energy is dominated by Na+ influx that the pump must reverse (10.1097/00004647-200110000-00001; 10.1242/jeb.017574). Per-AP cost can differ >10-fold between neuron types (10.1371/journal.pcbi.1000840). So "large neurons cost more" is a proxy assumption, not a law.
8. Networks with heavy-tailed degree tolerate random loss but collapse under targeted hub removal (10.1038/35019019; 10.1103/PhysRevLett.85.5468).

---

## 2. Shared operational definitions (fixed before simulation)

- **Stimulus set S.** Per connectome, a prespecified list of sensory populations driven by Poisson input. For FAFB, BANC and MCNS this includes at least sugar GRNs, water GRNs, bitter GRNs, Johnston's-organ/antennal mechanosensory neurons, and LC4/LC9 visual projection neurons. For MANC it is leg/wing sensory neuron classes. Each population is driven at 3 rates (e.g. 50, 100, 200 Hz). The list is frozen in code before the first insult run, and any change is logged as a deviation.
- **Readout set M.** Motor neurons, descending neurons (e.g. DNa02, giant fibre DNp01, MN9 where present) or VNC motor neurons for MANC. Response = mean firing rate during stimulation minus pre-stimulus rate.
- **Response matrix.** R ∈ ℝ^{|S|×|M|}, averaged over n ≥ 10 Poisson seeds.
- **Functional Identity Index (FII).** Pearson r between vec(R_condition) and vec(R_baseline), computed over entries where the baseline response > 0 or the condition response > 0. If the condition matrix is all zero, FII is undefined and set to 0. **Noise ceiling** FII_ceil = median split-half r between two disjoint baseline seed sets. Report FII and FII_norm = FII / FII_ceil. A secondary rank-based variant (Spearman) is reported for robustness.
- **Sensory-to-motor reachability (Reach).** Fraction of (s, m) pairs with at least one directed path from s to m through surviving neurons, using edges with ≥5 synapses (the FlyWire convention) and path length ≤ 6. A weighted variant uses the summed product of normalised weights.
- **Fraction spontaneously active (FSA).** Fraction of surviving + total neurons (both reported) firing ≥ 1 spike/s in a 1 s window under weak background drive, without stimulus.
- **Damage fraction f.** Fraction of neurons dead (structural removal) at the end of the protocol.
- **Recovered.** FII_norm ≥ 0.8 **and** Reach ≥ 0.8 × baseline, at the end of the post-reperfusion period.
- **Uncertainty.** 95% CIs by bootstrap over seeds (2,000 resamples) and, for removal curves, over ≥ 20 independent removal orders.
- **Multiplicity.** Holm correction across the six primary tests (one per hypothesis). Secondary analyses are labelled exploratory.

---

## 3. Hypotheses, rivals, predictions and falsification thresholds

### H1 — A sharp, percolation-like critical damage threshold exists

- **Candidate (H1).** As damage fraction f increases, FII and Reach do not decline linearly. They stay near baseline and then collapse over a narrow interval around f_c. Beyond f_c, no reperfusion protocol in the tested grid restores function ("recovered" criterion unmet).
- **Rivals.** H1-R1: graceful, approximately linear decline (distributed redundancy). H1-R2: exponential decay with no inflection (independent per-pathway failure). H1-R3: an apparent threshold produced by the readout definition (e.g. the 5-synapse edge cut-off or the firing threshold of a few readout neurons), not by network structure.
- **Discriminating predictions.** Under H1, a 4-parameter logistic fits FII(f) and Reach(f) better than linear and exponential models (ΔAIC ≥ 10). The 10–90% transition width w ≤ 0.10 in absolute f. Across-order variance of Reach peaks within ±0.05 of f_c, the finite-size susceptibility signature. Under R1/R2, w > 0.25 or the logistic is not preferred.
- **Controls.** (a) Degree-preserving rewired null graph (same in/out-degree sequence), which tests whether the threshold is structural. (b) Edge-threshold sensitivity: repeat Reach with ≥1 and ≥10 synapse cut-offs; a threshold that appears only at one cut-off supports R3. (c) Negative control: removing neurons with no path to any readout should leave FII unchanged (|ΔFII| < 0.05 up to the same f).
- **Support threshold.** ΔAIC(logistic − best alternative) ≤ −10 **and** w ≤ 0.10 **and** irrecoverability for f > f_c + 0.05 under all tested reperfusion protocols, in the FAFB primary analysis.
- **Falsified if.** w > 0.25, **or** ΔAIC > −2, **or** any protocol yields "recovered" at f ≥ f_c + 0.10.
- **Indeterminate.** 0.10 < w ≤ 0.25, or −10 < ΔAIC ≤ −2: report as "graded transition"; no claim of criticality.

### H2 — Vulnerability is non-random, and targeted loss is more destructive than random loss

- **Candidate (H2).** (a) Under a metabolic insult where ATP demand scales with neuron size (proxy: total synapse count or cable length where available), high-demand neurons and hubs die earlier. (b) Removing neurons in order of demand or degree destroys function at a lower damage fraction than random removal.
- **Rivals.** H2-R1: death order is determined by network position (input drive), not intrinsic demand. H2-R2: targeted ≈ random, because the fly connectome is not hub-dominated enough. H2-R3: the demand–death correlation is circular, since demand is defined from synapse count and synapse count also sets drive.
- **Discriminating predictions.** (a) Spearman ρ between demand proxy and time-of-death ≤ −0.3 (higher demand, earlier death). (b) Ratio f_c,targeted / f_c,random ≤ 0.7, with the bootstrap 95% CI excluding 1. Under R1, ρ with demand becomes ~0 once input drive is partialled out. Under R2, the ratio CI includes 1.
- **Controls.** (a) Demand decoupled from degree: rerun with uniform per-neuron demand; if the death order is unchanged, R1 is supported. (b) Synapse-count-matched random removal (random neurons drawn from the same degree bins). (c) Targeted removal on the rewired null, to separate degree-distribution effects from wiring specifics.
- **Support threshold.** Both (a) and (b) met in FAFB, and (a) survives partial correlation on input drive (|partial ρ| ≥ 0.2).
- **Falsified if.** ρ > −0.1 **and** the ratio CI includes 1. Mixed outcomes, where one criterion is met and the other is not, are reported as partial support.

### H3 — Staged/gradual reperfusion beats abrupt reperfusion

- **Candidate (H3).** Abrupt restoration of O2/ATP production after an insult near f_c triggers a synchronous re-excitation surge. In the model, that surge adds activity-dependent (excitotoxic) damage. A ramped or staged (stuttering) reperfusion lowers the peak surge and preserves more neurons and higher FII.
- **Rivals.** H3-R1: no difference, because damage is set during the insult. H3-R2: abrupt is better, because faster ATP recovery shortens time spent depolarised. H3-R3: any staged advantage comes from the confound of different total O2 delivery or total time.
- **Protocols (fixed).** Abrupt step to 100% supply. Linear ramp over T_r ∈ {0.5, 1, 2, 5} model time units. Postconditioning-like on/off cycles (inspired by 3×30 s/10 s and 10×10 s/10 s; 10.1002/jnr.21703), rescaled to model time. All protocols are matched for total integrated supply and total observation time (control for R3).
- **Discriminating predictions.** Under H3, the peak population firing rate during reperfusion is lower for staged protocols. FII_end(staged) − FII_end(abrupt) ≥ 0.05 and survival difference ≥ 2 percentage points, with CIs excluding 0, for insults in the window f_c ± 0.1 (or equivalent duration). R2 predicts the opposite sign.
- **Negative control.** With the excitotoxic damage term set to zero, staged and abrupt should not differ (|ΔFII| < 0.02). If they still differ, the advantage does not arise from the hypothesised mechanism.
- **Support threshold.** ΔFII ≥ 0.05 (CI excludes 0) in ≥ 2 of 3 insult severities near threshold, **and** the negative control is null.
- **Falsified if.** ΔFII ≤ 0, or the CI includes 0 at all severities, or the benefit persists unchanged with the excitotoxic term removed.

### H4 — Hypothermia widens the recovery window approximately by Q10

- **Candidate (H4).** Lowering temperature by ΔT scales metabolic demand, pump rate and damage kinetics by Q10^(ΔT/10). This lengthens the maximum insult duration that still meets the "recovered" criterion (D_max) by roughly Q10^(−ΔT/10). The prior Q10 range is 2–3 (SD propagation Q10 2.13–2.98 in flies, 10.1152/ajpregu.00040.2025; general metabolic Q10 ≈ 2.2–2.4, derived and flagged as unverified).
- **Circularity warning.** Q10 is an input to the model, so recovering exactly Q10 is partly by construction. The informative test is whether D_max scales **by** Q10 (proportional regime) or deviates from it because of threshold nonlinearities in the network.
- **Rivals.** H4-R1: no widening (window set by temperature-insensitive passive leak or structural damage). H4-R2: super-Q10 widening (effective Q10_window > 4), with the threshold amplifying small kinetic changes. H4-R3: sub-Q10 widening (1 < Q10_window < 1.5), because cooling also slows recovery processes (pump-mediated K+ clearance).
- **Design.** Temperatures 17, 20, 23, 26, 29 °C (bracketing the Rodríguez & Robertson 2012 values), with model Q10 ∈ {2.0, 2.5, 3.0}. Separate Q10 values for demand and for recovery are explored as a sensitivity analysis.
- **Discriminating predictions / thresholds.** Estimate Q10_window = (D_max(T1)/D_max(T2))^(10/(T2−T1)). Q10_window within [0.67, 1.5] × model Q10 → "approximately Q10" (H4 supported). Q10_window > 1.5 × model Q10 → R2. 1 < Q10_window < 0.67 × model Q10 → R3.
- **Falsified if.** The D_max ratio CI includes 1 across 17 vs 29 °C (no widening), **or** Q10_window falls outside [0.67, 1.5] × model Q10 (here H4's "≈ Q10" clause is rejected while widening is still reported).

### H5 — Functional identity decays before global activity and structure ("resuscitated but different")

- **Candidate (H5).** Across insult severities and after reperfusion, FII declines at lower damage than FSA or Reach. So states exist where the network is active and structurally connected but maps stimuli to responses differently from baseline.
- **Rivals.** H5-R1: all three metrics decline together, so identity loss is just loss of activity or structure. H5-R2: an apparent early FII decline reflects measurement noise (FII has lower reliability), not identity change. H5-R3: FII drops because responses are globally scaled, not re-mapped. Pearson r is scale-invariant, so this predicts FII stays high while magnitude drops. Magnitude-ratio is reported to check.
- **Metrics.** For each metric X ∈ {FII_norm, FSA/FSA_0, Reach/Reach_0}, D50_X = damage (or insult duration) at which X falls halfway from baseline to its floor.
- **Discriminating predictions.** D50_FII < min(D50_FSA, D50_Reach) by ≥ 0.05 in f (bootstrap CI of the difference excludes 0). A "resuscitated-but-different" region exists: FSA ≥ 0.8 × baseline **and** Reach ≥ 0.8 × baseline **and** FII_norm ≤ 0.5, covering ≥ 5% of the explored (severity × protocol × temperature) grid.
- **Controls.** (a) Noise ceiling: FII_norm uses split-half reliability, addressing R2. (b) Label-shuffled negative control: FII between baseline and baseline with readout labels permuted should be ≈ 0, setting the floor. (c) Rate-scaled baseline (all responses ×0.5) should give FII ≈ 1, confirming FII is insensitive to pure gain change (R3).
- **Support threshold.** Both the D50 ordering and the ≥5% regime criterion are met in FAFB.
- **Falsified if.** D50_FII ≥ D50_FSA or D50_FII ≥ D50_Reach (CI-based), **and** the regime covers < 1% of the grid. Coverage of 1–5% is reported as a "marginal regime".

### H6 — Findings replicate across FAFB, BANC, MANC and MCNS

- **Candidate (H6).** The qualitative outcomes of H1–H5 (direction of effect; threshold present or absent) are the same across connectomes, and normalised f_c values are similar.
- **Rivals.** H6-R1: dataset-specific results driven by proofreading completeness, neurotransmitter-prediction differences, sex, or brain-only vs CNS vs VNC-only scope. H6-R2: replication only for structural metrics (Reach), not dynamical ones (FII).
- **Comparability rules (fixed).** Identical simulator code and parameters. Stimulus and readout sets are matched by cell type where types exist in multiple datasets (cross-dataset typing, 10.1038/s41586-024-07686-5). MANC is analysed with VNC-appropriate S and M and compared only on direction of effect. Neurotransmitter sign is assigned by the same rule in every dataset.
- **Support threshold.** For each of H1–H5, the FAFB outcome category (supported / falsified / indeterminate) matches in ≥ 3 of 4 connectomes. For H1, f_c values lie within a factor of 1.5 of each other (max/min) among datasets where H1 is supported.
- **Falsified if.** ≤ 2 of 4 connectomes match for the majority of H1–H5, or f_c max/min > 2.
- **Planned deviation.** MAOL (optic lobe only) is not available locally and is excluded. If it is later obtained, it will be analysed as exploratory for H1/H2 only.

---

## 4. Analysis plan and decision rules

1. **Baseline validation first.** The FAFB LIF must reproduce Shiu et al. 2024 qualitative results (sugar-GRN drive activates MN9; 0 Hz baseline) before any insult run. If it fails, no hypothesis is tested; the failure is reported.
2. **Primary dataset.** FAFB v783. Replication datasets follow H6.
3. **Seeds and replication.** ≥ 10 Poisson seeds per condition; ≥ 20 removal orders for percolation curves. Report effect sizes with 95% bootstrap CIs; p-values are secondary and Holm-corrected across the six primary tests.
4. **Parameter provenance.** Every metabolic parameter used is drawn from `parameters.csv` and cited by DOI. Parameters flagged UNVERIFIED are sensitivity-tested over ≥ ±50% or a literature range before any conclusion relies on them.
5. **Sensitivity analyses (prespecified).** Wsyn ∈ {0.2, 0.275, 0.35} mV; damage-threshold rule ±50%; excitotoxic gain ×{0, 0.5, 1, 2}; edge cut-off {1, 5, 10} synapses. A hypothesis whose outcome category flips within these ranges is reported as "not robust".
6. **No outcome-dependent changes.** If definitions or thresholds must change after seeing results, both the original and the changed analysis are reported, and the change is logged with date and reason in `deviations.md`.
7. **Interpretation boundary.** Support for any hypothesis means "the model, under these assumptions, behaves this way". No result will be phrased as showing that biological brain death is or is not reversible.

## 5. Rival-prediction matrix (summary)

| H | Key observable | H predicts | Main rival predicts | Falsifying result |
|---|---|---|---|---|
| H1 | FII(f), Reach(f) | logistic, width ≤ 0.10, irrecoverable beyond f_c | linear/exponential, width > 0.25 | width > 0.25 or ΔAIC > −2 or recovery at f ≥ f_c+0.10 |
| H2 | ρ(demand, death time); f_c,targ/f_c,rand | ρ ≤ −0.3; ratio ≤ 0.7 | ρ ≈ 0; ratio ≈ 1 | ρ > −0.1 and ratio CI ∋ 1 |
| H3 | ΔFII staged−abrupt; peak surge | ≥ +0.05; lower surge | ≤ 0 | ΔFII ≤ 0, or benefit survives with excitotoxic term = 0 |
| H4 | Q10_window | within [0.67, 1.5] × model Q10 | ≈ 1 (no widening) or > 1.5× | D_max ratio CI ∋ 1, or outside band |
| H5 | D50 ordering; regime coverage | D50_FII lowest; ≥ 5% coverage | equal D50s; < 1% | D50_FII not lowest and coverage < 1% |
| H6 | outcome concordance | ≥ 3/4 datasets | ≤ 2/4 | ≤ 2/4 for most H, or f_c max/min > 2 |

## 6. Known limitations recorded in advance

- Connectome LIF models lack gap junctions, neuromodulation, glia, extracellular-space geometry and dendritic compartments; the metabolic layer is phenomenological.
- Neuron "size/demand" is approximated from connectome-derived proxies (synapse count, optionally cable length), and per-spike energy cost varies >10-fold between neuron types (10.1371/journal.pcbi.1000840).
- Fly SD is largely protective; mammalian timing (anoxic depolarisation 2–5 min; terminal SD ~3.9 min in humans, 10.1002/ana.25147) is used only for qualitative scaling comparisons.
- Several mammalian values commonly quoted (pump ≈ 50% of brain ATP; [K+]o ≈ 60 mM at anoxic depolarisation; ~5 min irreversibility rule) could not be verified from retrieved text and are flagged UNVERIFIED.


---

# Part 2 — Verdicts (copied from track reports)

### 4. Verdicts on preregistered predictions (structural parts only)

| hypothesis | criterion | result | verdict |
|---|---|---|---|
| **H1** sharp threshold (FAFB primary, Reach_N L6) | ΔAIC ≤ −10 **and** w ≤ 0.10 (irrecoverability part belongs to the simulator) | random: ΔAIC -830.1, w 0.671; metabolic area: ΔAIC -436.2, w 0.421; hub-targeted (fine): ΔAIC -35.7, w 0.231 | **Falsified** for random and metabolic loss (w > 0.25); **indeterminate / graded** for hub-targeted loss. The FAFB random-loss verdict holds at ≥10 synapses (w 0.960). The hub-targeted category flips (indeterminate → falsified at ≥10) and is therefore not robust. |
| **H2(b)** targeted more destructive | f_c,targeted / f_c,random ≤ 0.7, CI excludes 1 (FAFB) | metabolic area 0.346 [0.337, 0.356]; metabolic synapse 0.285 [0.278, 0.291]; hub 0.041 [0.041, 0.041] | **(b) met** in FAFB. H2-R2 rejected. (a), the death-time correlation, is for the simulator, so H2 overall is at least partial support pending (a). Caveat H2-R3: membrane area correlates with degree (Spearman ρ = 0.675), so the "metabolic" ordering is partly a hub ordering. |
| **H2 control (b)** degree-matched random | — | FAFB 0.044 [0.044, 0.045] ≈ hub-targeted | Destructiveness is explained by degree, not by hub identity. |
| **H2 control (c)** targeted on rewired null | — | FAFB null hub f_c 0.266 vs real 0.035 | Real wiring amplifies hub vulnerability beyond the degree sequence. |
| **H5-structural** | D50 ordering of structural metrics | Random loss: weighted influence and GSCC < Reach L4 < Reach L6 < Reach ∞ in 4/4 datasets. Targeted loss: influence ≈ L4 (no consistent order) | Structural precondition for "resuscitated but different" **present under random loss, not under targeted loss**. Dynamical H5 is untested here. |
| **H6** replication | outcome category matches FAFB in ≥3/4; f_c max/min ≤ 1.5 (supported datasets) / falsified if > 2 | H1 random: 4/4 falsified (f_c max/min 1.37, but MCNS and MANC midpoints extrapolated beyond 0.9). H1 hub-targeted: 1/4 match (BANC supported; MCNS and MANC falsified). H2(b): ratio ≤ 0.7 with CI excluding 1 in 4/4 for metabolic_syn and hub-targeted. f_c for metabolic_syn max/min 6.32; hub-targeted 23.35 (MANC coarse-grid f_c 0.388 in range; 25.36 if the extrapolated fine-grid 0.421 is used) | **Direction of effect replicates** (H1 absent under random loss; targeted ≫ random). **Magnitudes do not**: f_c max/min > 2 for targeted regimes, falsifying the quantitative clause of H6. The VNC-only MANC is the outlier (much more hub-tolerant). |

Holm correction across the six primary tests is not applied here, because the structural track produces bootstrap CIs, not p-values, for H1/H2. The metabolic and hub-targeted ratio CIs exclude 1 by a wide margin (largest upper bound 0.672, MANC metabolic_syn). Focal-lesion CIs reach 0.91–1.01 and are not treated as support.


### 8. Hypothesis tests (preregistration H1, H3, H4, H5; H2/H6 belong to the Structure track)
*[Figure: H4 and H1 — see [docs/figures/dynamics/](../figures/dynamics/)]*

*Fig. 4 — (a) D_max vs temperature with Q10 scaling predictions. (b) FII_norm vs dead fraction for non-runaway reset runs (n = 372).*

* **H1 (sharp percolation-like threshold): FALSIFIED.** Pooled B/C/H reset runs without runaway (n = 372): logistic ΔAIC vs the best alternative (linear) = +0.11 (> −2). 10–90 % width = 1.56 (> 0.25). FII declines approximately linearly with dead fraction. Recovered runs reached up to 0.325 dead (focal optic lobe), while global anoxia lost recovery at ~10 % dead. Under this battery, *where* neurons die matters more than how many.
* **H3 (staged/gradual reperfusion beats abrupt): FALSIFIED (direction reversed, rival H3-R2).** At 360 s, every staged protocol had lower FII_norm than naive: ramp30 Δ = -0.331 [-0.583, -0.130]; ramp60 -0.297 [-0.540, -0.092]; pulsed3x10 -0.221 [-0.460, -0.011]. Staged protocols also had *more* deaths (ramp60 +6.6 pp, ramp30 +4.0 pp; CIs exclude 0), because slower re-energization prolongs depolarization. The peak surge was not lower (ramp30 283 vs naive 265 Hz at 360 s). The only positive cell was ramp10 at 30 s (+0.46, both arms runaway). The negative control was not null (ramp30_noinj − naive_noinj: -0.23 / +0.14 / -0.28), so FII differences between runaway states are attractor noise, not excitotoxicity. The model has no spike-driven injury term, so H3's mechanism is absent by construction. The decisive variable is the block reset (Δ FII_norm vs naive = +1.06 / +0.87 / +0.49, all CIs exclude 0).
* **H4 (hypothermia widens the window by ≈Q10): SUPPORTED.** D_max 17 °C / 29 °C ratio CI = [2.29, 2.57] (excludes 1). Q10_window = 2.09 [1.99, 2.19] at model Q10 2.5; acceptance band [1.67, 3.75]. Sensitivity: model Q10 2.0 → 1.71 [1.60, 1.83] (band [1.33, 3.0]); Q10 3.0 → 2.44 [2.26, 2.64] (band [2.0, 4.5]). All three fall inside the 'approximately Q10' band, but sit consistently at Q10^0.8 (2.5^0.8 = 2.08, 2^0.8 = 1.74, 3^0.8 = 2.41). That is exactly what the Hill-4 hazard predicts (cumulative hazard ∝ rate⁴·D⁵ ⇒ D_max ∝ rate^−4/5). Per the prereg circularity warning, the widening is set by the damage rule, not by network nonlinearity.
* **H5 (identity lost before activity/structure): NOT SUPPORTED as preregistered; mixed.** With the block reset, D50 values coincide: FII_norm 374 s [363, 385], Act 368 s [365, 373], Reach 390 s [388, 391]. So D50_FII ≥ D50_Act, and the 'resuscitated-but-different' region excluding runaway covers 0.0 % of runs ⇒ falsified for reset recovery. With naive reperfusion, D50_FII = 42 s [21, 59] ≪ Reach 388 s, but the CI overlaps Act (291 s [21, 384]). The active + connected + FII_norm ≤0.5 region covers 22.1 % of runs (22.8 % of conditions), all of them runaway. 'Resuscitated but different' exists in this model only as the pathological runaway attractor; once it is cleared, identity tracks structure.
* Multiplicity: H1, H3 and H5 fail on effect-size/structure criteria, not on p-values, so Holm correction does not change any verdict. H4's CI excludes 1 by a wide margin.


### 6. Hypothesis and claim verdicts
* **Runaway attractor is robust to parameter uncertainty: SUPPORTED.** It occurred in 48/48 LHS sets and 288/288 runs; it was high-rate in 41/48.
* **20-s block reset generalizes: SUPPORTED at 30 s (48/48), CONDITIONAL at 180 s (34/48).** Failures are re-ignition (slow ATP recovery, low reserve) or structural loss (low J50).
* **H3 (graded beats abrupt reperfusion) with an excitotoxic term: still FALSIFIED for FII and for every block combination.** There is a weak, non-decisive H3-direction effect on death without a block for pulsed / 60-s ramp at k_exc = 3 (−0.8 to −1.6 pp). The H3 reversal from the dynamics track is not an artefact of the missing excitotoxicity term, within k_exc ≤ 3. k_exc = 10 is provisional (abrupt vs block only).
* **Attractor has an identifiable driver core: SUPPORTED.** The top 0.5 % sustain/ignite set is a single strongly connected excitatory core (81–96 %, vs 40 % hubs and 0.4 % random), centred on antennal-lobe local and projection neurons.
* **Targeted reset can replace whole-brain block: SUPPORTED after 30 s** (0.5 % driver block rescues 5/5 held-out seeds; hubs need 2 %, random fails at ≤20 %). **Largely NOT SUPPORTED after 180 s** (20 % needed; hubs and random fail at ≤20 %).
* **Bistability and hysteresis in the intact network: SUPPORTED.** Quiet and active states coexist at the same parameters. The basin boundary follows the analytic K⁺ threshold at δ ≥ 0. Hysteresis is ≥13 mV wide. The escape time is ≥10–20 s. The mechanism is ionic (K⁺) and survives removal of all synapses.



> TODO(PENDING-H6-dynamic): dynamic replication of H1/H3/H4/H5 on MCNS/BANC — see [PENDING.md](PENDING.md).
