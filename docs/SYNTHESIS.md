# SYNTHESIS — final claim-by-claim verdicts (v1.0, 2026-10-20)

Scope: every statement below concerns a connectome-constrained leaky integrate-and-fire (LIF) model (Shiu et al. 2024) with a phenomenological ATP / Na⁺-K⁺-pump / extracellular-K⁺ / injury layer. The metabolic constants are not fitted to *Drosophila* recordings. Nothing here is evidence about living tissue, clinical brain death or human resuscitation.

Sources (all under `D:\FRUITFLY\research\`):
- `lit/preregistration.md` (H1–H6, registered 2026-10-16)
- `structure/STRUCTURE_RESULTS.md` (four connectomes, percolation)
- `dynamics/DYNAMICS_RESULTS.md` (FAFB, 1,209 runs)
- `followup/mechanism/MECHANISM_RESULTS.md` (FAFB, 2,061 runs) and `followup/mechanism/raw/sign_res/` (6 FAFB sign-variant runs, analysed here)
- `followup/replication/REPLICATION_RESULTS.md` + `tables/v4_*.{csv,json}` (MCNS + FAFB-proxy, 211 runs)
- `followup/index/tables/model_results.json` and `INDEX_INTERIM.md` (recoverability index)
- `ara/logic/claims.md` (C01–C14)

Analysed simulation runs: 1,209 + 2,061 + 211 + 6 = **3,487** condition–seed runs.

Confidence scale: **high** = consistent across seeds, parameter sets and (where tested) connectomes, with CIs well clear of the criterion; **moderate** = consistent but tested in one dataset, few seeds, or with a known confound; **low** = n ≤ 3, untested replication, or a known artefact risk.

---

## 0. Revised headline

**Abrupt reoxygenation can trap the model brain in a self-sustained, K⁺-coupled high-activity state, and a 10–20 s total activity block releases it. Whether the trap forms depends on how long the anoxia lasted and, for short anoxia, on transmitter sign.**

- **Long anoxia (180 s): robust.** The runaway occurs in 48/48 Latin-hypercube parameter sets (FAFB) and in MCNS 5/5. It persists when antennal-lobe (AL) local neurons are set inhibitory (MCNS FII_norm −0.135, FSA 0.897; FAFB-proxy −0.125). The total block still rescues (MCNS 0.989; FAFB-proxy 0.965).
- **Short anoxia (30 s): sign-sensitive.**
  - With the default sign rule (unknown/monoamine → excitatory), it occurs in FAFB 5/5, MCNS 5/5 and LHS 48/48.
  - With the 147 excitatory-assigned AL local neurons (of 431 in the AL-LN class/type mask) set inhibitory, the high-rate runaway disappears:
    - FAFB full model, variant A: 3/3 seeds no runaway (late rate 0.035–0.17 Hz/neuron, FSA ≤ 0.005, FII 0.915–0.936);
    - MCNS: FII_norm 1.000 [0.996, 1.002];
    - FAFB-proxy: FII_norm 0.859.
- **Release timing matters.** A closed-loop release, gated on K⁺ < 6 mM and mean ATP e > 0.9, was triggered after 11 s; the network re-ignited (MCNS runaway 5/5, 11.3 Hz/neuron; FAFB-proxy 40.2 Hz/neuron).

**Flag (file contradicts the task summary):** "the short-anoxia runaway disappears" is exact for the FAFB full model (variant A, 3/3 below both runaway thresholds). It is not exact for MCNS or the FAFB-proxy.
- In both, a low-rate oscillation remains: MCNS late rate 3.0 Hz/neuron, FSA 0.022; FAFB-proxy 0.82 Hz/neuron.
- This still meets the preregistered runaway rule (>0.5 Hz/neuron), so both are classed "runaway" 5/5 in `v4_condition_summary.csv`.
- Identity is nevertheless preserved (FII_norm 1.000 and 0.859).
- Correct wording: *the high-rate runaway and the identity loss disappear at 30 s.* In MCNS and the FAFB-proxy a weak low-rate oscillation persists that the preregistered classifier still counts as runaway.

Mechanistic reading (model-internal interpretation, not separately tested):
- The dissection arm shows the *maintained* high-rate state is ionic: it survives removal of all synapses and collapses with a K⁺ clamp.
- *Entry* into it after short anoxia needs synaptic recruitment (removing excitatory synapses cut the 30-s runaway to 14 Hz).
- The AL local-neuron loop supplies that recruitment under the default signs.
- After long anoxia the K⁺ load alone (K⁺ at the 80-mM ceiling) places the network in the active basin.

---

## 1. Preregistered hypotheses

| H | Verdict | Key evidence (file → numbers) | Confidence |
|---|---|---|---|
| **H1** sharp percolation-like threshold | **Falsified** (random structural loss and dynamical damage); hub-targeted structural loss graded/indeterminate | STRUCTURE §3.2/§4: FAFB random w = 0.671, ΔAIC −830.1 (w > 0.25 bound) in 4/4 datasets; still falsified at ≥10 synapses (w 0.960). Hub-targeted FAFB w 0.231 (indeterminate); only BANC meets sharpness (w 0.086, ΔAIC −84.7); category flips with edge cut-off. DYNAMICS §8: FII vs dead fraction, ΔAIC +0.11 vs linear, width 1.56 (n = 372) | High (random/dynamical); low (hub-targeted category not robust) |
| **H2** targeted vulnerability | **Partial**: (b) supported, (a) not tested | STRUCTURE §4: FAFB f_c ratio area 0.346 [0.337, 0.356], synapse 0.285 [0.278, 0.291], hub 0.041 [0.041, 0.041]; ratio ≤ 0.7 with CI excluding 1 in 4/4 connectomes. Degree-matched random reproduces hub curve (0.044). H2(a) death-time vs demand correlation never computed | High for (b); (a) untested |
| **H3** graded > abrupt reperfusion | **Falsified (direction reversed)** | DYNAMICS §8: 360 s ramp30 Δ −0.331 [−0.583, −0.130], ramp60 −0.297, pulsed3x10 −0.221; ramp60 +6.6 pp dead. MECHANISM §3: with excitotoxic term no graded protocol beat abrupt + block at k_exc ≤ 3 (Δ FII_norm −0.045 to +0.006); weak H3-direction death effect without block (pulsed −1.6 pp [−1.8, −1.4] at 180 s, k_exc 3) while still in runaway. k_exc 10 abrupt-vs-block only | High within model (k_exc ≤ 3); provisional at k_exc 10; tested only under default signs |
| **H4** hypothermia widens window ≈ Q10 | **Supported**, but set by the injury rule | DYNAMICS §8: D_max 359/271/149/124 s at 17/22/29/33 °C; Q10_window 2.09 [1.99, 2.19] at model Q10 2.5 (≈ Q10^0.8, predicted by Hill-4 hazard). REPLICATION §2d: MCNS 17 − 29 °C after 180 s + block +0.185 [+0.123, +0.260] (direction replicates) | High that it holds; low informativeness (circular by construction) |
| **H5** identity lost before activity/structure | **Not supported as preregistered**; dissociation exists only inside the runaway state | DYNAMICS §8: with reset D50 FII 374 s [363, 385] ≈ Act 368 s ≈ Reach 390 s; regime coverage 0.0 %. Naive: 22.1 % of runs in active + connected + FII ≤ 0.5, all runaway | High (FAFB) |
| **H6** replication across connectomes | **Partial** | Structure: direction replicates 4/4 (H1 absent under random loss; targeted ≫ random); magnitudes do not (f_c max/min 6.32 synapse-weighted, 23.35 hub). Dynamics: MCNS replicates runaway, block rescue, hub-protection efficiency (84× vs FAFB ~30×) and temperature direction (REPLICATION §3). BANC dynamics present on the Modal Volume (146 files) but not analysed; MANC not simulated dynamically. The preregistered "≥3 of 4 connectomes" rule cannot be evaluated for dynamics (2 of 4 tested: FAFB, MCNS) | Moderate (structure); moderate (dynamics, single replication dataset; sex and dataset confounded) |

---

## 2. ARA claims C01–C14

### C01 — Naive reperfusion ignites a self-sustained runaway attractor
**Verdict: partial (revised).** Supported after ≥180 s anoxia. Sign-sensitive after 30 s anoxia. The "≥7 s" threshold is established only under the default sign rule.

Evidence under the default signs:
- DYNAMICS §3: 30/30 FAFB runs from 7 to 480 s at 22 °C ran away.
- MECHANISM §2: 48/48 LHS sets and 288/288 runs, at both 30 and 180 s.
- REPLICATION §2a: MCNS 5/5 at 30, 180 and 360 s and at 17/22/29 °C.

Evidence with AL-LN signs corrected:
- 30 s (`sign_res/S/*.json`; REPLICATION sign section): FAFB variant A 0/3 runaway; MCNS FII_norm 1.000 with a 3.0 Hz/neuron residual oscillation; FAFB-proxy 0.859 with 0.82 Hz/neuron.
- 180 s: runaway persists (MCNS −0.135, FAFB-proxy −0.125).

Under corrected signs, the duration at which runaway first appears lies between 30 and 180 s and was not localised. The LHS did not vary AL-LN sign; it varied glutamate sign, which was the strongest determinant of runaway intensity (PRCC −0.71 to −0.73).

**Confidence:** high for ≥180 s; moderate for the sign sensitivity at 30 s (n = 3 FAFB full model, n = 5 MCNS and FAFB-proxy; variants B and C not run).

### C02 — A total activity block of ≥10–20 s restores functional identity
**Verdict: supported** at 30 s; **conditional** at 180 s.

Evidence:
- DYNAMICS §4 (FAFB): FII_norm 1.001 at 30 s and 0.953 at 180 s; 2- and 5-s blocks fail; 10 s works at 30/180 s but not reliably at 360 s (0.428); a block started 5 or 20 s late still works.
- MECHANISM §2: 48/48 parameter sets rescued at 30 s; 34/48 = 0.71 [0.58, 0.83] at 180 s. The failures are re-ignition (11 sets; slow ATP recovery, τ_ATP PRCC +0.55, Holm p = 0.02) and structural loss (3 sets).
- REPLICATION: MCNS reset minus naive +0.948 [+0.916, +0.976] at 30 s; 5/5 recovered at 30 and 180 s.
- The block still rescues with AL-LN signs corrected: MCNS 0.989–1.000, FAFB-proxy 0.965–1.002.
- Closed-loop release at 11 s is insufficient: MCNS 0.926 vs 20-s block, Δ −0.071 [−0.081, −0.058], runaway 5/5; FAFB-proxy Δ −0.674 [−0.742, −0.624].

**Confidence:** high.

### C03 — No sharp percolation threshold under random loss
**Verdict: supported.** See H1. FAFB f_c 0.800, w 0.671; w > 0.25 in all four datasets.

**Confidence:** high.

### C04 — Hub- and demand-targeted loss is far more destructive than random loss
**Verdict: supported (direction); magnitudes differ.** See H2/H6. FAFB hub f_c 0.033 (fine grid), area-weighted 0.277, random 0.800. The rewired null is less fragile (hub f_c 0.266 vs 0.035). MANC, which is VNC-only, is hub-tolerant (f_c 0.388).

**Confidence:** high.

### C05 — Hub protection is far more efficient per neuron than whole-brain protection
**Verdict: supported and replicated.**
- DYNAMICS §7 (FAFB): top 1 % gives 0.578 → 0.703 [0.668, 0.743]; 0.090 vs 0.003 FII per 1,000 protected (~30×). Random sets give no gain.
- REPLICATION §2c (MCNS): top 1 % gives +0.201 [+0.127, +0.278], 84× per neuron vs whole CNS. Random 1 %: +0.010 [−0.089, +0.114].
- Useless without the block: all naive arms ran away.

**Confidence:** high (two connectomes).

### C06 — Graded reperfusion does not beat abrupt reperfusion plus block
**Verdict: supported** (k_exc ≤ 3); **provisional** at k_exc 10. See H3.

**Confidence:** high within model.

### C07 — Hypothermia widens the recoverability window
**Verdict: supported** (by construction of the injury rule). See H4.

**Confidence:** high that it holds; magnitude not a prediction.

### C08 — Simulated brain death requires near-total structural death
**Verdict: supported (model, FAFB).** DYNAMICS §3: 92 silent-and-unresponsive runs, all with ≥97.9 % of neurons dead. No silent-but-intact state occurred in 1,209 runs.

The index track found R_eff = 0 in every brain-death run (INDEX_INTERIM). Not tested in MCNS: the replication grid stopped at 360 s, which gave ~25 % dead.

**Confidence:** high for FAFB. This is a property of the model's injury rule and lack of other failure modes, not a biological finding.

### C09 — The runaway is ionic (K⁺) bistability, not synaptic reverberation
**Verdict: supported (model) for the maintained state. Entry after short anoxia is synaptically gated and sign-sensitive.**

MECHANISM §5 (default signs):
- All synapses off: 109/134 Hz (runaway maintained).
- K⁺ clamp or k_spk = 0: 2.4/2.3 Hz (collapse).
- Intact-network coexistence: silent 5/5; ≥9 mM kick gives 65–129 Hz.
- Hysteresis ≥13 mV.
- Excitatory synapses off cut the 30-s runaway to 14 Hz.

The sign-variant result (C11) shows that, under the default signs, the excitatory AL-LN loop is what carries the network into the ionic basin after 30 s. Bistability kicks under the corrected signs were planned (SY arm) but not run.

**Confidence:** high for the dissection; moderate for the generality of bistability across sign assignments.

### C10 — A 0.5 % driver-core block can replace a whole-brain block after short anoxia
**Verdict: not established (likely sign artefact) at 30 s; not supported at 180 s.**
- MECHANISM §4: under the default signs, the 0.5 % sustain set (696 neurons) rescues 5/5 held-out seeds (FII_norm 0.947).
- But 279 of the 696 sit in the AL, centred on the local-neuron types re-signed in C11.
- With those neurons inhibitory there is no high-rate 30-s runaway to target: FAFB variant A neurons ≥20 Hz = 0.05–0.28 % vs 64–70 % under the default signs (`sign_res/SN`, `raw/results/N`).
- Targeted blocks under corrected signs (ST arm) were not run.
- At 180 s, 20 % of the brain must be blocked (10 % fails 0/5).

**Confidence:** low that a small driver core exists in the real fly; the 180-s negative result is moderate.

### C11 — The AL driver core may be an artefact of transmitter sign calls
**Verdict: supported for short anoxia; refuted as an explanation of long-anoxia runaway.**

Evidence:
- Intact FAFB R_eff 19.6, falling to 7.4 with confident ACh/GABA/Glu calls only (`index/data/fafb_intact_index.json`: Reff_mf 19.62, Reff_mf_known 7.38).
- 30 s, AL-LN inhibitory: FAFB variant A 0/3 runaway; MCNS Δ FII_norm vs default signs +0.949 [+0.918, +0.977]; FAFB-proxy +0.714 [+0.632, +0.801].
- 180 s: MCNS Δ −0.991 [−1.119, −0.883]; FAFB-proxy −0.413 [−0.622, −0.192]. The runaway persists and becomes more continuous.

Preregistered falsification criterion: "runaway persists in sign-corrected variants with n ≥ 5". It is met at 180 s; at 30 s the high-rate state is absent. Variants B (low-confidence NT calls inhibitory) and C (no-call neurons silent) were not run.

**Confidence:** moderate.

### C12 — Runaway and block-reset findings replicate in MCNS
**Verdict: supported** (upgraded from n = 3 pilot to n = 5 full sweep). REPLICATION §3:
- (a) Runaway 5/5 at 30/180/360 s and 17/22/29 °C. Continuous at 30 s; episodic bursts at ≥180 s, so the 180-s naive FII_norm of 0.856 was measured in a burst gap.
- (b) Block rescue 5/5.
- (c) Hub protection 84×.
- (d) Temperature direction +0.185.
- Robust at 0.7× gain: 30-s naive still runaway 3/3.

BANC (146 result files on the Modal Volume) was not harvested or analysed.

**Confidence:** moderate (one replication dataset; MCNS is male and FAFB female, so sex and dataset are confounded).

### C13 — Connectome indices predict post-reset recovery but not naive runaway
**Verdict: supported (cross-validated; FAFB only; not written up as INDEX_RESULTS.md).** Source: `index/tables/model_results.json`, 1,118 of 1,209 dynamics runs (91 unprocessed), condition-grouped CV.
- **Post-reset FII_norm (non-runaway runs, n = 562):** excitatory bootstrap core size `ecore100` CV R² 0.964 [0.940, 0.978]; reciprocal-pair survival 0.961; simulator Reach 0.938.
  - Gain over dead fraction alone: ΔR² +0.066 [+0.025, +0.122].
  - Leave-arm-out transfer is weak for the ageing arm: ecore100 R² −0.50 on arm E; combined model 0.59.
- **Recovered vs not after reset (n = 582):** ecore100 CV AUC 0.999.
- **Runaway after naive reperfusion among viable networks (n = 482, 457 runaway):** best single-index AUC 0.53 (exc_in_frac); combined random forest CV AUC 0.61 [0.38, 0.76]. This is chance level: protocol and temperature, not wiring, decide entry.
- **Caveat:** indices were computed on default-sign graphs, which include the mis-signed AL local neurons.

**Confidence:** moderate.

### C14 — Metabolic arrest alone does not give hours-scale anoxia tolerance
**Verdict: supported (provisional extrapolation).** DYNAMICS §2: coma arrest (demand 5 %) changed the J slope from 0.996 to 0.886 per s. Extrapolated median-neuron LT50 is ≈1.51 h at J50 14,400 s, against a fly LT50 of ~8 h. Reaching ~8 h would need J50 ≈ 7.7×10⁴ s.

**Confidence:** low–moderate (extrapolation).

---

## 3. What was not done (explicit)

1. **BANC dynamic replication:** 146 result files exist on `flysim-data:/results_rep_v4/BANC*` but were never harvested or analysed (job ef7d8554 recorded as failed; cost unknown, ≤ $1.76).
2. **Sign variants B and C** and the FAFB variant-A arms beyond 30-s naive were not run. The missing variant-A arms are the intact baseline, reset, 180 s, bistability kicks and targeted blocks. Variant-A FII is therefore relative to the default-sign intact baseline (FII_norm ≈ FII/0.9940 = 0.920–0.942).
3. **Mechanism follow-up (605 runs) not executed:** 40/60-s blocks for all LHS sets, the full k_exc = 10 H3 grid, finer targeted-set brackets, and a −30 mV hysteresis ramp.
4. **Recoverability index:**
   - 91 of 1,209 runs not processed;
   - the BANC/MCNS/MANC recoverability map was only started (`cross_all.log`: FAFB random 0–0.6 only);
   - INDEX_RESULTS.md and index figures were not written.
5. **H2(a)** (death-time vs demand) was never computed.
6. **Preregistered seed count** (≥10 per condition) was not met (3–5).

---

## 4. Brain death and resuscitation: what this model can and cannot say

**Can say (as model predictions):**
- In a network built from a complete fly wiring diagram, failure to recover after anoxia can be a *state* rather than *damage*. A reversible ionic trap can block return of function while almost every neuron is alive.
- Clearing that trap with a time-limited, total activity block returns the original stimulus–response map, minus dead neurons. The block must outlast energy recovery: an 11-s K⁺/ATP-gated release re-ignited the trap.
- "Silent, unresponsive and unrecoverable" required ≥97.9 % neuron death in this model.
- Which neurons survive (hubs) matters more than how many.

**Cannot say:**
- Whether any of this occurs in living flies or mammals. The model has no glia beyond a clearance rate, no Ca²⁺/Na⁺ dynamics, no neuromodulation, no oxidative or inflammatory injury, no plasticity and no gene expression.
- Anything about clinical death by neurologic criteria, organ donation or patient care.
- Whether "original pathways persist" in a biological sense. FII measures a seven-probe stimulus → descending/motor map and is blind to lesions outside those pathways (32.5 % of the brain lost in the right optic lobe left FII ≈ 1.0).
- Whether the short-anoxia trap is real. It depends on transmitter-sign predictions for ~147 AL local neurons.
- Absolute times. Model seconds are compressed 50×, and injury constants are mammal-like, chosen for simulability.

Resuscitation of a "brain-dead" creature has not been demonstrated empirically. This work does not demonstrate it either.
