---
title: "A whole-brain Drosophila connectome model predicts a duration- and sign-dependent potassium trap after anoxia and a reset-based route to functional recovery"
author: "Srihari Muralikrishnan (Kalaris Labs; ORCID 0009-0004-0220-0361)"
correspondence: "hello@kalarislabs.com"
status: "Preprint v1.0 (2026-10-20). Not peer reviewed."
format: "bioRxiv / eLife-style research article"
---

# A whole-brain *Drosophila* connectome model predicts a duration- and sign-dependent potassium trap after anoxia and a reset-based route to functional recovery

**Srihari Muralikrishnan**¹

¹ Kalaris Labs (http://kalarislabs.com/). ORCID: 0009-0004-0220-0361. GitHub: @sriharideveloper.
Correspondence: hello@kalarislabs.com

**Article type:** Research article (computational neuroscience). **Status:** preprint v1.0, not peer reviewed.
**Claim type:** all results are predictions of a simplified model: a connectome-constrained leaky integrate-and-fire network with a phenomenological metabolic layer. The model was not fitted to living tissue, and no result is evidence about the reversibility of clinical brain death.

---

## Abstract

Clinical series and ex vivo perfusion experiments address whether a brain can be brought back after energy failure, and whether it returns as the same circuit. They cannot identify which network properties set the limit.

We built a whole-brain spiking model of the adult female *Drosophila* brain (FlyWire FAFB v783; 139,255 neurons, 3,732,460 connections). We added ATP, Na⁺/K⁺-pump, extracellular K⁺ and injury dynamics to a published leaky integrate-and-fire connectome model and tested six preregistered hypotheses. The analysis covered 3,487 simulated condition–seed runs; structural analyses used four connectomes and dynamic replication used the male CNS (MCNS). The model reproduced the anoxic spreading-depolarization latency of real flies (12.5 vs 14.2 s).

**Runaway after reoxygenation.** After 180 s of anoxia, abrupt reoxygenation drove the network into a self-sustained high-rate state in every parameter set tested (48/48), in MCNS, and under a corrected transmitter-sign assignment. The maintained state is ionic: it persisted with all synapses removed and vanished when activity-dependent K⁺ release was removed. Quiet and active states coexisted in the intact network, with hysteresis at least 13 mV wide.

**Sign dependence at short durations.** After 30 s of anoxia, the high-rate runaway depended on the transmitter sign of antennal-lobe local neurons that the default sign rule treats as excitatory (147 in FAFB, 219 in MCNS). With these neurons inhibitory, no runaway occurred in the FAFB model (0/3 seeds). The stimulus–response map was preserved in MCNS (identity index 1.000), leaving only a weak residual oscillation.

**Reset rescue.** A total activity block of 10–20 s after reperfusion restored the stimulus–response map at both durations and under both sign assignments. Releasing the block as soon as K⁺ and ATP had normalised (11 s) let the network re-ignite.

**Other results.** Graded reperfusion did not help. Remaining loss after the reset came from neuron death, and protecting the 1% most connected neurons was 30-fold (FAFB) to 84-fold (MCNS) more efficient per neuron than whole-brain protection. Silent, unresponsive "brain death" occurred only when ≥97.9% of neurons had died. Across four connectomes, targeted loss of hub or high-demand neurons cut sensorimotor reachability at a third or less of the damage required by random loss.

These predictions can be tested in flies with optogenetic silencing during reoxygenation.

**Keywords:** connectome; *Drosophila*; anoxia; spreading depolarization; reperfusion; bistability; extracellular potassium; brain death; resuscitation; neuroprotection

## Significance statement

It is unclear whether a brain that stops working after oxygen deprivation can be restarted as the same circuit, and what network property decides this. We simulated this question on a complete fly brain wiring diagram, using a spiking model with a simple energy and potassium layer. In the model, the main obstacle was usually not dead neurons. Reoxygenation could lock the network into a self-sustaining, potassium-driven hyperactive state. A brief, complete pause of activity released it and restored the original stimulus–response map, provided the pause outlasted energy recovery. After long oxygen loss this trap was robust. After short loss it depended on predicted neurotransmitter identities. These are model predictions and can be tested in flies with optogenetic silencing during reoxygenation.

---

## 1. Introduction

Death by neurologic criteria is defined clinically as the irreversible loss of all brain function, operationalised as coma, absent brainstem reflexes and apnoea after confounders are excluded [@greer2020; @greer2023]. Irreversibility is the contested term. Philosophical analyses index it to what current technology can reverse [@gligorov2023], and experimental work has repeatedly moved that boundary. Cat brains recovered neuronal function after periods of complete ischaemia longer than the classical few-minute rule allowed [@hossmann1973; @hossmann1973b]. In the BrainEx experiments, pig brains perfused from 4 h post-mortem regained cellular, vascular and spontaneous synaptic activity, yet showed no global electrocorticographic activity [@vrselja2019]. OrganEx extended cellular recovery to the whole body after 1 h of warm ischaemia [@andrijevic2022]. These results raise ethical questions [@farahany2019; @long2025] and a scientific one. If cells can be restored but organised network activity does not return, what network-level process prevents it, and is it reversible?

Mammalian data constrain the cellular side of this question. Electrical failure and ion-pump failure occur at distinct perfusion thresholds [@astrup1981; @hossmann1994]. Anoxic depolarization begins 2–5 min after severe ischaemia, and in dying human patients terminal spreading depolarization (SD) follows the final perfusion drop by a median of 3.9 min [@dreier2018]. SD is a continuum, from reversible waves in the penumbra to persistent depolarization in the core [@dreier2011; @dreier2015; @hartings2016]. Reperfusion can itself injure [@eltzschig2011; @choi1988], and in rodents gradual reperfusion and postconditioning reduce infarct size [@gao2008; @zhao2009]. Clinically, mild hypothermia after cardiac arrest improves neurological outcome [@haca2002; @bernard2002]. What is missing is a whole-network account of how these processes interact with the wiring that carries function.

Insects provide a natural comparison. Adult *Drosophila* enter anoxic coma within 30 s, survive several hours of anoxia and recover [@krishnan1997; @campbell2018; @callier2015]. Their brains generate SD with K⁺ surges similar to mammalian anoxic depolarization [@spong2016; @armstrong2011; @andersen2018], but in insects SD acts as a protective shutdown [@rodgers2010; @robertson2020]. The insect CNS still has an irreversible end state. Repeated ouabain comas end in a terminal extracellular K⁺ plateau and death [@armstrong2011], and repetitive anoxia at warm temperatures degrades K⁺ homeostasis [@rodriguez2012]. SD temperature dependence has been measured [@andersen2025], as has the role of pump kinetics [@andersen2022] and glial K⁺ channels [@andersen2026].

Connectomes now make it possible to ask these questions of a complete nervous system. The FlyWire reconstruction of the adult female brain [@zheng2018; @dorkenwald2024; @schlegel2024] has been turned into a leaky integrate-and-fire (LIF) model with a single free parameter. That model predicted sensorimotor pathways that were then validated in flies [@shiu2024]. Connectome-constrained models also predict activity in the fly visual system [@lappalainen2024]. Further connectomes exist for the male ventral nerve cord (MANC; [@takemura2024; @marin2024; @cheong2026]), the female brain and nerve cord (BANC; [@bates2026]) and the male CNS (MCNS; [@berg2026]). They allow the same questions to be asked of independent wiring diagrams.

Here we extend the Shiu et al. model with a phenomenological metabolic and ionic layer and use it as a proxy brain for anoxia, "brain death" and resuscitation. Before running any insult simulation, we preregistered six hypotheses:
- H1: a sharp percolation-like damage threshold.
- H2: targeted vulnerability of high-demand and hub neurons.
- H3: a benefit of graded over abrupt reperfusion.
- H4: widening of the recoverable window by hypothermia in proportion to Q10.
- H5: loss of functional identity before loss of activity or structure.
- H6: replication across connectomes.

Three of these were falsified (H1, H3, and H5 as preregistered); H2 and H6 were partly supported and H4 supported.

The main finding was not anticipated. Reoxygenation can place the model in a self-sustained, K⁺-coupled attractor, and a timed activity reset releases the network from it. Recovery is then limited only by neuron death. After long anoxia the trap is robust to parameters, connectome and transmitter-sign assignment. After short anoxia it depends on the sign assigned to a small set of antennal-lobe neurons.

---

## 2. Results

### 2.1 A calibrated whole-brain LIF model with an energy and ion layer

We reimplemented the Shiu et al. LIF model [@shiu2024] in NumPy/SciPy (V_rest = V_reset = −52 mV, V_th = −45 mV, τ_m = 20 ms, τ_syn = 5 ms, refractory 2.2 ms, delay 1.8 ms). It ran on the Princeton FAFB v783 connection table (≥5 synapses per pair; 3,732,460 edges, 50,666,648 synapses). Signs follow the presynaptic predicted transmitter: GABA and glutamate inhibitory; acetylcholine, monoamines and unknown transmitters excitatory. We call this the default sign rule.

**Weight calibration.** With the published synaptic weight w_syn = 0.275 mV this matrix was supercritical. After a 150-ms sugar-receptor stimulus the active population kept growing (538 → 1,805 → 2,721 neurons at 0–150, 200–400 and 1,000–1,200 ms), whereas Shiu's own matrix returned to silence (326 → 33 → 0). We therefore used w_syn = 0.185 mV. This gave 344 → 0 → 0 active neurons and an MN9 response of 107 Hz (97 Hz on Shiu's matrix) to 150-Hz sugar drive (Fig. 1A shows the response curves at the published weight). At this weight:
- sugar input at 50/100/200 Hz drove MN9 at 3.3/117.7/148.3 Hz;
- LC4 drove the giant fibre DNp01 (L/R) at 213/229 Hz;
- LC9 drove DNg100 at 116 Hz.

**Slow layer.** The added layer (Methods) tracks:
- per-neuron ATP and pump-dependent depolarization;
- extracellular K⁺ in 79 neuropil compartments coupled by a geometric contact adjacency;
- regenerative K⁺ release (SD) and depolarization block;
- a Ca²⁺-load proxy driving permanent death;
- reversible synaptic failure.

At full supply the layer left the LIF unchanged. With activity-dependent K⁺ release off, spike counts were bit-identical to the pure LIF model; with it on, the rate correlation was 0.999 (MN9 141 vs 140 Hz).

**Calibration against fly data.** After calibration (τ_ATP = 30 s, J50 = 600 s; Methods), anoxic SD began 12.5 model-s after anoxia onset at 22 °C. Wild-type flies show anoxic SD at 14.2 ± 0.3 s [@andersen2018]. Peak extracellular K⁺ during 30–60 s of anoxia was 45.2 mM, inside the 44–58 mM range recorded in fly brain [@andersen2018; @armstrong2011]. Fig. 1B shows the layer's behaviour at the pre-calibration constants. SD onset scaled with τ_ATP (9.0–10.0 s at 20 s; 14.0–16.5 s at 40 s). It moved by only 2–3 s across a fourfold change in the slow-time compression factor (Methods).

![Figure 1](figures/Figure1_model_validation.png)

**Figure 1. Model validation.** (A) Sensorimotor responses of the LIF model on the Princeton ≥5-synapse FAFB v783 matrix (blue) and on Shiu et al.'s connectivity (grey), at the published w_syn = 0.275 mV (5 trials × 1 s; mean ± s.d.). a, sugar GRNs → MN9; b, LC4 → giant fibre (DNp01); c, LC9 → walking DNs; d, turn presets lateralize DNa02. The calibrated w_syn = 0.185 mV used for all insult simulations lowers these curves (values in text). (B) Metabolic-layer unit tests at the default (pre-calibration) constants, speedup 50. a, raster: anoxia (S = 0, grey) silences the brain and naive reoxygenation re-ignites it; b, population rate with and without a 20-s post-reperfusion block; c, mean ATP, depolarized fraction and synaptic release; d, extracellular K⁺ per neuropil (no reset); e, f, focal ischaemia of the right optic lobe: K⁺ time courses coloured by distance to the core, and peak K⁺ vs centroid distance.

### 2.2 Metabolic vulnerability targets hubs, and the effect replicates across four connectomes

We first asked whether connectome structure alone has a damage threshold. We measured strict sensory-to-sink reachability as interneurons were removed (Fig. 2A). Reachability is the fraction of (sensory neuron, descending/motor neuron) pairs joined by a ≤6-hop path through surviving interneurons. Baseline reach was 0.8687 (FAFB), 0.7291 (BANC), 0.9019 (MCNS) and 0.9754 (MANC).

**Random loss was graceful in every dataset.** Reach stayed ≥0.8 of baseline up to f = 0.50 (FAFB), 0.45 (BANC), 0.65 (MCNS) and 0.85 (MANC). The logistic transition width (FAFB w = 0.671) exceeded the preregistered falsification bound of 0.25, so H1 was falsified for random structural loss in all four datasets. The verdict held at a ≥10-synapse edge cut-off (w = 0.960).

**Loss weighted by metabolic demand was far more destructive.** In FAFB:
- removing neurons with probability proportional to membrane area gave f_c = 0.277 (ratio to random 0.346 [95% CI 0.337, 0.356]);
- synapse-count weighting gave f_c = 0.228 (ratio 0.285 [0.278, 0.291]);
- hub-targeted removal halved reach after removing 7.3% of interneurons (fine-grid f_c 0.033; ratio 0.041 [0.041, 0.041]) (Fig. 2B).

The preregistered H2(b) criterion (ratio ≤0.7, CI excluding 1) was met in all four connectomes for synapse-weighted and hub-targeted loss.

Two controls located the effect:
- Degree-bin-matched random removal reproduced the hub-targeted curve (ratio 0.044 [0.044, 0.045]). Vulnerability is therefore carried by degree, not by the identity of particular hubs.
- A degree-preserving rewired null was much less fragile (hub-targeted f_c 0.266 vs 0.035 in the real graph). Real sensorimotor routes are therefore concentrated on hubs beyond what the degree sequence implies.

Magnitudes did not replicate. The f_c max/min ratio across datasets was 6.32 for synapse-weighted and 23.35 for hub-targeted loss, because the VNC-only MANC was much more hub-tolerant (f_c 0.388).

**No single cell type acted as a structural keystone.** Removing every neuron of any one evaluated type reduced reach by at most 0.88% (FAFB, CB1231), 1.03% (BANC), 0.86% (MCNS) and 0.00% (MANC) (Fig. S2). Graded type-level importance did replicate (Spearman ρ of influence loss 0.65–0.91 across dataset pairs).

Under random loss, weighted influence and the giant strongly connected component fell to 50% before binary reach did, in all four datasets (FAFB D50 0.497 and 0.506 vs 0.759; Fig. S1). This is a structural precondition for a network that keeps its routes but processes signals differently.

![Figure 2](figures/Figure2_structural_vulnerability.png)

**Figure 2. Structural vulnerability across four connectomes (FAFB, BANC, MCNS, MANC).** (A) Strict sensory→sink reach (≤6 hops) vs fraction of interneurons removed under random, metabolic (area- or synapse-weighted), hub-targeted, degree-matched and focal regimes; band = 2.5–97.5 percentile across removal orders. (B) Targeted/random ratio of logistic midpoints and D50 with 95% bootstrap CIs; dashed line = preregistered H2 cut-off 0.7.

### 2.3 Reoxygenation drives the network into a self-sustained runaway state

**Protocol.** We simulated whole-brain anoxia of duration D at 17–33 °C (1,209 condition–seed runs, 274 conditions; default sign rule). Each run was followed by:
- 60 s of reperfusion protocol;
- 20 s of free observation;
- a seven-probe stimulus battery read out at all 1,415 descending and motor neurons.

**Identity metric.** Functional identity was scored by the preregistered functional identity index (FII). FII is the Pearson correlation of the post-protocol stimulus × readout response matrix with the intact baseline, normalised by the noise ceiling (0.9940). FII is insensitive to gain: a 0.5×-scaled baseline gives r = 1.000, and a readout-shuffled baseline gives −0.220.

**Naive reoxygenation failed at every insult.** At 22 °C, abrupt ("naive") reoxygenation ignited a self-sustained high-activity state in all 30 runs from 7 to 480 s of anoxia (6 durations × 5 seeds; Fig. 3B). This included 7 s, which is shorter than SD onset. In this state:
- 13–49% of neurons fired spontaneously;
- late firing averaged 105–119 Hz per neuron;
- FII_norm was ≈0.

The runaway appeared at time-compression factors 25, 50 and 100 (FII_norm −0.10, −0.06, 0.065 after 30 s), so it is not a compression artefact. Sham protocols gave FII_norm 1.001. At 17 °C, 7 s of anoxia did not ignite it (FII_norm 0.988). Section 2.6 shows that the short-duration part of this result depends on transmitter-sign assignment.

**Two boundaries (Fig. 3A).**
- *Functional boundary.* Under the default signs it lies at ~7 s of anoxia, beyond which naive reperfusion produces runaway.
- *Structural boundary.* Even when the runaway was prevented (Section 2.4), neuron death removed identity beyond a recoverable duration D_max. D_max was 359 [342, 381], 271 [261, 278], 149 [145, 153] and 124 [112, 129] model-s at 17, 22, 29 and 33 °C.

**Brain death required near-total neuron death.** Silent, unresponsive "brain death" was defined as a spontaneously active fraction <0.01 and no readout response >8 Hz after the protocol. It occurred in 92 runs, always with ≥97.9% of neurons dead. The model never produced a silent but structurally intact unresponsive state: every surviving network either responded to probes or ran away. In this model, therefore, irreversible functional death required near-total structural death.

**Other insults.**
- Partial supply at 25% of normal was tolerated indefinitely (960 s: no SD, no deaths).
- At 10% supply, SD began at 17 s and 97.9% of neurons were dead at 960 s.
- Hyperthermia without hypoxia (33/37/40 °C for 300 s) had no effect.
- Focal ischaemia of the right optic lobe killed 32.5% of the brain at 960 s but left FII_norm at 0.997–1.000, because the probe battery barely uses that lobe. Identity in this model is defined relative to the battery.

![Figure 3](figures/Figure3_brain_death_map_runaway.png)

**Figure 3. Brain-death map and the reperfusion runaway (calibrated set; default sign rule; speedup 50).** (A) FII_norm with naive reperfusion and with a 20-s block reset, and dead fraction, by anoxia duration × temperature; letters give the modal outcome class (n = 5 seeds per cell). (B) 180-s anoxia, seed 0: naive reperfusion (red) vs a 20-s total activity block at reperfusion (blue).

### 2.4 A brief total activity reset rescues function, with a therapeutic window and temperature dependence

**Protocol comparison.** We compared 32 reperfusion protocols at 30, 180 and 360 s of anoxia (5 seeds each; Fig. 4A). Every protocol without a total activity block of ≥10 s failed at all three durations. The failures included:
- −10 and −20 mV hyperpolarization;
- two- to fourfold GABAergic gain;
- SD-release blockade and threefold K⁺ clearance;
- ramped (10/30/60 s) or pulsed (3×10 s, 6×5 s) reperfusion;
- cooling to 12 or 17 °C during reperfusion;
- low-level noise and spike-frequency adaptation (the best of these, at FII_norm 0.73).

**The block.** A 20-s block of all spiking and synaptic transmission at reperfusion is the model analogue of a pharmacological isoelectric state. It recovered FII_norm 1.001, 0.953 and 0.578 at 30, 180 and 360 s.
- Blocks of 2 and 5 s failed.
- A 10-s block worked at 30 and 180 s but not reliably at 360 s (0.428).
- Blocks of 20–40 s worked at every duration.
- The block also worked when started 5 or 20 s after reperfusion, once the runaway was established (FII_norm 0.841–0.846 at 360 s, the same as an immediate block). The attractor can therefore be erased after it has formed.

**After the reset, the remaining loss was structural.** At 360 s, 24.8% of neurons were dead and FII_norm was 0.578. A counterfactual that revived all dead neurons before the reset gave 1.001 at every duration. Adding clearance, noise, adaptation, SD blockade or staged release to the block changed FII_norm by ≤0.01. Staged reawakening (sensory → central → descending, or the reverse) made no difference (0.846 vs 0.844).

**H3 was falsified, with the direction reversed.** At 360 s, ramped and pulsed reperfusion gave lower FII_norm than abrupt reperfusion:
- ramp 30 s: Δ = −0.331 [−0.583, −0.130];
- ramp 60 s: −0.297 [−0.540, −0.092];
- pulsed 3×10 s: −0.221 [−0.460, −0.011].

They also caused more deaths (ramp 60 s +6.6 percentage points), because slower re-energization prolongs depolarization. The baseline model has no spike-driven injury term, so we added one in a follow-up (Section 2.7; Fig. S3). With it, no graded protocol beat abrupt reperfusion plus a block at any tested gain or duration. Without a block, pulsed reperfusion reduced death slightly at the highest calibrated gain (−1.6 pp [−1.8, −1.4] at 180 s), but the network stayed in runaway.

**H4 was supported.** Cooling widened the recoverable window. The D_max ratio 17 °C/29 °C had CI [2.29, 2.57], giving Q10_window = 2.09 [1.99, 2.19] at model Q10 2.5. This lies inside the preregistered "approximately Q10" band. The same pattern held at model Q10 2.0 (1.71) and 3.0 (2.44), with each value close to Q10^0.8. That exponent is what the Hill-4 injury hazard predicts (Methods). The widening is therefore set by the damage rule, as the preregistration warned, and not by network nonlinearity (Fig. 4B).

Timing of cooling mattered:
- Cooling to 12 °C before and during the insult gave 0.6% dead and FII_norm 0.977 after 360 s, the strongest tolerance intervention tested.
- Cooling applied only during reperfusion was harmful (runaway; FII_norm −0.30 to −0.13).
- Reducing demand to 20% delayed SD to 33 s and gave FII_norm 0.704 at 360 s.
- Preconditioning (15 s anoxia plus reset 60 s earlier) gave no protection, as expected for a model without gene expression or plasticity.

**H5 was not supported as preregistered.** With the reset, the D50 values of identity, evoked activity and reach coincided: FII_norm 374 s [363, 385]; evoked-neuron count 368 s [365, 373]; reach 390 s [388, 391]. The "resuscitated but different" region (active, connected, FII_norm ≤0.5) covered 0.0% of non-runaway runs. A dissociation between identity and activity did exist, covering 22.1% of runs, all of them in the runaway state. In this model, "resuscitated but different" exists only as the pathological attractor; once that is cleared, identity tracks structure. Pooled over 372 non-runaway reset runs, FII declined approximately linearly with dead fraction (logistic ΔAIC +0.11 vs linear; width 1.56), which also falsifies H1 for dynamical damage.

![Figure 4](figures/Figure4_resuscitation_temperature.png)

**Figure 4. Resuscitation protocols and temperature.** (A) FII_norm per protocol after 30-, 180- and 360-s anoxia (n = 5 seeds; 95% bootstrap CI). (B) Left: D_max vs temperature with Q10 scaling predictions. Right: FII_norm vs dead fraction for non-runaway reset runs (n = 372).

### 2.5 The maintained runaway is an ionic attractor with bistability and hysteresis

**Dissection.** To find what sustains the runaway, we applied manipulations from 20 to 60 s after abrupt reperfusion (5 seeds; Fig. 5A-c). Rates below are given after 30 and 180 s of anoxia; unmanipulated rates were 105 and 121 Hz per neuron.
- Removing all synaptic transmission did not stop it (109 and 134 Hz per neuron during the manipulation).
- Clamping extracellular K⁺ at 5 mM collapsed it (2.4 / 2.3 Hz), as did removing spike-driven K⁺ release (2.4 / 2.3 Hz).
- Removing regenerative SD release did not abolish it.
- Synaptic excitation contributed after short anoxia: removing excitatory synapses cut the 30-s runaway to 14 Hz.
- Releasing the K⁺ clamp left a third, low-rate (≈4 Hz) synaptic reverberation with FII_norm 0.71–0.72.

**Bistability in the intact network** (Fig. 5A-a,b). The same mechanism produced bistability in a network that had never been insulted. With identical parameters, the unperturbed network stayed silent (5/5 seeds; K⁺ 5.0 mM). A 2-s K⁺ clamp to ≥9 mM drove it into a self-sustained state (5/5; 65–129 Hz per neuron; K⁺ ≈15 mM).
- The minimal ignition kick tracked the analytic single-neuron K⁺ threshold, K_th = K0·exp((V_th − V_rest − δ)/(g_K·RT/F)). For tonic offsets δ ≥ 0 it sat one grid step above it (simulated 6–9 mM vs analytic 5.3–7.9 mM for δ = +6 to 0 mV).
- Synaptic gain between 0.5× and 1.5× did not move the boundary.
- Spike kicks without a K⁺ rise produced only the low-rate state.

**Hysteresis.** Ramping δ upward ignited the silent network at +7.15 mV in 3/3 seeds, matching V_th − V_rest = 7 mV. Ramping back down, the active state persisted to −6 mV (252–276 Hz per neuron), sustained by falling ATP (0.94 → 0.73) and rising K⁺ (9.5 → 20 mM). The hysteresis loop is therefore at least 13 mV wide; its lower edge was not reached. A total block quenched a kicked runaway only if it lasted ≥10 model-s at δ = 0 (1, 2 and 5 s failed) or ≥20 s at δ = +4 mV. These durations match those needed after anoxia.

**Interpretation within the model.** These results fit a positive feedback loop: spiking releases K⁺, K⁺ depolarizes neurons through the Nernst term, and depolarization raises spiking. The loop creates two stable states, separated by a K⁺ threshold set by the gap between rest and spike threshold:
- quiet (K⁺ 5 mM);
- active (K⁺ ≈10–20 mM).

Anoxia does not create the active state; on reoxygenation it supplies the perturbation. As neurons repolarize, a K⁺ load of 45–80 mM passes through the window between K_th and depolarization block. Whether that perturbation lands in the active basin also depends on synaptic recruitment after short anoxia (Section 2.6). The reset holds the network silent until K⁺ and ATP relax back into the quiet basin.

**Robustness to parameters.** We varied ten uncertain constants over wide ranges in 48 Latin-hypercube parameter sets (784 runs; Fig. 5B). Abrupt reperfusion produced runaway in 48/48 sets (288/288 runs), of the high-rate type in 41/48. The 20-s block rescued 48/48 sets after 30 s of anoxia but only 34/48 (0.71 [0.58, 0.83]) after 180 s.
- Eleven failing sets re-ignited after the block ended. These had slow ATP recovery (median τ_ATP 44 vs 26 s; PRCC +0.55, Holm p = 0.02) and low reserve.
- Three failed through neuron death (J50 ≈ 325 s).

The existence of the attractor is robust across the explored parameter box, but the adequacy of a fixed 20-s reset is not. We therefore tested a closed-loop release in MCNS and the FAFB-proxy (Section 2.9): it ended the block when the maximum compartment K⁺ fell below 6 mM and mean ATP exceeded 0.9. The criterion was met after 11 s in every run, and the network re-ignited (Section 2.9). Normalised K⁺ and ATP are therefore not a sufficient release signal in this model. The planned 40/60-s blocks for all 48 parameter sets were not run (Modal spend limit; Supplement S9).

![Figure 5](figures/Figure5_ionic_bistability_robustness.png)

**Figure 5. Ionic bistability and robustness (FAFB, default sign rule).** (A) a, intact network: state 30–40 s after a 2-s K⁺ clamp kick vs tonic offset δ (5 seeds at δ = 0, 3 otherwise); blue line = analytic intrinsic threshold K_th; b, slow ramp of δ, −6 → +10 → −6 mV over 2 × 160 s (3 seeds); c, manipulations applied 20–60 s after abrupt reperfusion (n = 5); d, quench of a kicked runaway by a total block of duration B (n = 3). (B) a, set-mean FII_norm in 48 Latin-hypercube parameter sets × 3 seeds (diamond = calibrated set; dashed line = 0.8); b, PRCC of set-level outcomes with each parameter (* Holm p < 0.05); c, outcome of the 20-s block after 180-s anoxia in the (τ_ATP, reserve) plane.

### 2.6 Transmitter-sign robustness: short-anoxia runaway depends on antennal-lobe local-neuron sign

The default sign rule maps unknown and monoaminergic transmitters to +1. Many antennal-lobe (AL) local neurons are GABAergic in flies, yet in FAFB they enter the model as excitatory. A connectome index computed on the intact FAFB graph flagged this. The mean-field excitatory reproduction number R_eff was 19.6 with all neurons and fell to 7.4 when only confidently predicted ACh/GABA/Glu neurons were counted (Section 2.10).

**Variant tested.** We set every neuron of class ALLN or of type lLN/vLN/v2LN/ilLN/il2LN to inhibitory. In FAFB this mask contains 431 neurons. Of these, 284 were already inhibitory and 147 changed sign (90 predicted cholinergic, 57 predicted monoaminergic). The same rule flipped 147 neurons in the FAFB-proxy model and 219 in MCNS. We tested the corrected signs:
- in the FAFB full model (variant A; 30-s anoxia, naive reperfusion; 3 seeds, run locally);
- in MCNS and the FAFB-proxy model (30 and 180 s, naive and reset; 5 seeds each; Table 1, Fig. 8c).

**Table 1. Effect of setting antennal-lobe local neurons inhibitory (22 °C).** FII_norm with 95% bootstrap CI where n = 5. Late rate = mean firing in the last part of free observation (Hz per neuron; `obs_late_Hz`, the metric used by the runaway rule). For the FAFB default-sign row it is the dynamics-track arm-C value (n = 5); the separate anatomy runs report a different metric, the mean rate over the whole observation window (91.3–104.3 Hz per neuron; Supplement Table S6h-1). Runaway class = late rate > 0.5 Hz per neuron or spontaneous fraction > 0.05 (preregistered rule).

| model | anoxia | protocol | default signs: FII_norm (late rate; runaway) | AL local neurons inhibitory: FII_norm (late rate; runaway) | Δ FII_norm [95% CI] |
|---|---|---|---|---|---|
| FAFB full model (variant A) | 30 s | naive | −0.058 [−0.180, 0.152] (111.0 Hz [103.2, 119.2]; 5/5) | FII 0.915–0.936 (0.035–0.17 Hz; **0/3**) | not computed (n = 3) |
| MCNS | 30 s | naive | 0.050 [0.022, 0.078] (90.9 Hz; 5/5) | **1.000 [0.996, 1.002]** (3.0 Hz; 5/5) | +0.949 [+0.918, +0.977] |
| FAFB-proxy | 30 s | naive | 0.146 [0.057, 0.228] (47.1 Hz; 5/5) | 0.859 [0.852, 0.867] (0.82 Hz; 5/5) | +0.714 [+0.632, +0.801] |
| MCNS | 180 s | naive | 0.856 [0.810, 0.901] (127.8 Hz, episodic; 5/5) | −0.135 [−0.234, −0.041] (120.5 Hz, FSA 0.897; 5/5) | −0.991 [−1.119, −0.883] |
| FAFB-proxy | 180 s | naive | 0.288 [0.086, 0.505] (124.1 Hz; 5/5) | −0.125 [−0.155, −0.093] (126.6 Hz; 5/5) | −0.413 [−0.622, −0.192] |
| MCNS | 30 / 180 s | 20-s block | 0.999 / 0.997 (recovered 5/5) | 1.000 / 0.989 (recovered 5/5) | — |
| FAFB-proxy | 30 / 180 s | 20-s block | 1.002 / 0.971 (recovered 5/5) | 1.002 / 0.965 (recovered 5/5) | — |

**Short anoxia (30 s): the high-rate runaway and the identity loss disappeared.**
- In the FAFB full model all three seeds stayed below both runaway thresholds (spontaneous fraction ≤0.005).
- Only 0.05–0.28% of neurons fired ≥20 Hz, against 64–70% under the default signs.
- Raw FII was 0.915–0.936. It was computed against the default-sign intact baseline, because no variant-A baseline was run.
- In MCNS and the FAFB-proxy, identity was preserved, but a weak oscillation persisted (3.0 and 0.82 Hz per neuron). The preregistered classifier still counts this as runaway, so the runaway did not vanish entirely in these two models.

**Long anoxia (180 s): the runaway persisted under corrected signs.** It became continuous rather than episodic in MCNS (spontaneous fraction 0.897), and the total block still rescued.

**Interpretation within the model.** After short anoxia, the excitatory AL local-neuron loop supplies the synaptic recruitment that carries the network into the ionic basin. After long anoxia, the K⁺ load alone is enough.

Consequences for other results:
- Under corrected signs, the duration at which naive reperfusion first produces runaway lies between 30 and 180 s. It was not localised further.
- The "~7 s" functional boundary of Section 2.3 holds only for the default signs.
- The driver-core result of Section 2.7 inherits this caveat.
- Sign variants assigning inhibition to all low-confidence calls (variant B) or silencing neurons without a call (variant C) were prepared but not run.

### 2.7 A recurrent excitatory driver core and targeted reset (default signs)

**Ignition.** After 30 s of anoxia, ignition was near-synchronous: population rate rose from ~0 to ~100 Hz per neuron within ≈1.5 model-s, starting ≈4 model-s after reoxygenation. Only 13–20% of neurons were intrinsically suprathreshold at the end of reperfusion, so most firing was recruited (Fig. 6a).

**Driver set.** We ranked excitatory neurons by their contribution to sustaining the state (rate × excitatory synapses onto active neurons) using training seeds, and tested targeted 20-s blocks on held-out seeds. The top 0.5% "sustain" set (696 neurons) formed a single strongly connected excitatory core. 96% of it lay in one strongly connected component (31 excitatory edges per node), against 40% for top-degree hubs and 0.4% for random sets of the same size.

**Targeted blocks.**
- After 30 s of anoxia, blocking this 0.5% set rescued 5/5 held-out seeds (FII_norm 0.947 [0.943, 0.951]). Hubs needed 2% (2,785 neurons), and random sets failed up to 20%.
- After 180 s of anoxia, 20% of the brain had to be blocked (bracket 10–20%), consistent with near-ceiling K⁺ everywhere.

**This result is not established for the real fly.** 279 of the 696 driver neurons sit in the antennal lobes, centred on local neurons lLN1_bc and lLN2P and on DA1/VA1v/VM5d/DL2v projection neurons. These are the populations re-signed in Section 2.6, and with the corrected signs there was no high-rate 30-s runaway to target. Targeted blocks under the corrected signs were not run.

**Excitotoxic reperfusion term** (Fig. S3). We added an injury term that accrues only when neurons fire above 50 Hz in raised K⁺. It increased death under abrupt reperfusion: after 180 s of anoxia, 0.9%, 1.1%, 3.5% and 39.5% dead at gains 0.3, 1, 3 and 10, against 0.9–1.0% with the block. The block's protective margin grew with gain (Δ dead −38.6 pp at gain 10). In this model the reset is also the best protection against excitotoxic injury, because it removes the surge that causes it. The gain-10 result covers only the abrupt-vs-block comparison and is provisional.

![Figure 6](figures/Figure6_driver_core_targeted_reset.png)

**Figure 6. Attractor anatomy and targeted reset (FAFB, default sign rule; see Section 2.6 for the sign caveat).** (a) Population rate after reoxygenation (5 seeds thin, mean thick). (b, c) FII_norm on held-out seeds 10–14 vs size of the set blocked for 20 s from reperfusion onset; black = whole-brain block, dotted = no block. (d) Fraction of each set in its largest strongly connected component of the excitatory subgraph.

### 2.8 Hub neuroprotection and ageing

Remaining loss after the reset was structural, and structural vulnerability concentrates on hubs (Section 2.2). We therefore tested protecting neuron subsets from injury during 360 s of anoxia followed by the reset (FAFB; Fig. 7b).
- Protecting the top 1% of neurons by degree (1,393 neurons) raised FII_norm from 0.578 to 0.703 [0.668, 0.743].
- Whole-brain protection reached 1.001, so per protected neuron hub protection was about 30× more efficient (0.090 vs 0.003 FII per 1,000 neurons).
- The top 10% reached 0.906, above the recovery criterion, although 22.2% of neurons still died.
- Random sets of the same sizes gave no gain (0.568–0.590), even though random 10% prevented a similar number of deaths. What mattered was which neurons were saved.
- Without the reset, all protected arms stayed in runaway (FII_norm −0.25 to 0.34).

MCNS replicated this effect (Section 2.9).

**Ageing.** An ageing parameter combined synapse loss with reduced energy reserve and K⁺ clearance. It lowered FII_norm roughly linearly with synapse loss (Fig. 7a).
- Losing synapses in proportion to demand was more damaging than random loss at every level (30% lost: 0.668 vs 0.745; 15%: 0.803 vs 0.856).
- Demand-weighted ageing crossed the 0.8 recovery line at ~15% synapse loss, and random ageing at ~25%.
- Partial supply at 30% triggered SD only at the oldest level.
- Boosting reserve or clearance, or mild hypometabolism, did not restore identity in the oldest network (all 0.668). The age deficit in this model is in the wiring, and metabolic support does not reach it.

![Figure 7](figures/Figure7_hub_protection_ageing.png)

**Figure 7. Ageing and targeted neuroprotection (FAFB).** (a) FII_norm vs fraction of synapses lost, demand-weighted vs random (n = 3 seeds per point). (b) FII_norm after 360-s anoxia with injury protection of top-degree vs random neuron sets of 1–10%, and whole brain.

### 2.9 Dynamic replication in the male CNS connectome

We rebuilt the model on MCNS v1.0 (166,700 neurons, 6,242,085 edges at ≥5 synapses). w_syn was set at the same distance below the reverberation threshold as in FAFB (0.0538 mV). Two substitutions were needed because MCNS lacks the corresponding FAFB data:
- demand was proxied by total synapse count;
- K⁺ coupling used a connection-based Jaccard 9-nearest-neighbour neuropil graph.

To separate dataset from proxy effects, FAFB was re-run with the same proxies ("FAFB-proxy"). Identity used the responsive probes; sugar → MN9 is silent in MCNS and was excluded. The analysis covered 211 runs, with 5 seeds per core condition (Fig. 8).

**All four tested FAFB dynamic results replicated in MCNS:**
- *Runaway.* Naive reperfusion produced runaway in 5/5 runs at 30, 180 and 360 s and at 17, 22 and 29 °C. At 30 s, FII_norm was 0.050 [0.022, 0.078] with a late rate of 90.9 Hz per neuron. After ≥180 s the MCNS runaway was episodic: all-region bursts of up to ~410 Hz per neuron were separated by silent gaps. The 180-s naive FII_norm (0.856) was taken in such a gap and understates the pathology; the VNC motor pattern was lost (VNC FII_norm 0.045).
- *Block rescue.* The 20-s block recovered 5/5 runs at 30 and 180 s. Reset minus naive at 30 s was +0.948 [+0.916, +0.976]. Sham runs were unaffected (1.001).
- *Hub protection.* After 360 s plus block, protecting the top 1% (1,667 neurons) gained +0.201 [+0.127, +0.278], 84× more efficient per neuron than whole-CNS protection. Random 1% gained +0.010 [−0.089, +0.114].
- *Temperature.* After 180 s plus block, 17 °C minus 29 °C was +0.185 [+0.123, +0.260].

**Robustness and new observations.**
- The 30-s runaway persisted at 0.7× the reverberation-threshold gain (3/3).
- After naive reoxygenation the VNC was the first region to exceed 1 Hz per neuron (1.5 s, vs 2.5 s for the central brain and 4.0 s for the optic lobe). It was also the most persistent part of the 30-s attractor (spontaneous fraction 0.43).

**Closed-loop release was insufficient.** Ending the block when K⁺ and ATP had normalised (11 s) left residual activity. In MCNS: 11.3 Hz per neuron, runaway 5/5, FII_norm 0.926; Δ vs 20-s block −0.071 [−0.081, −0.058]. In the FAFB-proxy: 40.2 Hz per neuron, FII_norm 0.296; Δ −0.674 [−0.742, −0.624]. Some slow state not captured by K⁺ and ATP must still relax before release; which one was not tested.

**Not analysed.** A BANC sweep (146 result files) completed on the compute volume but was not retrieved or analysed. MANC was not simulated dynamically. MCNS is male and FAFB female, so dataset and sex are confounded.

![Figure 8](figures/Figure8_mcns_replication_sign.png)

**Figure 8. Cross-connectome replication, sign correction and release rule.** FII_norm on responsive probes; dots = runs (n = 5 seeds; FAFB reference arms in b used 3), bars = mean with 95% bootstrap CI; dashed line = recovered threshold 0.8. Grey = FAFB (dynamics track: geometric K⁺ graph, area demand); lavender = FAFB re-run with the CNS-model proxies; orange = MCNS. (a) 30- and 180-s anoxia, naive vs 20-s block. (b) 360-s anoxia + block with injury protection of random or top-degree sets. (c) Left pair: 30-s naive reperfusion with the default signs ("Shiu") vs antennal-lobe local neurons inhibitory ("AL-LN −1"). Right pair: 180-s anoxia, closed-loop release (K⁺ < 6 mM and ATP > 0.9; 11 s) vs fixed 20-s block.

### 2.10 A connectome-only recoverability index

We asked whether statistics computed on the surviving wiring alone, without simulation, predict the dynamical outcome. We computed 19 indices on the surviving FAFB graph of 1,118 of the 1,209 dynamics runs (91 runs unprocessed). Prediction was evaluated with condition-grouped 5-fold cross-validation (4–20 repeats depending on the model) and cluster-bootstrap CIs. All indices used the default sign rule.

**Intact network.** The intact FAFB excitatory network is far above criticality. Using the model's own neuron parameters, its mean-field reproduction number (extra spikes per spike) is R_eff ≈ 19.6, or 17.0 with inhibition included. Counting only confidently typed ACh/GABA/Glu neurons, R_eff falls to 7.4, still above 1.

**Post-reset recovery was predictable from wiring.** For non-runaway reset runs (n = 562), the best single index was the size of the excitatory bootstrap core. This is the maximal set of surviving neurons each of which receives enough excitatory synapses from other members to stay active at 100 Hz. It predicted FII_norm with cross-validated R² 0.964 [0.940, 0.978], an improvement over dead fraction alone of ΔR² +0.066 [+0.025, +0.122]. It separated recovered from non-recovered reset runs with cross-validated AUC 0.999 (n = 582).

Transfer to held-out experimental arms was good for tolerance (arm D) and hub protection (arm H) but poor for ageing (single-index R² −0.50; combined model 0.59). Wiring metrics trained on acute death therefore do not fully capture synapse-loss ageing.

**Runaway was not predictable from wiring.** Among viable networks after naive reperfusion (n = 482, 457 runaway), no index beat chance: best single-index AUC 0.53; combined random forest AUC 0.61 [0.38, 0.76]. Every brain-death run had R_eff = 0.

The structural index thus predicts how much of the original circuit can return, but not whether the network will be trapped. Whether it is trapped depends on protocol, temperature and (Section 2.6) sign assignment.

These results have not been written up as a standalone report. The planned BANC/MCNS/MANC recoverability maps were started but not completed.

### 2.11 Summary of preregistered hypotheses

| Hypothesis | Preregistered criterion | Result | Verdict |
|---|---|---|---|
| H1 sharp damage threshold | ΔAIC ≤ −10 and width ≤ 0.10 | Structural random loss: FAFB w 0.671 (4/4 datasets w > 0.25); dynamical: ΔAIC +0.11, width 1.56 | **Falsified** (hub-targeted structural loss: graded, not robust to edge cut-off) |
| H2 targeted vulnerability | (a) ρ(demand, death time) ≤ −0.3; (b) f_c ratio ≤ 0.7, CI excludes 1 | (b) area 0.346 [0.337, 0.356], hub 0.041 in FAFB; met in 4/4 connectomes. (a) not computed | **Partial: (b) supported; (a) untested** |
| H3 graded > abrupt reperfusion | ΔFII ≥ 0.05, negative control null | ΔFII −0.221 to −0.331 at 360 s; persists with excitotoxic term (k_exc ≤ 3) | **Falsified (reversed)** |
| H4 hypothermia ≈ Q10 | Q10_window in [0.67, 1.5] × model Q10 | 2.09 [1.99, 2.19] at Q10 2.5 (≈ Q10^0.8); MCNS direction replicates (+0.185) | **Supported** (set by the damage rule) |
| H5 identity lost first | D50_FII lowest; regime ≥5% of grid | Reset: D50s coincide, 0.0% coverage; runaway only: 22.1% | **Not supported** for recovery; present only in runaway |
| H6 replication | ≥3/4 connectomes match | Structural direction 4/4, magnitudes differ (f_c max/min 6.32–23.35); dynamics: MCNS replicates runaway, rescue, hub protection, temperature; BANC not analysed, MANC not simulated | **Partial** |

---

## 3. Discussion

### 3.1 What the model says

In a model built from a complete fly brain wiring diagram, the obstacle to recovery after anoxia was a state, not damage. Reoxygenation could leave the network in a self-sustained, K⁺-coupled high-activity attractor that coexists with the normal quiet state.

Ion handling maintains the attractor: it survives removal of every synapse and vanishes without activity-dependent K⁺ release. Entering it is a separate question. After 180 s of anoxia, entry was robust to parameters (48/48 sets), connectome (FAFB and MCNS) and transmitter-sign assignment. After 30 s of anoxia, entry depended on synaptic recruitment through antennal-lobe local neurons (147 in FAFB, 219 in MCNS) whose sign the default rule gets wrong. With these neurons inhibitory, the high-rate state did not form and identity was preserved.

The headline is therefore duration-dependent. The trap is robust after long insults, but its short-insult form is sensitive to sign predictions and should not be generalised.

**The reset.** A 10–20 s total block let K⁺ and ATP relax back into the quiet basin and restored the original stimulus–response map. This held at both durations, in both connectomes and under both sign assignments. After the reset, the brain that returned was the original one, minus the neurons that had died.

The line between "recoverable" and "dead" was set by two separate processes:
- a reversible dynamical trap determined whether function returned at all;
- slow injury kinetics determined how much of the original circuit was left to return.

**Failed expectations.** Three preregistered expectations failed, and each failure is informative.
- *H1.* No sharp percolation threshold appeared under random or dynamical damage. Identity fell roughly in proportion to dead fraction, and which neurons died mattered more than how many: hub-targeted loss collapsed reach at a small fraction of the random-loss damage, and protecting 1% of hubs was 30–84× more efficient than whole-brain protection.
- *H3.* Graded reperfusion was worse, not better, because it lengthened the depolarized period and did not remove the attractor.
- *H5.* "Resuscitated but different" networks appeared only as the pathological runaway state. They were not a stable alternative identity.

### 3.2 Relation to mammalian physiology

The ingredients of the trap are familiar. Extracellular K⁺ accumulation driving further depolarization is the core positive feedback of spreading depolarization [@dreier2011; @hartings2016]. Mammalian cortex after cardiac arrest often passes through burst-suppression, periodic discharges and seizures, which are associated with poor outcome. Our model does not reproduce these patterns, and we do not claim it does. The episodic bursting seen in MCNS after long anoxia resembles burst patterns only superficially.

The model does suggest a specific question for post-arrest electrophysiology: do some post-anoxic hyperexcitable states persist because the ionic environment sits in a second stable state, rather than because of irreversible injury? If so, a time-limited, deep suppression of activity, ended only after energy has recovered, would be a mechanistically targeted intervention. The model also predicts that the timing of release matters:
- a block that ended before ATP recovered re-ignited the attractor in 11 of 48 parameter sets;
- a release triggered as soon as K⁺ and ATP normalised (11 s) re-ignited it in MCNS and the FAFB-proxy.

**Temperature.** Temperature results agree with clinical and fly data in direction. Cooling before and during the insult was the strongest protective manipulation, consistent with fly K⁺ stability at 17 °C during repetitive anoxia [@rodriguez2012] and with clinical hypothermia after cardiac arrest [@haca2002; @bernard2002; @erecinska2003]. The size of the effect (Q10^0.8) follows from our injury rule and should not be read as a prediction. Cooling applied only during reperfusion was harmful, because slower pump-mediated K⁺ clearance deepened the trap. This prediction can be tested.

**Graded reperfusion.** The rodent postconditioning literature reports benefit from graded reperfusion [@gao2008; @zhao2009], which our model does not reproduce. The model lacks the oxidative, inflammatory and gene-expression mechanisms thought to mediate that benefit [@eltzschig2011; @caraballo2014; @vigne2009]. Our negative H3 result is therefore specific to the electrical and ionic component of reperfusion.

### 3.3 Relation to post-mortem restoration

BrainEx restored cellular function without global electrocorticographic activity [@vrselja2019]. That experiment used an activity-suppressing perfusate by design, so it does not test our model. Our results suggest that restoring cells and restoring organised network function are separate problems, with a reversible network state between them.

In our model, a silent unresponsive state never arose while most neurons were alive. In vivo, such a state would require mechanisms the model lacks, such as persistent synaptic failure, glial dysfunction or neuromodulatory collapse. Whether original pathways persist after resuscitation depends here on neuron survival, and above all on survival of hubs.

### 3.4 Brain death and resuscitation: what this model can and cannot say

**What it can say, as model predictions.**
- In a network built from a complete fly wiring diagram, failure to recover after anoxia can be a reversible ionic state rather than structural damage.
- Clearing that state with a time-limited total activity block returns the original stimulus–response map, minus dead neurons.
- The block has to outlast energy recovery.
- Silent, unresponsive and unrecoverable network failure required ≥97.9% neuron death.
- Which neurons survive matters more than how many.

**What it cannot say.**
- *Living tissue.* It cannot say whether any of this happens in living flies or mammals. The model has no glia beyond a clearance rate, no Ca²⁺ or Na⁺ dynamics, no neuromodulation, no oxidative or inflammatory injury, no plasticity and no gene expression.
- *Clinical questions.* It cannot address clinical death by neurologic criteria, organ donation practice or the care of any patient.
- *Persistence of pathways.* It cannot say whether "original pathways persist" in a biological sense. FII measures a seven-probe stimulus → descending/motor map and is blind to lesions outside those pathways (losing a third of the brain in the right optic lobe left FII ≈ 1.0).
- *Short-insult trap.* It cannot say whether the short-insult trap exists at all, since it hinges on transmitter-sign predictions for a small population.
- *Absolute timing.* It does not give absolute times: model seconds are compressed 50-fold, and injury constants were chosen for simulability, not fitted to flies.

Resuscitation of a brain-dead organism has not been demonstrated empirically, and this work does not demonstrate it either. It identifies, in a model, one class of network-level obstacle and one class of intervention that could be tested experimentally.

### 3.5 Testable predictions in *Drosophila*

1. **Reset rescue.** Consider flies exposed to N₂ anoxia beyond the duration that normally gives delayed recovery. They should recover faster, with fewer deficits, if neuronal activity is silenced pan-neuronally for a defined window at reoxygenation. Possible methods are GtACR1 optogenetic silencing, or shibire^ts at restrictive temperature applied with care for its thermal confound. Readouts: time to righting, proboscis extension to sucrose, and giant-fibre escape latency.
2. **Duration dependence.** Under the corrected-sign model, a reset should matter more after long anoxia than after short anoxia. Silencing should therefore show a larger benefit at longer insults.
3. **K⁺ trap signature.** K⁺-selective microelectrodes in the fly brain [@armstrong2011; @rodriguez2012] should show a post-reoxygenation K⁺ plateau of ~10–20 mM with hyperactivity in poorly recovering flies. Silencing should abolish it.
4. **Timing window.** Silencing started 5–20 s after reoxygenation should still rescue. Silencing that ends before K⁺ clearance and energy recovery should fail.
5. **Cooling timing.** Pre-insult cooling should protect, and cooling only during reoxygenation should not.
6. **Hub protection.** Protecting large, high-degree neurons, for example by targeted expression of protective factors such as glial Hsp70 [@armstrong2011], should preserve more sensorimotor function per protected cell than random protection.

### 3.6 Translational caution

The fly brain is anoxia-tolerant and treats SD as a protective shutdown [@robertson2020; @campbell2018]. Our injury parameters are mammal-like in timescale, chosen so that damage occurs within simulable durations. The extrapolated median-neuron LT50 under arrest was ~1.5 h, while real adult flies have an LT50 of 8.3 h [@callier2015]. Because of time compression (Methods), model seconds are not biological seconds. Nothing here bears directly on human brain death determination, organ donation practice or the treatment of individual patients.

---

## 4. Limitations

1. **Unfitted constants.** τ_ATP, J50, the injury hazard, g_K coupling, and the K⁺ release and clearance rates are phenomenological and not fitted to *Drosophila* recordings. SD latency and K⁺ peaks were matched; ATP kinetics (fly whole-body ATP 35% at 15 min; [@campbell2018]) were not. The Latin-hypercube analysis shows that the attractor exists across a 10-parameter box, while reset adequacy depends on τ_ATP and reserve.
2. **LIF abstraction.** The model has point neurons and no voltage-gated K⁺ or Na⁺ conductances beyond threshold and reset. It has no gap junctions, dendrites or neuromodulation. Real K⁺-dependent conductances could abolish or create bistability.
3. **No glia, Ca²⁺ or Na⁺ dynamics.** Glial K⁺ buffering [@andersen2026; @armstrong2011] is represented only as a clearance rate. Ca²⁺ load is a proxy integral, and excitotoxicity is a K⁺-based proxy.
4. **Time compression.** The slow layer runs 50× faster than spiking time (speedup 50). Outcome classes were invariant from 25× to 100×, but absolute escape times and therapeutic windows in biological units are uncertain.
5. **Sign assignment.** Transmitter signs come from predictions, with glutamate inhibitory and unknown or monoaminergic transmitters excitatory. Glutamate sign was the strongest determinant of runaway intensity (PRCC −0.71). Correcting antennal-lobe local neurons abolished the high-rate 30-s runaway (Section 2.6). This correction was tested with 3 seeds in the full FAFB model and 5 seeds in MCNS and the FAFB-proxy, at 30 and 180 s only. Broader sign variants (all low-confidence calls inhibitory; no-call neurons silent) were not run. Most FAFB results (Sections 2.3–2.8 and 2.10) use the default signs.
6. **Two-dataset dynamics.** Dynamic results come from FAFB and MCNS, and MCNS required proxies for demand and K⁺ coupling. A BANC sweep finished but was not analysed. MANC was not simulated dynamically. Sex and dataset are confounded.
7. **Readout dependence.** FII is defined by a seven-probe battery read at descending and motor neurons. Large lesions outside the battery's pathways leave FII intact. A single 500-ms evaluation window can fall in a silent gap of an episodic runaway (MCNS 180 s).
8. **Statistical design below preregistration.** We used 3–5 seeds per condition (preregistered ≥10) and one stimulus rate per probe (preregistered 3). Several structural regimes had fewer than 20 removal orders. Several between-seed CIs are degenerate because outcomes are nearly deterministic.
9. **No spontaneous activity.** The Shiu model is silent at rest, so the preregistered spontaneously-active-fraction metric was replaced by an evoked-neuron count.
10. **Incomplete follow-ups.** The following were not executed because the compute budget was reached:
    - a 605-run mechanism follow-up: longer blocks for all parameter sets, the full high-gain excitotoxic grid, finer targeted-set brackets and a wider hysteresis ramp;
    - the BANC analysis;
    - the cross-connectome recoverability maps;
    - H2(a), the correlation of death time with demand, which was never computed.

---

## 5. Methods

### 5.1 Data

The following connectomes were downloaded from FlyWire Codex (https://codex.flywire.ai/) [@dorkenwald2024; @schlegel2024; @bates2026; @takemura2024; @marin2024; @berg2026]:
- FAFB v783 (FlyWire; 139,255 neurons; Princeton synapse table aggregated over neuropils, ≥5 synapses per pair; 3,732,460 pairs and 50,666,648 synapses);
- BANC v888 (158,262 neurons; 3,037,361 pairs, release floored at ≥3 synapses);
- MANC v1.2.1 (23,665 neurons; 5,305,638 pairs);
- MCNS v1.0 (166,700 neurons; 6,242,118 pairs at ≥5).

Pair counts matched the Codex totals. MAOL v1.1 [@nern2025] was not available locally and was not analysed. Neuron annotations (super-class, class, cell type, side, predicted transmitter), FAFB cell size statistics and FlyWire skeletons were used as distributed with each release. Glia, trachea and non-neuronal objects were excluded (BANC 473 objects; MANC 15 glia).

### 5.2 Fast layer

The Shiu et al. LIF model [@shiu2024] was reimplemented in NumPy/SciPy (`flysim.py`; brian2 and numba were not installable on the host).
- **Dynamics:** dv/dt = (v₀ − v + g)/τ_m and dg/dt = −g/τ_syn, with exact integration at dt = 0.1 ms.
- **Synapses:** a presynaptic spike adds w_syn·sign·n_syn to g after a 1.8-ms delay.
- **Refractory period:** 2.2 ms, during which v and g are frozen.
- **Input:** Poisson input adds 68.75 mV to v.
- **Integration:** the Brian2 schedule order was preserved. Integration is event-driven over an active set, which gave identical results to full integration.
- **Weights:** FAFB used w_syn = 0.185 mV (Section 2.1). For MCNS and BANC, w_syn was matched by scaling each dataset's reverberation threshold by the FAFB ratio 0.185/0.2135 (MCNS 0.0538, BANC 0.2227).

### 5.3 Slow metabolic and ionic layer

The slow layer was updated every 5 ms of spiking time, with all slow rates multiplied by a speedup factor (50 unless stated). Times are reported in model seconds.

- **ATP** a_i: da/dt = [m·b_i·Q_s·S·(1−a)/(K_adp+1−a) − d_i·a]/τ_ATP, with reserve m = 4 × resting demand and a₀ = 0.95.
  - Demand: d_i = b_i·Q_d·(1 + r_i/r_half + k_K·(K_i−K₀)/K₀), with b_i = clip((area_i/median)^β, 0.25, 4), β = 0.5 and r_half = 100 Hz.
  - τ_ATP = 30 s (calibrated).
  - Temperature: Q_d = 2.5^((T−22)/10); supply Q10 1.5; pump/clearance Q10 1.8; injury Q10 2.5.
- **Pump failure**: dV_atp = 40·((e_ok − e)/e_ok)^1.5 for e < e_ok = 0.75 (zero at full supply).
- **Extracellular K⁺** in 79 neuropil compartments (K₀ = 5 mM), with four terms:
  - spike-driven release (k_spk);
  - regenerative release dep·R·(1 − K/K_max), with K_max = 80 mM and a depleting release store (τ 8 s; recovery τ 30 s);
  - pump-dependent leak and clearance (k_clr = 1 s⁻¹ × pump, half glial);
  - diffusion over a geometric adjacency (neuropil contact area / centroid distance, from 80,215,790 synapse positions).

  Neuron K_i is the synapse-weighted mixture of its compartments. dV_K = g_K·(RT/F)·ln(K/K₀), with g_K = 0.6, and total depolarization is soft-capped at 50 mV.
- **SD and depolarization block**: dep = σ((v_rest,eff − V_sd)/2), with V_sd = −30 ± 2 mV. Block h → 1 − σ((v_rest,eff + 32)/2), raising threshold by 100·(1 − h) mV.
- **Injury and death**: dJ/dt = dep·(1 − e) − k_rep·e·(1 − dep)·J. Death hazard h_max·J⁴/(J⁴ + J50⁴), with h_max = 1/30 s⁻¹ and J50 = 600 s (calibrated; 60 s killed 90.6% after 120 s and was rejected). Death is permanent. Under this hazard the cumulative hazard scales as rate⁴·D⁵, so D_max ∝ rate^(−4/5), which explains Q10_window ≈ Q10^0.8.
- **Synaptic failure** (reversible): release fraction clip((e − 0.3)/0.3, 0, 1).
- **Excitotoxic term** (follow-up only): dJ_i/dt += q·k_exc·max(r_i − r₀, 0)/r₀·Kex_i/(Kex_i + K₀), with r₀ = 50 Hz and k_exc ∈ {0, 0.3, 1, 3, 10}.

Parameter provenance and literature anchors are listed in the Supplement (Table S1, from the project parameter table of 81 values, 6 flagged unverified).

### 5.4 Protocols and readouts

**Run timeline:** 2 s normoxia; anoxia (supply S = 0) for D s; 60 s of reperfusion protocol; 20 s free observation (S = 1, no interventions); evaluation at speedup 1.

**Probe battery:**
- sugar GRNs, bitter GRNs, and Johnston's-organ auditory and wind/gravity neurons (150 Hz);
- bilateral LC4 and LC9 (100 Hz);
- left turn (LC9_L + DNa02_L, 180 Hz).

Each probe ran for 250 ms from a copy of the post-protocol state, with readout at all 1,415 descending and motor neurons. Response = probe rate − rate in a 500-ms unstimulated window.

**Interventions:**
- `block_all` (total spike and synaptic block);
- hyperpolarizing offsets, GABA gain, SD-release scaling and clearance scaling;
- ramps and pulses of supply, and temperature schedules;
- noise and adaptation;
- per-neuron injury protection and counterfactual revival;
- in MCNS and the FAFB-proxy, a closed-loop block released when the maximum compartment K⁺ < 6 mM and mean ATP e > 0.9.

The protocol list is in Table S5b.

### 5.5 Metrics and outcome classes

- **FII** = Pearson r of vec(R_condition) vs vec(R_baseline) over entries > 0 in either, with baseline = mean of 10 intact seeds.
- **Noise ceiling** = median leave-one-out single-seed r (FAFB 0.9940; MCNS 0.9972; FAFB-proxy 0.9922); FII_norm = FII/ceiling.
- **Reach** = fraction of (probe, readout) pairs joined by ≤6 hops through surviving neurons.
- **Outcome classes:**
  - runaway: late rate > 0.5 Hz per neuron or spontaneous fraction > 0.05;
  - brain death: spontaneous fraction < 0.01 and no readout response > 8 Hz;
  - recovered: FII_norm ≥ 0.8 and Reach ≥ 0.8 × baseline;
  - structural loss: ≥10% dead, not recovered;
  - impaired: otherwise.
- **D_max** = duration at which mean FII_norm crosses 0.8 (log-interpolated).

### 5.6 Structural percolation

**Reach computation.** Strict reach used 128 sensory source neurons stratified by modality, with sinks = descending + motor (+ efferent) neurons and relays only through removable interneurons. It was computed with uint64 bitset propagation.

**Removal regimes:**
- random;
- demand-weighted (Efraimidis–Spirakis order, weights ∝ membrane area or synapse count);
- hub-targeted by total degree;
- degree-bin-matched random (40 log-degree bins);
- focal neuropil lesions.

**Fits and statistics.** Fits compared 4-parameter logistic, linear and exponential models by Gaussian AIC. Ratios used bootstrap over removal orders (2,000 resamples for the primary metric).

**Controls:**
- degree-preserving directed rewiring (FAFB, MANC);
- edge cut-offs ≥10 (all datasets) and ≥1 (MANC);
- removal of off-path neurons (Δ = 0 in all datasets).

Keystone analysis evaluated top-flow-centrality cell types (Supplement S4).

### 5.7 Mechanism experiments

The mechanism track comprised 2,061 runs in total:
- Latin-hypercube robustness (48 sets × 10 parameters, 784 runs);
- the excitotoxic grid (12 protocols × 3 durations × 5 seeds at k_exc 1 and 3);
- ignition anatomy with 5-ms first-spike resolution;
- driver-set ranking on seeds 0–4, with targeted blocks tested on held-out seeds 10–14;
- dissection (synaptic removal, K⁺ clamp, k_spk = 0, k_rel = 0);
- bistability (2-s K⁺ clamp kicks, tonic offsets δ, slow δ ramps in both directions, and block quench).

### 5.8 Transmitter-sign variant

Variant A set the sign of every neuron with class ALLN or primary type matching lLN, vLN, v2LN, ilLN or il2LN to −1. In FAFB this covers 431 neurons, of which 147 changed sign. In the FAFB full model, variant A was run for 30-s anoxia with naive reperfusion: 3 seeds of the probe protocol plus 3 anatomy runs with first-spike recording. These runs were made locally with independent processes. In MCNS and the FAFB-proxy, the equivalent variant ("AL-LN −1") was run with its own 5-seed intact baseline for 30 and 180 s, naive and reset (25 runs each).

### 5.9 Cross-connectome dynamic replication

**MCNS and BANC models.** These were built with `flysim_cns.py`, using:
- the same LIF and slow-layer code;
- demand proxied by total synapse count;
- a connection-based Jaccard 9-nearest-neighbour neuropil graph for K⁺ coupling (a dense shared-neuron adjacency was numerically unstable and was replaced);
- w_syn matched to the FAFB distance below the reverberation threshold.

**Probes and readouts.** Probes were matched by cell type. A probe counted as responsive if ≥3 brain readouts exceeded 8 Hz in the baseline mean. VNC identity used VNC motor and efferent readouts with the LC4 probe.

**Runs.** The MCNS sweep (106 main runs, 25 sign-variant runs, 15 low-gain runs) and the FAFB-proxy control (40 + 25 runs) ran on cloud CPU containers. Core contrasts used 5 seeds, the MCNS baseline 10, and sham and sensitivity arms 3.

### 5.10 Recoverability index

Nineteen graph indices were computed on the surviving FAFB graph of each processed dynamics run:
- mean-field reproduction numbers (all neurons and confident-transmitter neurons), and a signed leading eigenvalue;
- excitatory/inhibitory synapse balance;
- hub survival;
- rich-club density;
- reciprocal and 3-cycle survival;
- giant excitatory SCC;
- excitatory bootstrap-core size;
- simulator Reach.

Prediction used condition-grouped 5-fold cross-validation (4–20 repeats depending on the model), with cluster-bootstrap CIs, as well as logistic and random-forest models and leave-arm-out transfer.

### 5.11 Statistics

Bootstrap 95% CIs (2,000–4,000 resamples) were computed over seeds, removal orders or parameter sets. Protocol and variant contrasts used a paired-by-seed bootstrap. PRCC used Holm correction over 250 tests. The six primary tests were judged on preregistered effect-size criteria, and Holm correction across them changes no verdict.

### 5.12 Preregistration and deviations

Hypotheses, definitions and thresholds were registered on 2026-10-16, before any insult simulation (Supplement S1). Deviations are listed in Supplement S9. They include:
- seed counts and single-rate probes;
- replacement of the spontaneous-activity metric;
- J50/τ_ATP calibration;
- unmatched integrated supply in H3 protocols;
- reduced removal orders;
- the unexecuted follow-ups.

### 5.13 Compute

Dynamics, mechanism and replication sweeps ran on Modal CPU containers:
- dynamics: 4 × 32 CPU, 24.3 CPU-h;
- mechanism: 2 × 32 CPU, 34–37 min wall;
- replication: 32-CPU containers, 217 s and 269 s wall for the analysed jobs.

Structural analyses, the FAFB sign variant and the recoverability index ran locally on a 16-core Windows host.

---

## Data availability

Connectome data (FAFB v783, BANC v888, MANC v1.2.1, MCNS v1.0) are available from FlyWire Codex (https://codex.flywire.ai/) under the terms described at https://codex.flywire.ai/about_codex and https://flywire.ai/principles.html; they are not redistributed. All derived result tables (per-condition summaries, per-run metrics, contrasts) and figure source data will be deposited with the archived code release (Zenodo; DOI to be minted at publication). Raw per-run outputs (neuron-level spike counts and alive masks) will be deposited where size permits. The 146 BANC replication result files remain on the compute volume and were not analysed for this version.

## Code availability

All code (the simulator `flysim.py` and its variants, sweep drivers, and analysis and figure scripts) will be released at https://github.com/KalarisLabs/fly-brain-resuscitation (planned; not yet public at the time of this version) under an open-source licence, with an archived Zenodo release (DOI to be minted). The build scripts that generate this manuscript and its supplementary tables directly from the result files are included.

## AI-use disclosure

This study was carried out with substantial assistance from an AI system, Claude (Anthropic), under the direction of the author. The agentic tooling came from Kalaris Labs: Research Agent Skills (https://github.com/KalarisLabs/research-agent-skills) and Myriad (https://github.com/KalarisLabs/myriad). Aryan Kumar Jha (Kalaris Labs) is credited for AI (Claude) development work.

The AI system:
- retrieved and verified literature (all 82 bibliography entries were checked against Crossref records);
- wrote and ran simulation and analysis code;
- drafted text and figures.

The author set the questions and scope, directed the research and is responsible for all content. The AI system is not an author. All numerical results were copied programmatically or verbatim from analysis outputs, not generated by the language model.

## Author contributions (CRediT)

**Srihari Muralikrishnan** (sole author; Kalaris Labs; ORCID 0009-0004-0220-0361): Conceptualization, Methodology, Software (AI-assisted), Formal analysis, Investigation, Resources (primary compute), Data curation, Visualization, Writing – original draft (AI-assisted), Writing – review & editing, Supervision, Project administration.

**Non-author contributor.** Sayan Chowdhury (Founder & CEO, Kalaris Labs): Resources (secondary compute provision), institutional support.

## Competing interests

The author is a co-founder of Kalaris Labs. Sayan Chowdhury, who provided secondary compute and institutional support, is Founder & CEO of Kalaris Labs. No other competing interests are declared.

## Funding

This work received no external funding. It was self-funded by the author, who paid the primary compute costs. Secondary compute was provided in kind by Sayan Chowdhury (Kalaris Labs).

## Ethics statement

This is a purely computational study. No animals, human participants or human data were used. All simulations used publicly released connectome reconstructions.

## Acknowledgments

We thank Sayan Chowdhury (Founder & CEO, Kalaris Labs) for secondary compute provision and institutional support, and Aryan Kumar Jha (Kalaris Labs) for AI development work.

We acknowledge the FlyWire Consortium, the Princeton FlyWire team, and the teams that produced the BANC, MANC and male CNS connectomes. They made their data openly available through FlyWire Codex (https://codex.flywire.ai/about_codex; https://codex.flywire.ai/about_flywire). Use of these data follows the FlyWire principles (https://flywire.ai/principles.html), and we cite the dataset publications [@dorkenwald2024; @schlegel2024; @bates2026; @takemura2024; @marin2024; @berg2026].

The LIF model builds on Shiu et al. [@shiu2024]. The Kalaris Labs Myriad and Research Agent Skills were used throughout.

---

## References

{{REFERENCES}}
