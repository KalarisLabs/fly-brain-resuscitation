# Results summary

> Model predictions only. Every number below is copied verbatim from [`SYNTHESIS.md`](SYNTHESIS.md) or [`results/KEY_CLAIMS.json`](../results/KEY_CLAIMS.json), the single source of truth for this project.

Claim set version 1.0 (2026-10-20). Usage rule: Downstream messaging must use only these claims, with their numbers and confidence qualifiers; always include K15. Paths are relative to D:\FRUITFLY.

## Key claims

| id | claim | number | confidence |
|---|---|---|---|
| K01 | Headline: in a connectome-constrained LIF model with a phenomenological K+/ATP layer, abrupt reoxygenation can trap the model brain in a self-sustained K+-coupled high-activity state; a 10-20 s total activity block releases it. The trap is robust after long (180 s) anoxia and sign-sensitive after short (30 s) anoxia. | 180 s: runaway in 48/48 parameter sets (FAFB) and 5/5 MCNS runs; 30 s with antennal-lobe local neurons inhibitory: 0/3 FAFB runaway | high (180 s); moderate (30 s sign sensitivity) |
| K02 | Under the default transmitter-sign rule, naive (abrupt) reoxygenation produced the runaway after as little as 7 s of anoxia at 22 C in the FAFB model. | 30/30 runs from 7 to 480 s anoxia; FII_norm ~0; late firing 105-119 Hz/neuron | high for default signs; NOT valid under corrected antennal-lobe signs (see K03) - do not use '7 s' without this qualifier |
| K03 | Setting antennal-lobe local neurons inhibitory (147 sign changes in FAFB, 219 in MCNS) removed the high-rate runaway and the identity loss after 30 s anoxia. MCNS and FAFB-proxy kept a weak low-rate oscillation that the preregistered classifier still counts as runaway. | FAFB variant A 0/3 runaway (late 0.035-0.17 Hz/neuron, raw FII 0.915-0.936); MCNS FII_norm 1.000 [0.996, 1.002] (3.0 Hz/neuron residual), delta +0.949 [+0.918, +0.977]; FAFB-proxy 0.859 (0.82 Hz/neuron) | moderate (n = 3 FAFB, n = 5 MCNS/proxy; variants B, C not run) |
| K04 | After 180 s anoxia the runaway persisted with corrected antennal-lobe signs, and the 20-s total block still rescued identity. | MCNS naive FII_norm -0.135 (FSA 0.897) vs block 0.989; FAFB-proxy naive -0.125 vs block 0.965 | moderate-high |
| K05 | A 20-s total activity block at reperfusion restored the stimulus-response map; 2-5 s blocks failed; the block worked even when started 20 s late. | FAFB FII_norm 1.001 / 0.953 at 30 / 180 s anoxia; rescue in 48/48 parameter sets at 30 s, 34/48 (0.71 [0.58, 0.83]) at 180 s | high |
| K06 | Releasing the block as soon as K+ and ATP had normalised was not enough: the network re-ignited. | closed-loop criterion met after 11 s; MCNS runaway 5/5 (11.3 Hz/neuron), delta vs 20-s block -0.071 [-0.081, -0.058]; FAFB-proxy delta -0.674 [-0.742, -0.624] | moderate (MCNS + proxy only, 180 s, default signs) |
| K07 | The maintained runaway is ionic, not synaptic: it survives removal of all synapses and collapses when activity-dependent K+ accumulation is removed; quiet and active states coexist in the intact network. | all synapses off: 109/134 Hz/neuron; K+ clamp 5 mM: 2.4/2.3 Hz/neuron; hysteresis >= 13 mV | high within model (default signs) |
| K08 | Simulated 'brain death' (silent and unresponsive) occurred only with near-total neuron death; the model never produced a silent but structurally intact brain. | 92 runs, all with >= 97.9% of neurons dead (of 1,209 FAFB runs) | high (model property; FAFB only) |
| K09 | After the reset, remaining identity loss was due to neuron death; protecting the top 1% most-connected neurons was far more efficient per neuron than whole-brain protection, in two connectomes. | FAFB ~30x (0.578 -> 0.703); MCNS 84x (+0.201 [+0.127, +0.278]); random sets no gain | high |
| K10 | Cooling before/during the insult widened the recoverable window, but the size of the effect is set by the model's injury rule; cooling only during reperfusion was harmful. | D_max 359 / 271 / 149 / 124 model-s at 17 / 22 / 29 / 33 C; Q10_window 2.09 [1.99, 2.19] | high that it holds; magnitude is by construction, not a prediction |
| K11 | Graded (ramped or pulsed) reperfusion did not beat abrupt reperfusion plus a block (preregistered H3 falsified, direction reversed). | 360 s: ramp30 delta FII_norm -0.331 [-0.583, -0.130]; ramp60 +6.6 pp dead | high within model (no oxidative/inflammatory mechanisms modelled) |
| K12 | Structurally, random neuron loss degrades sensorimotor reachability gradually (no sharp threshold), whereas hub- or demand-targeted loss is far more destructive, in all four connectomes. | FAFB random f_c 0.800 (width 0.671); area-weighted 0.277 (ratio 0.346); hub 0.033 (ratio 0.041); direction replicates 4/4 | high (direction); magnitudes differ across datasets (f_c max/min up to 23.35) |
| K13 | Wiring statistics predict how much identity returns after the reset, but not whether the network gets trapped. | post-reset FII: cross-validated R2 0.964 [0.940, 0.978] (excitatory core size); runaway among viable networks: AUC 0.53 (chance) | moderate (FAFB only, 1,118 of 1,209 runs, default signs, no standalone report) |
| K14 | The small 'driver core' (0.5% of neurons) that rescued after 30 s anoxia under default signs is not established: it is centred on the antennal-lobe neurons whose sign correction removes the 30-s runaway. | 279 of 696 driver neurons in the antennal lobes; after 180 s, 20% of the brain had to be blocked | low (likely sign artefact) - do not message as a finding |
| K15 | Scope: all results are predictions of a simplified simulation (connectome-constrained LIF + phenomenological metabolism, constants not fitted to Drosophila). They are not evidence about living tissue, clinical brain death or human resuscitation; resuscitation of a brain-dead organism has not been demonstrated. | 3,487 analysed simulated condition-seed runs; dynamics in 2 connectomes (FAFB, MCNS); BANC dynamics not analysed | n/a (mandatory caveat for all messaging) |

## Revised headline (from SYNTHESIS.md, section 0)

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

The full claim-by-claim verdicts (preregistered hypotheses H1-H5, claims C01-C14) are in [`SYNTHESIS.md`](SYNTHESIS.md); the track reports are in [`reports/`](reports/).
