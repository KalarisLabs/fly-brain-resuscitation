# Limitations

> Model predictions only. Every number below is copied verbatim from [`SYNTHESIS.md`](SYNTHESIS.md) or [`results/KEY_CLAIMS.json`](../results/KEY_CLAIMS.json), the single source of truth for this project.

## Model limitations (paper, section 4)

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

## Translational caution (paper, section 3.6)

The fly brain is anoxia-tolerant and treats SD as a protective shutdown [@robertson2020; @campbell2018]. Our injury parameters are mammal-like in timescale, chosen so that damage occurs within simulable durations. The extrapolated median-neuron LT50 under arrest was ~1.5 h, while real adult flies have an LT50 of 8.3 h [@callier2015]. Because of time compression (Methods), model seconds are not biological seconds. Nothing here bears directly on human brain death determination, organ donation practice or the treatment of individual patients.

## What was not done (SYNTHESIS.md, section 3)

1. **BANC dynamic replication:** 146 result files exist on `flysim-data:/results_rep_v4/BANC*` but were never harvested or analysed (job ef7d8554 recorded as failed; cost unknown, ≤ $1.76).
2. **Sign variants B and C** and the FAFB variant-A arms beyond 30-s naive were not run. The missing variant-A arms are the intact baseline, reset, 180 s, bistability kicks and targeted blocks. Variant-A FII is therefore relative to the default-sign intact baseline (FII_norm ≈ FII/0.9940 = 0.920–0.942).
3. **Mechanism follow-up (605 runs) not executed:** 40/60-s blocks for all LHS sets, the full k_exc = 10 H3 grid, finer targeted-set brackets, and a −30 mV hysteresis ramp.
4. **Recoverability index:**
   - 91 of 1,209 runs not processed;
   - the BANC/MCNS/MANC recoverability map was only started (`cross_all.log`: FAFB random 0–0.6 only);
   - INDEX_RESULTS.md and index figures were not written.
5. **H2(a)** (death-time vs demand) was never computed.
6. **Preregistered seed count** (≥10 per condition) was not met (3–5).

## What this model can and cannot say (SYNTHESIS.md, section 4)

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

Citation keys of the form `[@key]` resolve in [`paper/references.bib`](../paper/references.bib).
