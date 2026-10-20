# Real-fly validation protocol (proposed; not executed)

> Status: **design only.** No experiment in living *Drosophila* has been run for this project. Every model number quoted here is copied from [`results/KEY_CLAIMS.json`](../results/KEY_CLAIMS.json) (claim ids K01-K15) or [`SYNTHESIS.md`](SYNTHESIS.md). Model seconds are compressed 50x relative to biology and the injury constants are mammal-like (SYNTHESIS.md, section 4), so **no model duration or concentration should be used as an experimental set-point**. Durations, sample sizes and thresholds below must be set by pilot experiments and a power analysis.

The protocol operationalises the six testable predictions of the paper (section 3.5, "Testable predictions in *Drosophila*"). Its purpose is to falsify the model, not to demonstrate resuscitation. Resuscitation of a brain-dead organism has not been demonstrated empirically, and this work does not demonstrate it (K15).

## 1. Predictions under test

| # | prediction (paper section 3.5) | model basis | outcome that would falsify it |
|---|---|---|---|
| P1 | **Reset rescue.** After N2 anoxia longer than the duration that normally gives delayed recovery, flies recover faster and with fewer deficits if neuronal activity is silenced pan-neuronally for a defined window at reoxygenation. | K01, K05 | Silenced flies do not recover faster than unsilenced controls with matched genotype and temperature history. |
| P2 | **Duration dependence.** Under the corrected-sign model, the benefit of a reset is larger after long anoxia than after short anoxia. | K03, K04 | The silencing benefit is equal or larger after short insults. |
| P3 | **K+ trap signature.** Poorly recovering flies show a sustained post-reoxygenation rise in extracellular K+ with hyperactivity, and silencing abolishes it. | K07 | No sustained K+ elevation in poor recoverers, or silencing leaves it unchanged. |
| P4 | **Timing window.** Silencing started late after reoxygenation still rescues (in the model, the 20-s block worked even when started 20 s late). Silencing released as soon as K+ and ATP have normalised fails (in the model the network re-ignited). | K05, K06 | Late-onset silencing fails, or early release is as effective as prolonged silencing. |
| P5 | **Cooling timing.** Pre-insult cooling protects; cooling only during reoxygenation does not. | K10 | Reoxygenation-only cooling protects as well as pre-insult cooling. |
| P6 | **Hub protection.** Protecting large, high-degree neurons (for example targeted expression of protective factors such as glial Hsp70) preserves more sensorimotor function per protected cell than random protection. | K09 | Hub-targeted and size-matched random protection give the same benefit per cell. |

## 2. Animals and insult

- Adult flies of one sex per experiment (the model results come from a female (FAFB) and a male (MCNS) connectome; sex and dataset are confounded, so test both sexes separately).
- Controlled-temperature rearing; record age, genotype, rearing temperature and time of day.
- Anoxia: pure N2 in a gas-exchange chamber with an O2 sensor. Establish, in a pilot, the duration-recovery curve (time to righting vs anoxia duration). Choose a "short" and a "long" duration from that curve: the long one must give delayed or incomplete recovery in untreated flies.
- Reoxygenation: abrupt switch to air (the condition in which the model predicts the trap).

## 3. Interventions

| arm | method | control |
|---|---|---|
| Pan-neuronal silencing at reoxygenation (P1, P2, P4) | GtACR1 (pan-neuronal GAL4 > UAS-GtACR1) with all-trans-retinal and light; window length and onset delay as factors | same genotype, no light; no-retinal light control |
| Thermogenetic silencing (P1, alternative) | shibire^ts at restrictive temperature | temperature-matched wild-type, to separate the thermal confound noted in the paper |
| Early release (P4) | silencing ended as soon as the K+ signal (P3 assay) returns to baseline | silencing held for the full window |
| Cooling (P5) | pre-insult cooling vs cooling only during reoxygenation | constant temperature |
| Hub protection (P6) | targeted expression of a protective factor in large/high-degree neuron classes (selected from connectome degree, `results/structure/*_type_flow_centrality.csv`) | expression in a size-matched random driver set |

## 4. Readouts

1. Behaviour: time to righting; proboscis extension to sucrose; giant-fibre escape latency (paper section 3.5). These map onto the model's stimulus -> descending/motor readout (FII).
2. Physiology: K+-selective microelectrode recordings in the brain during and after reoxygenation (paper refs. armstrong2011, rodriguez2012), with field-potential or calcium-imaging readout of hyperactivity.
3. Survival and gross deficits at 24 h.

## 5. Design and analysis

- Preregister hypotheses, arms, primary readout (time to righting) and stopping rules before data collection.
- Randomise flies to arms and blind the scorer to genotype and treatment.
- Set sample sizes by power analysis on pilot variance. (The simulation itself used 3-5 seeds per condition, below its own preregistered 10, SYNTHESIS.md section 3, so it gives no basis for effect-size assumptions.)
- Analyse with survival or time-to-event models (righting) and mixed models with day/batch as random effects; report effect sizes with 95% CIs.

## 6. Interpretation guard-rails

- A positive P1 result would show only that post-anoxic silencing helps flies; it would not show that the in-model ionic trap is the mechanism unless P3 and P4 also hold.
- The short-anoxia trap depends on transmitter-sign predictions for about 147 antennal-lobe local neurons (K03, K14). Short-insult results are therefore the most likely to diverge from the model.
- Flies are anoxia-tolerant and treat spreading depolarisation as protective (paper section 3.6). A null result in flies would not by itself refute the mechanism in other species, and a positive one would not transfer to mammals or humans.
- Nothing here bears on clinical brain-death determination, organ donation or patient care.
