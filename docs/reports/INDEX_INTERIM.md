# Recoverability index — INTERIM findings (index track, verbatim from sub-agent report; not cross-validated)

Source: index sub-agent 07de2ee0 report, 2026-10-18. Data: D:\FRUITFLY\research\followup\index\data\run_indices.csv (~695 of 1,209 runs processed).

**Findings so far** (correlations only, no cross-validation yet; dynamic results are from the ~695 processed runs)
- **Intact FAFB is far above criticality:** using the model's own neuron parameters, the excitatory network's reproduction number (extra spikes per spike) is R_eff ≈ 19.6, or 17.0 with inhibition included.
- **That excitatory capacity is partly a sign artefact:** the strongest excitatory loop sits on about 99 antennal-lobe local neurons. These are inhibitory in real flies, but the model's sign rule counts them as excitatory. Counting only neurons with a confident ACh/GABA/Glu call, R_eff drops to 7.4, still above 1.
- **Identity after the reset protocol:** three indices tie as the best predictors of FII_norm (rank correlation 0.86–0.88): the simulator's sensory-to-motor Reach, survival of the top-1 % hub neurons, and the confident-transmitter R_eff. Dead fraction alone gives −0.83.
- **Runaway after naive reperfusion:** no connectome index predicts it. 209 of 245 runaway runs have the intact structure (R_eff 19.6), the same as runs that recovered. Temperature and protocol decide the outcome. R_eff > 1 is necessary for runaway but does not predict it.
- **Brain death:** every brain-death run so far has R_eff = 0. The only runaway-class runs with R_eff below 1 are the 98.7 %-dead ones, and they fire weakly (about 5.5 Hz).

**Deviations to declare**
- About 514 runs were not processed.
- No cross-validated models or confidence intervals yet.
- Steps 3 and 4 are not done: the BANC/MCNS/MANC recoverability map, and the other-species proposal (outline only).
- Ageing synapse loss is mostly not applied.
- The signed "self-sustaining core" index is capped (40 iterations, 200 Hz only) because it never converges.
- Modal was used for mask extraction only.
- `INDEX_RESULTS.md` and the figures are not produced.