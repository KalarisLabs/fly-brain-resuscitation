# Follow-up track 2 — cross-connectome DYNAMIC replication (MCNS v1.0; BANC v888 not analysed)

**Status (2026-10-19).** The v4 sweep ran on Modal (company workspace). Two jobs were **harvested and analysed**:
* smoke `14de3ae1-d0bf-412c-8841-c1a279754447`: MCNS 360-s hub arm, 45 runs, 217 s;
* job A `a9ef9894-0c76-4de2-89c9-6879f59167ff`: MCNS remaining main arms 61 + AL-sign-corrected 25 + w_syn sensitivity 15; FAFB-proxy main 40 + AL-corrected 25; 269 s.

That is **211 analysed runs**: MCNS 146 (106 main + 25 + 15) and FAFB-proxy 65.

**BANC job `ef7d8554-1bbc-4a61-a027-11bb646177dc` is recorded as failed with 0 files harvested.** A Volume listing taken before the pause showed all 146 BANC result files present on `flysim-data:/results_rep_v4/BANC*`. They were not retrieved: that needs a Modal sandbox, and new Modal work was stopped. **BANC is therefore not analysed.**

Seeds:
* Core contrasts (baseline excepted): 5 seeds (0–4, sequential).
* Baseline: 10 seeds for MCNS, 5 for the FAFB proxy and the AL-corrected arms.
* Sham and the w_syn-sensitivity arm: 3 seeds.

Bootstrap 95 % CIs use 2,000 resamples (contrasts: 4,000). An earlier local MCNS pilot (3 seeds; D30 naive/reset) gave the same answer (FII_norm 0.060 vs 1.001); it is superseded and kept in `raw/local_v3/`.
All statements are about this model (Shiu LIF + phenomenological metabolic layer, constants not fitted to *Drosophila*); none is evidence about real brains.

## 1. Models and calibration (unchanged from earlier report; summary)
* **MCNS:** 166,700 neurons, 6,242,085 edges (≥5 syn), 89.9 M synapses, 62.5 % excitatory synapses.
* **BANC:** 157,789 neurons, 3,036,600 edges (≥3 syn). Built, but its dynamics were not analysed.
* **Loader:** `flysim_cns.py`. flysim.py and flysim_dyn.py are untouched.
* **Signs:** GABA/Glu −1 (Shiu); HIST −1; unknown +1.
* **Substitutions:**
  - demand proxy = total synapse count (no size data);
  - K⁺ adjacency = connection-based Jaccard 9-NN (no geometry). It replaced the dense shared-neuron adjacency, which is numerically unstable (§6).
* **w_syn:** the same distance below the reverberation threshold as FAFB (0.185 = 0.867 × FAFB any-probe threshold 0.2135). This gives **MCNS 0.0538** (BANC 0.2227). Thresholds per probe are in `raw/calib/`.
* **FII:**
  - Preregistered matrix FII over brain readouts (descending neurons + brain motor neurons, 1,439) and **responsive probes only**. A probe counts as responsive if ≥3 brain readouts exceed 8 Hz in the baseline mean. MCNS responsive probes: bitter, JO_aud, JO_wg, LC4, LC9, turnL; sugar→MN9 is silent in MCNS and is excluded.
  - Normalised by the leave-one-out noise ceiling (MCNS 0.9972; FAFB-proxy 0.9922).
  - The VNC FII uses VNC motor/efferent readouts (802); only LC4 drives them at baseline.
* **Outcome classes** as in the FAFB dynamics track (runaway = late observation rate >0.5 Hz/neuron or FSA >0.05).

## 2. Results
{{artifact:f34e654a-e7b0-45ec-adda-3e7176836271}}

*Fig. 1 — FII_norm (responsive probes) per condition. Dots are runs (n = 5 seeds; FAFB reference arms H used 3). Bars are means with 95 % bootstrap CIs. FAFB values come from the dynamics track (geometric K⁺ graph, area demand, all 7 probes). "FAFB with CNS-model proxies" is FAFB re-run with synapse-count demand and the Jaccard-9NN graph. Dashed line = recovered threshold (0.8). Panel c, left pair: 30-s naive with Shiu signs vs antennal-lobe local neurons set inhibitory. Right pair: 180-s anoxia, closed-loop release vs fixed 20-s block.*

{{artifact:06758b52-e814-4468-92ed-91528037296a}}

*Fig. 2 — MCNS population rate by region (seed 0). Grey = anoxia (180-s run clipped), blue = 20-s total block, dotted line = start of the 20-s free observation. In (c), the evaluation window happened to fall in a silent gap between bursts.*

| condition | **MCNS** FII_norm [95 % CI] · dead · outcomes | FAFB-proxy | FAFB (dynamics track) |
|---|---|---|---|
| 30 s anoxia, naive | 0.050 [0.022, 0.078] · dead 0.000 · runaway 5/5 | 0.146 [0.057, 0.228] · dead 0.000 · runaway 5/5 | -0.058 [-0.180, 0.153] · runaway 5/5 |
| 30 s, 20-s block | 0.999 [0.994, 1.001] · dead 0.000 · recovered 5/5 | 1.002 [1.001, 1.003] · dead 0.000 · recovered 5/5 | 1.001 [1.000, 1.001] · recovered 5/5 |
| 180 s, naive | 0.856 [0.810, 0.901] · dead 0.009 · runaway 5/5 | 0.288 [0.086, 0.505] · dead 0.009 · runaway 5/5 | 0.084 [-0.118, 0.410] · runaway 5/5 |
| 180 s, 20-s block | 0.997 [0.992, 1.001] · dead 0.009 · recovered 5/5 | 0.971 [0.961, 0.984] · dead 0.009 · recovered 5/5 | 0.953 [0.936, 0.967] · recovered 5/5 |
| 180 s @17 °C, naive | 0.915 [0.874, 0.943] · dead 0.001 · runaway 5/5 | – | – |
| 180 s @17 °C, block | 0.992 [0.984, 0.998] · dead 0.001 · recovered 5/5 | – | – |
| 180 s @29 °C, naive | 0.483 [0.374, 0.592] · dead 0.128 · runaway 5/5 | – | – |
| 180 s @29 °C, block | 0.806 [0.738, 0.868] · dead 0.127 · recovered 3/5, structural_loss 2/5 | – | – |
| 360 s, naive | -0.033 [-0.275, 0.255] · dead 0.251 · runaway 5/5 | – | 0.085 [-0.120, 0.311] · runaway 5/5 |
| 360 s, block | 0.759 [0.691, 0.833] · dead 0.246 · structural_loss 5/5 | 0.613 [0.562, 0.671] · dead 0.246 · structural_loss 5/5 | 0.578 [0.559, 0.602] · structural_loss 5/5 |
| 360 s hub 1 %, naive | 0.235 [0.042, 0.396] · dead 0.249 · runaway 5/5 | – | -0.218 [-0.286, -0.093] · runaway 3/3 |
| 360 s hub 1 %, block | 0.960 [0.948, 0.973] · dead 0.244 · structural_loss 5/5 | 0.707 [0.678, 0.737] · dead 0.244 · structural_loss 5/5 | 0.703 [0.668, 0.743] · structural_loss 3/3 |
| 360 s hub 10 %, naive | -0.137 [-0.197, -0.071] · dead 0.223 · runaway 5/5 | – | 0.342 [0.300, 0.378] · runaway 3/3 |
| 360 s hub 10 %, block | 0.977 [0.972, 0.980] · dead 0.219 · recovered 5/5 | – | 0.906 [0.845, 0.944] · recovered 3/3 |
| 360 s random 1 %, block | 0.769 [0.704, 0.837] · dead 0.244 · structural_loss 5/5 | – | 0.568 [0.530, 0.591] · structural_loss 3/3 |
| 360 s random 10 %, block | 0.813 [0.775, 0.851] · dead 0.222 · structural_loss 5/5 | – | 0.588 [0.563, 0.637] · structural_loss 3/3 |
| 360 s whole-brain, block | 0.999 [0.994, 1.001] · dead 0.000 · recovered 5/5 | – | 1.001 [1.001, 1.001] · recovered 3/3 |
| 180 s, closed-loop release | 0.926 [0.919, 0.939] · dead 0.007 · runaway 5/5 | 0.296 [0.228, 0.346] · dead 0.007 · runaway 5/5 | – |
| sham (0 s), naive | 1.001 [1.000, 1.001] · dead 0.000 · recovered 3/3 | – | – |
| sham (0 s), block | 1.001 [1.000, 1.001] · dead 0.000 · recovered 3/3 | – | – |

### (a) Runaway attractor after naive reperfusion — REPLICATES in MCNS
* **30-s anoxia, naive: runaway in 5/5 runs, FII_norm 0.050 [0.022, 0.078]** (FAFB −0.058 [−0.180, 0.153]).
  - Late observation rate is 90.9 Hz/neuron (CB 130, VNC 137, OL 62).
  - At evaluation the network is still spontaneously active: FSA 0.084, of which VNC 0.43, central brain 0.03, optic lobe 0.04.
* **180-s naive: runaway in 5/5 by the late-rate criterion (127.8 Hz/neuron), but episodic.** The network fires recurrent all-region bursts of up to ~410 Hz/neuron separated by silent gaps of several seconds (Fig. 2c).
  - The evaluation window fell in a gap: FSA 0.001 and brain FII_norm 0.856 [0.810, 0.901].
  - VNC FII_norm is 0.045, so the VNC motor pattern is lost.
  - FII taken in a silent gap understates the pathology. The run is classed runaway because of the bursts.
* The same episodic pattern occurs in the FAFB proxy at 180 s (FII_norm 0.288 [0.086, 0.505], FSA 0.044).
* **Robust to gain:** at 0.7 × the reverberation threshold (w_syn 0.0434, 3 seeds) the 30-s naive case is still runaway in 3/3 (FII_norm 0.096).

### (b) 20-s block reset rescue — REPLICATES
* Reset minus naive, in FII_norm:
  - 30 s: +0.948 [+0.916, +0.976];
  - 180 s: +0.142 [+0.094, +0.187].
  The 180-s difference is small only because naive FII was measured in a burst gap.
* Outcomes: reset 5/5 recovered at both durations; naive 5/5 runaway.
* The sham (0 s anoxia) is recovered with or without the block, so neither the protocol nor the block alters identity on its own. The v1 instability (§6) is gone.
* SD onset is 12.5 s (FAFB 12.5 s), K⁺ peak 44.1 mM at 30 s (FAFB 45.2), and repolarisation takes 4.5 s after reperfusion (FAFB 4.5).

### (c) Hub-protection efficiency — REPLICATES (stronger in MCNS)
Setting: 360-s anoxia + 20-s block, injury_scale 0 on the protected set. The unprotected value is 0.759 [0.691, 0.833] (FAFB 0.578; MCNS loses less identity at the same ~25 % neuron death).

| protected set (MCNS, n neurons) | gain in FII_norm vs none | gain per 1,000 protected |
|---|---|---|
| top-degree 1 % (1,667) | +0.201 [+0.127, +0.278] | 0.1208 |
| top-degree 10 % (16,670) | +0.218 [+0.143, +0.293] | 0.0131 |
| random 1 % | +0.010 [-0.089, +0.114] | 0.0061 |
| random 10 % | +0.054 [-0.028, +0.140] | 0.0032 |
| whole CNS (166,700) | +0.239 [+0.165, +0.315] | 0.0014 |

* **Per protected neuron, top-1 % hub protection is 84× more efficient than whole-CNS protection** (FAFB ~30×). It recovers 84 % of the whole-CNS gain with 1 % of the neurons.
* Random sets of equal size give gains whose CIs include 0.
* Protection is useless without the block: hub-protected naive runs are runaway in 5/5.
* FAFB-proxy hub 1 %: +0.093 [+0.032, +0.152], close to the FAFB value (0.703 − 0.578 = 0.125). Proxies preserve this effect.

### (d) Temperature (180 s, block reset) — REPLICATES in direction
* 17 °C: 0.992 [0.984, 0.998]; 22 °C: 0.997 [0.992, 1.001]; 29 °C: 0.806 [0.738, 0.868] (12.7 % dead; 3/5 recovered, 2/5 structural loss).
* 17 °C minus 29 °C: +0.185 [+0.123, +0.260].
* FAFB had no 180-s cell in the temperature grid. FAFB's recoverability limit (D_max) at 29 °C is 149 s, so 180 s at 29 °C lies beyond that limit in FAFB. MCNS sits right at the boundary (0.806), i.e. it is somewhat more tolerant.
* Naive runs are runaway at every temperature. At 29 °C the naive late rate is lower (23 Hz/neuron vs 128 at 22 °C), matching FAFB's weaker hot runaway.

### VNC vs brain (new in this CNS dataset)
* After naive reoxygenation the VNC is the first region to exceed 1 Hz/neuron: 1.5 s after reperfusion, vs 2.5 s for the central brain and 4.0 s for the optic lobe (30 s anoxia; all 5 seeds identical at 0.5-s resolution; 3.5 / 4.0 / 7.0 s at 180 s).
* The central brain leads the large oscillations (Fig. 2a).
* The VNC is the most persistent part of the 30-s attractor: FSA 0.43 vs ≤0.04 elsewhere at evaluation.
* VNC motor output pattern (LC4 probe): FII_norm −0.231 after 30-s naive and 0.998 after the block.
* VNC dead fraction tracks the whole CNS (360 s: ~25 %).

### Sign-corrected variant (antennal-lobe local neurons inhibitory; 219 MCNS / 147 FAFB neurons flipped)
* **30-s naive:** the correction largely abolishes identity loss.
  - MCNS FII_norm 1.000 [0.996, 1.002]; difference vs Shiu signs +0.949 [+0.918, +0.977].
  - FAFB-proxy 0.859; difference +0.714 [+0.632, +0.801].
  - A low-level oscillation (~3–5 Hz/neuron bursts) remains, so the runaway criterion (>0.5 Hz) is still met.
* **180-s naive:** the runaway persists and becomes continuous.
  - MCNS −0.135, FSA 0.897; difference -0.991 [-1.119, -0.883].
  - FAFB-proxy −0.125; difference -0.413 [-0.622, -0.192].
* **Block reset still rescues** (MCNS 0.989–1.000; proxy 0.965–1.002).
* Interpretation (model-level): the AL-LN sign artefact flagged by the index track drives the short-anoxia runaway. After ≥180 s anoxia the attractor no longer depends on it. The core finding "naive reoxygenation → runaway, block → rescue" survives the sign correction for long insults.

### Closed-loop reset (block until max K⁺ < 6 mM and mean ATP e > 0.9)
* The criterion was met after **11 s** in every run (K⁺ 5.55 mM, e 0.901).
* MCNS FII_norm 0.926 [0.919, 0.939], but residual activity remains (11 Hz/neuron, FSA 0.055; runaway 5/5). Difference vs 20-s block: -0.071 [-0.081, -0.058].
* FAFB-proxy: 0.296; difference -0.674 [-0.742, -0.624].
* **Normalised K⁺ and ATP are not a sufficient release signal.** The network still re-ignites when released at that point. A fixed ≥20-s block is needed (FAFB: 10 s was inconsistent and 20 s reliable). Hidden slow state not captured by K⁺/ATP is implicated; untested which.

## 3. Verdicts
| test | verdict (MCNS; BANC not analysed) |
|---|---|
| (a) runaway attractor after naive reperfusion | **replicates** (5/5 at 30, 180, 360 s and 17/22/29 °C; continuous at 30 s, episodic at ≥180 s) |
| (b) 20-s block reset rescue | **replicates** (5/5 recovered at 30 and 180 s) |
| (c) hub-protection efficiency | **replicates, stronger** (~84× per neuron vs FAFB ~30×; random controls null) |
| (d) temperature effect | **replicates in direction** (17 °C − 29 °C = +0.185); magnitude not comparable to FAFB at matched D |
| VNC vs brain | the VNC crosses 1 Hz/neuron first and is the most persistent part of the attractor |

## 4. Limitations
* **Single CNS dataset.** BANC results exist on the Volume but are unanalysed; MCNS is male, FAFB female, so sex and dataset are confounded.
* **Seeds:** 5 per condition (3 for sham and sensitivity). Seeds vary Poisson input and SD-threshold jitter only; the connectome is fixed.
* **Proxies and probes:** proxy substitutions (synapse-count demand, connection-based K⁺ graph) differ from FAFB. The FAFB-proxy control shows the core effects survive them, but quantitative values shift (e.g. proxy 30-s naive 0.146 vs FAFB −0.058). Sugar→MN9 is silent in MCNS and was excluded from FII; VNC FII rests on a single probe (LC4).
* **FII timing:** a single 500-ms evaluation window cannot distinguish an episodic attractor in a silent gap from recovery. FII_norm for 180-s naive (0.856) must be read together with the outcome class.
* **Model constants:** all metabolic constants are phenomenological and taken from the FAFB calibration (τ_ATP 30 s, J50 600 s, q10_injury 2.5, speedup 50). They were not re-calibrated per connectome.
* **Compute accounting:** BANC job cost is unknown (status probe lost; ≤ $1.76 by the container cap). Logged in `research/compute_ledger.csv`.

## 5. Files
* **Code:** `flysim_cns.py`, `build_cns.py`, `calib.py`, `rsweep.py`, `analyze_rep.py` (local pilot); the v4 analysis cell logic is reproduced in `tables/` outputs.
* **Tables:**
  - `tables/v4_all_runs.csv` (211 runs);
  - `tables/v4_condition_summary.csv`;
  - `tables/v4_contrasts.json` (difference, CI low, CI high[, per-1,000 gain, n protected]);
  - `tables/v4_meta.json`;
  - `tables/fafb_reference.csv`.
* **Raw:** `raw/v4/` (harvested json + light npz), `raw/calib/`, `raw/local_v3/`.
* **Figures:** `figures/fig_rep1_identity.png`, `figures/fig_rep2_region_traces.png`.
* **Modal job scripts:** `modal_jobs/`.
