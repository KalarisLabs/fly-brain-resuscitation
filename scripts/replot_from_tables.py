"""Table-derived re-plots (written for this repository; not the original interactive figure code).

Reads the staged tables under $RESEARCH_ROOT (run scripts/reproduce_figures.py first, which stages results/) and
writes outputs/replots/{R1_brain_death_map,R2_protocol_ranking,R3_lhs_robustness,R4_cross_connectome}.png.
Values are plotted exactly as stored; nothing is recomputed.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
sys.path.insert(0, os.path.join(ROOT, "src"))
import fbr_paths as P  # noqa: E402

OUT = os.path.join(P.OUT, "replots")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 200})


def tab(*parts):
    return pd.read_csv(os.path.join(P.OUT, *parts))


def save(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print("wrote", p)


# R1: brain-death map (dynamics B table): FII_norm by insult temperature x anoxia duration, per reperfusion arm
B = tab("dynamics", "tables", "B_brain_death_map.csv")
arms = sorted(B.reperf.unique())
fig, axs = plt.subplots(1, len(arms), figsize=(3.4 * len(arms), 2.6), squeeze=False, gridspec_kw=dict(wspace=0.15))
for ax, a in zip(axs[0], arms):
    M = B[B.reperf == a].pivot_table(index="T_ins", columns="D", values="FII_norm")
    im = ax.imshow(M.values, aspect="auto", origin="lower", vmin=-0.3, vmax=1.0, cmap="viridis")
    ax.set_xticks(range(M.shape[1]), [f"{c:g}" for c in M.columns], rotation=90)
    ax.set_yticks(range(M.shape[0]), [f"{r:g}" for r in M.index])
    ax.set_xlabel("anoxia duration D (model s)"); ax.set_title(f"reperfusion: {a}")
axs[0][0].set_ylabel("insult temperature (C)")
fig.colorbar(im, ax=axs[0].tolist(), label="FII_norm")
save(fig, "R1_brain_death_map.png")

# R2: reperfusion-protocol ranking (dynamics C table), FII_norm at 30 / 180 / 360 s anoxia with bootstrap CI
C = tab("dynamics", "tables", "C_resuscitation_ranking.csv").sort_values("FII_norm_mean")
fig, ax = plt.subplots(figsize=(4.2, 0.18 * len(C) + 0.8))
y = np.arange(len(C))
for k, (D, c) in enumerate([(30, "#0072b2"), (180, "#e69f00"), (360, "#cc79a7")]):
    m, lo, hi = C[f"FII_norm_D{D}"], C[f"lo_D{D}"], C[f"hi_D{D}"]
    ax.errorbar(m, y + (k - 1) * 0.22, xerr=[m - lo, hi - m], fmt="o", ms=2.5, lw=0.8, color=c, label=f"D = {D} s")
ax.set_yticks(y, C.protocol); ax.set_xlabel("FII_norm (mean, 95% CI)"); ax.axvline(1, lw=0.5, color="k")
ax.legend(frameon=False, loc="lower right")
save(fig, "R2_protocol_ranking.png")

# R3: Latin-hypercube robustness (mechanism L table): fraction of parameter sets with runaway / rescue
L = tab("followup", "mechanism", "tables", "L_fractions.csv")
fig, ax = plt.subplots(figsize=(3.6, 2.4))
lab = [f"D{d} {r}" for d, r in zip(L.D, L.reperf)]
x = np.arange(len(L))
for off, (col, lo, hi, c, nm) in zip([-0.18, 0.18], [("frac_sets_runaway", "run_lo", "run_hi", "#b2182b", "runaway"),
                                                     ("frac_sets_rescued", "res_lo", "res_hi", "#2166ac", "rescued")]):
    ax.bar(x + off, L[col], 0.36, color=c, label=nm,
           yerr=[(L[col] - L[lo]).clip(lower=0), (L[hi] - L[col]).clip(lower=0)], capsize=2)
ax.set_xticks(x, lab); ax.set_ylabel(f"fraction of {int(L.n_sets.max())} parameter sets"); ax.set_ylim(0, 1.05)
ax.legend(frameon=False)
save(fig, "R3_lhs_robustness.png")

# R4: cross-connectome dynamic replication (replication v4 table): FII_norm (responsive probes) per condition
V = tab("followup", "replication", "tables", "v4_condition_summary.csv")
V = V.dropna(subset=["FIIn_resp"]).sort_values(["label", "condition"])
fig, ax = plt.subplots(figsize=(4.6, 0.16 * len(V) + 0.8))
y = np.arange(len(V))
ax.errorbar(V.FIIn_resp, y, xerr=[(V.FIIn_resp - V.FIIn_resp_lo).clip(lower=0), (V.FIIn_resp_hi - V.FIIn_resp).clip(lower=0)],
            fmt="o", ms=2.5, lw=0.8, color="#333333")
ax.set_yticks(y, [f"{a} | {b}" for a, b in zip(V.label, V.condition)], fontsize=6)
ax.axvline(1, lw=0.5, color="k"); ax.axvline(0, lw=0.5, color="k", ls=":"); ax.set_xlabel("FII_norm (responsive probes, 95% CI)")
save(fig, "R4_cross_connectome.png")
