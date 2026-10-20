"""Figures for the Structure track. Run inside the python kernel after figure-style is loaded
(apply_figure_style, panel_letter available) via exec(open(...).read())."""
# --- repo bootstrap: make src/ importable and resolve data/output roots (see src/fbr_paths.py) ---
import os as _os, sys as _sys
_REPO = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), _os.pardir, _os.pardir))
for _d in ("src",):
    _p = _os.path.join(_REPO, _d)
    if _os.path.isdir(_p) and _p not in _sys.path:
        _sys.path.insert(1, _p)
import fbr_paths as _P  # noqa: E402
# --- end bootstrap ---
import os
try:
    apply_figure_style  # provided by the authoring kernel when exec()-ed there
except NameError:  # standalone fallback
    def apply_figure_style():
        pass

    def panel_letter(ax, letter):
        ax.text(-0.12, 1.05, letter, transform=ax.transAxes, fontweight="bold", va="bottom")
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
ROOT = os.path.join(_P.OUT, 'structure'); T = os.path.join(ROOT, 'tables'); FIG = os.path.join(ROOT, 'figures')
DSS = ['FAFB', 'BANC', 'MCNS', 'MANC']
DSLAB = {'FAFB': 'FAFB · female brain', 'BANC': 'BANC · female CNS', 'MCNS': 'MCNS · male CNS',
         'MANC': 'MANC · male VNC'}
REG = ['random', 'metabolic_area', 'metabolic_syn', 'focal', 'degree_matched_random', 'degree_finegrid']
REGLAB = {'random': 'random', 'metabolic_area': 'metabolic (membrane area)', 'metabolic_syn': 'metabolic (synapse count)',
          'degree': 'hub-targeted (degree)', 'degree_finegrid': 'hub-targeted (degree)', 'degree_matched_random': 'degree-bin-matched random', 'focal': 'focal neuropil lesion'}
REGCOL = {'random': '#4d4d4d', 'metabolic_area': '#e69f00', 'metabolic_syn': '#b35806', 'degree': '#0072b2', 'degree_finegrid': '#0072b2',
          'degree_matched_random': '#56b4e9', 'focal': '#cc79a7'}
cur = pd.read_csv(os.path.join(T, 'curves_normalised_all.csv.gz'))
fits = pd.read_csv(os.path.join(T, 'fits_all.csv'))
rat = pd.read_csv(os.path.join(T, 'ratios_all.csv'))
apply_figure_style()


def band(ax, sub, m, col, lab, ls='-', lw=1.4, alpha=0.18):
    g = sub.groupby('f')[m]
    mu = g.mean(); lo = g.quantile(0.025); hi = g.quantile(0.975)
    ax.fill_between(mu.index, lo.values, hi.values, color=col, alpha=alpha, lw=0)
    ax.plot(mu.index, mu.values, color=col, lw=lw, ls=ls, marker='o', ms=2.2, label=lab)


# ---------- Figure 1: primary metric, all regimes, all datasets ----------
fig, axs = plt.subplots(1, 4, figsize=(7.2, 2.5), sharey=True, gridspec_kw=dict(wspace=0.12))
for ax, ds in zip(axs, DSS):
    sub = cur[(cur.dataset == ds) & (cur.graph == 'thr5')]
    for rg in REG:
        s = sub[sub.regime == rg]
        if len(s): band(ax, s, 'reachN_L6_n', REGCOL[rg], REGLAB[rg])
    ax.set_title(DSLAB[ds], loc='left', pad=8)
    ax.axhline(0.5, color='0.75', lw=0.6, ls=':', zorder=0)
    ax.set_xlim(-0.02, 0.92); ax.set_ylim(-0.03, 1.08)
axs[0].set_ylabel('sensory→motor reach (≤6 hops)\nrelative to intact')
h, l = axs[0].get_legend_handles_labels()
fig.supxlabel('fraction of interneurons removed, f', y=-0.06)
fig.legend(h, l, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.27), handlelength=1.4, columnspacing=1.2)
for ax, L in zip(axs, 'abcd'): panel_letter(ax, L)
f1 = os.path.join(FIG, 'fig1_reach_vs_damage.png'); fig.savefig(f1, bbox_inches='tight'); fig.savefig(f1.replace('.png', '.pdf'), bbox_inches='tight')
plt.close(fig)

# ---------- Figure 2: targeted / random ratio (forest) ----------
r = rat[(rat.graph == 'thr5') & (rat.metric == 'reachN_L6')].copy()
order_rg = ['degree_finegrid', 'degree_matched_random', 'metabolic_syn', 'metabolic_area', 'focal']
fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.0), sharey=True, gridspec_kw=dict(wspace=0.08))
y = 0; ypos = {}
for rg in order_rg:
    for ds in DSS:
        if len(r[(r.regime == rg) & (r.dataset == ds)]): ypos[(rg, ds)] = y; y += 1
    y += 0.6
for ax, (col, lo, hi, ttl) in zip(axs, [('ratio_fc', 'ratio_fc_lo', 'ratio_fc_hi', 'logistic midpoint $f_c$'),
                                         ('ratio_D50', 'ratio_D50_lo', 'ratio_D50_hi', 'D50')]):
    for (rg, ds), yy in ypos.items():
        q = r[(r.regime == rg) & (r.dataset == ds)].iloc[0]
        open_ = (col == 'ratio_fc' and not q.random_fc_in_range)
        if col == 'ratio_D50' and not np.isfinite(q[col]):
            ub = q['ratio_D50_upper_bound_if_random_censored']
            if np.isfinite(ub):
                ax.plot(ub, yy, marker=4, ms=5, color=REGCOL[rg], mfc='white')   # '<' bound: random never fell to 50%
            continue
        ax.errorbar(q[col], yy, xerr=[[max(q[col] - q[lo], 0)], [max(q[hi] - q[col], 0)]], fmt='o', ms=3.4,
                    color=REGCOL[rg], mfc='white' if open_ else REGCOL[rg], ecolor=REGCOL[rg], elinewidth=1, capsize=1.5)
    ax.axvline(1.0, color='0.4', lw=0.8); ax.axvline(0.7, color='0.4', lw=0.8, ls='--')
    ax.set_xlabel(f'targeted / random ({ttl})'); ax.set_xlim(0, 1.15)
axs[0].set_yticks(list(ypos.values())); axs[0].set_yticklabels([f'{REGLAB[rg]} · {ds}' for (rg, ds) in ypos])
axs[0].invert_yaxis()
axs[0].text(0.69, -0.9, 'H2 cut-off 0.7', ha='right', va='bottom', fontsize=6, color='0.3')
axs[1].text(0.99, -0.9, 'no effect', ha='right', va='bottom', fontsize=6, color='0.3')
for ax, L in zip(axs, 'ab'): panel_letter(ax, L)
f2 = os.path.join(FIG, 'fig2_targeted_vs_random_ratio.png'); fig.savefig(f2, bbox_inches='tight'); fig.savefig(f2.replace('.png', '.pdf'), bbox_inches='tight')
plt.close(fig)

# ---------- Figure 3: FAFB metric hierarchy (random) + controls ----------
from matplotlib.lines import Line2D
MET = [('influence_L6_n', 'weighted influence ≤6 hops', '#756bb1'), ('gscc_n', 'giant strongly connected component', '#636363'),
       ('reachN_L4_n', 'reach ≤4 hops', '#9ecae1'), ('reachN_L6_n', 'reach ≤6 hops (primary)', '#3182bd'),
       ('reachN_Linf_n', 'reach, any length (≤30 hops)', '#08519c'), ('redund_L6_n', 'path redundancy (log walks ≤6)', '#e6550d'),
       ('reachG_L6_n', 'modality-level reach ≤6 hops', '#31a354')]
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw=dict(wspace=0.32))
sub = cur[(cur.dataset == 'FAFB') & (cur.graph == 'thr5') & (cur.regime == 'random')]
for m, lab, col in MET:
    band(axs[0], sub.dropna(subset=[m]), m, col, lab)
axs[0].set_title('FAFB, random loss', loc='left', pad=8)
axs[0].set_ylabel('relative to intact')
axs[0].legend(frameon=False, fontsize=5.5, loc='lower left', handlelength=1.2)
LS = {'thr5': '-', 'null0': '--', 'thr10': ':', 'thr1': '-.'}
for ax, ds in zip(axs[1:], ['FAFB', 'MANC']):
    sub = cur[(cur.dataset == ds)]
    for gr in LS:
        for rg in ['random', 'degree']:
            s_ = sub[(sub.graph == gr) & (sub.regime == rg)]
            if len(s_): band(ax, s_, 'reachN_L6_n', REGCOL[rg], None, ls=LS[gr], lw=1.1, alpha=0.08)
    ax.set_title(f'{ds} controls', loc='left', pad=8)
    ax.set_ylabel('reach ≤6 hops, relative to intact')
for ax in axs: ax.set_xlabel('fraction of interneurons removed, f')
hs = [Line2D([], [], color=REGCOL['random'], lw=1.4), Line2D([], [], color=REGCOL['degree'], lw=1.4),
      Line2D([], [], color='k', ls='-'), Line2D([], [], color='k', ls='--'), Line2D([], [], color='k', ls=':'), Line2D([], [], color='k', ls='-.')]
ls_ = ['random loss', 'hub-targeted loss', 'real graph, ≥5 synapses', 'degree-preserving rewired null', 'real graph, ≥10 synapses', 'real graph, ≥1 synapse (MANC)']
fig.legend(hs, ls_, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.62, -0.16), handlelength=2.2)
for ax, L in zip(axs, 'abc'): panel_letter(ax, L)
f3 = os.path.join(FIG, 'fig3_metrics_and_controls.png'); fig.savefig(f3, bbox_inches='tight'); fig.savefig(f3.replace('.png', '.pdf'), bbox_inches='tight')
plt.close(fig)

# ---------- Figure 4: keystone analysis ----------
from scipy.stats import spearmanr
K = {}
for ds in DSS:
    d = pd.read_csv(os.path.join(T, f'{ds}_keystone.csv')).drop_duplicates('type'); d['key'] = d.type.astype(str).str.lower(); K[ds] = d
fig, axs = plt.subplots(1, 4, figsize=(7.4, 2.3), gridspec_kw=dict(wspace=0.55, width_ratios=[1.3, 1, 1, 1]))
ax = axs[0]
rng = np.random.default_rng(0)
for i, ds in enumerate(DSS):
    v = 100 * K[ds].rel_reachN_L6.values
    ax.scatter(i + rng.uniform(-0.18, 0.18, v.size), v, s=5, color='#3182bd', alpha=0.5, lw=0)
    ax.plot([i - 0.25, i + 0.25], [np.median(v)] * 2, color='k', lw=1)
ax.set_xticks(range(4)); ax.set_xticklabels(DSS)
ax.set_ylabel('reach lost when one whole\ncell type is removed (%)')
ax.set_title('single-type removal', loc='left', pad=8)
pairs = [('BANC', 'MCNS'), ('FAFB', 'BANC'), ('MCNS', 'MANC')]
for ax, (a_, b_) in zip(axs[1:], pairs):
    m = K[a_].set_index('key')[['rel_influence_L6']].join(K[b_].set_index('key')[['rel_influence_L6']], lsuffix='_a', rsuffix='_b', how='inner')
    x = 100 * m.rel_influence_L6_a; y_ = 100 * m.rel_influence_L6_b
    rho, pv = spearmanr(x, y_)
    ax.scatter(x, y_, s=6, color='#756bb1', alpha=0.7, lw=0)
    lim = max(x.max(), y_.max()) * 1.08
    ax.plot([0, lim], [0, lim], color='0.7', lw=0.6, ls=':')
    ax.set_xlim(-0.02 * lim, lim); ax.set_ylim(-0.02 * lim, lim)
    ax.set_xlabel(f'{a_} influence lost (%)'); ax.set_ylabel(f'{b_} influence lost (%)')
    ax.text(0.96, 0.04, f'ρ = {rho:.2f}\nn = {len(m)} types', transform=ax.transAxes, va='bottom', ha='right', fontsize=6)
for ax, L in zip(axs, 'abcd'): panel_letter(ax, L)
f4 = os.path.join(FIG, 'fig4_keystone_types.png'); fig.savefig(f4, bbox_inches='tight'); fig.savefig(f4.replace('.png', '.pdf'), bbox_inches='tight')
plt.close(fig)
print(f1, f2, f3, f4)
