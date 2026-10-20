"""flysim_cns — loaders that turn MCNS v1.0 / BANC v888 into flysim_dyn.Model objects (replication track).

Does not modify flysim.py / flysim_dyn.py; it only builds a flysim_dyn.Model from the caches written by
build_cns.py.  Conventions (declared in REPLICATION_RESULTS.md):
  sign      presynaptic predicted NT (BANC: verified NT used only where no prediction exists);
            GABA, GLUT -1 (Shiu), HIST -1 (histamine-gated Cl- channels; Xenova demo uses 0 — see hist_sign),
            ACH +1, DA/SER/OCT/TYR +1 (Shiu modulatory), unknown +1 (Shiu/FAFB convention).
  demand    metabolic demand proxy = total synapse count (no size/area data in either release),
            used through flysim demand_base: clip((x/median)^0.5, 0.25, 4).
  adjacency connection-based: Jaccard overlap of the neuron sets of two neuropils, 9 nearest neighbours per
            neuropil (union-symmetrised), row-max normalised, rescaled to FAFB geometric median row sum (2.62).
            (The dense shared-neuron B^T B adjacency makes the explicit K+ diffusion step unstable:
            k_diff*dt*lambda_max = 2.9-4.8 > 2; jknn9 gives 0.69-0.78 vs FAFB geometric 0.89. On FAFB the jknn9
            edges recover geometric contacts with precision 0.59 / recall 0.68, chance 0.135.)
  compartments  neuropils with <0.1 % of synapse mass merged into <region>_MISC pools.
"""
import os, sys
import numpy as np, pandas as pd, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (HERE, os.path.join(HERE, "..", "..", "dynamics"), os.path.join(os.environ.get("RESEARCH_ROOT", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs")), 'dynamics')):
    if os.path.exists(os.path.join(_p, "flysim_dyn.py")) and _p not in sys.path:
        sys.path.insert(0, _p)
import flysim_dyn as fs

CNS_CACHE = os.environ.get("CNS_CACHE", os.path.join(os.environ.get("RESEARCH_ROOT", os.path.join(os.path.dirname(HERE), "outputs")), "followup", "replication", "cache"))

# calibrated w_syn (mV/synapse) per connectome — filled from calibration (see REPLICATION_RESULTS.md §2)
CAL_WSYN = {"FAFB": 0.185, "MCNS": 0.0538, "BANC": 0.2227}   # = 0.8667 x any-probe reverberation threshold (FAFB-matched)


def check_stability(adj, k_diff=0.5, dts=0.25):
    L = adj - np.diag(adj.sum(1)); lam = float(np.abs(np.linalg.eigvalsh(0.5 * (L + L.T))).max())
    if k_diff * dts * lam >= 1.5:
        raise ValueError(f"K+ diffusion step unstable/marginal: k_diff*dt*lambda = {k_diff*dts*lam:.2f}")
    return lam


def sign_vector(nt, hist_sign=-1, unknown_sign=1, modulatory_sign=1):
    nt = pd.Series(nt).fillna("UNKNOWN").astype(str).str.upper().values
    s = np.ones(len(nt), np.int8)
    s[np.isin(nt, ["GABA", "GLUT"])] = -1
    s[nt == "HIST"] = hist_sign
    s[np.isin(nt, ["DA", "SER", "OCT", "TYR"])] = modulatory_sign
    s[nt == "UNKNOWN"] = unknown_sign
    return s


ALLN_CLASS = {"MCNS": "ALLN", "BANC": "antennal_lobe_local_neuron", "FAFB": "ALLN"}


def alln_mask(neurons, ds):
    """antennal-lobe local neurons (GABAergic in real flies; many predicted ACh/Glu/unknown)."""
    return (neurons["class"].astype(str) == ALLN_CLASS[ds]).values


def load_cns(ds, w_syn=None, hist_sign=-1, unknown_sign=1, cache=None, alln_inhib=False):
    """ds: 'MCNS' | 'BANC'. Returns flysim_dyn.Model. demand_mode for metab must be 'syn_total'."""
    import time
    t0 = time.time()
    cache = cache or CNS_CACHE
    t = pd.read_csv(os.path.join(cache, f"cns_{ds}_neurons.csv.gz"), low_memory=False)
    A = sp.load_npz(os.path.join(cache, f"cns_{ds}_csr.npz")).tocsr(); A.sort_indices()
    N = A.shape[0]
    sign = sign_vector(t.nt.values, hist_sign, unknown_sign)
    n_flip = 0
    if alln_inhib:
        am = alln_mask(t, ds); n_flip = int((sign[am] > 0).sum()); sign[am] = -1
    deg = np.diff(A.indptr)
    w = A.data.astype(np.float32) * np.repeat(sign, deg).astype(np.float32)
    Mraw = sp.load_npz(os.path.join(cache, f"cns_{ds}_neuron_neuropil.npz")).tocsr().astype(np.float32)
    rs = np.asarray(Mraw.sum(1)).ravel()
    Mn = (sp.diags((1.0 / np.where(rs > 0, rs, 1)).astype(np.float32)) @ Mraw).tocsr().astype(np.float32)
    nps = pd.read_csv(os.path.join(cache, f"cns_{ds}_neuropils.csv")).neuropil.tolist()
    adj = np.load(os.path.join(cache, f"cns_{ds}_adjacency_jknn9.npy"))
    check_stability(adj)
    if w_syn is None:
        w_syn = CAL_WSYN.get(ds, 0.185)
    m = fs.Model(N=N, topk=fs._topk_compartments(Mraw, 3), w_syn=float(w_syn), MnT=Mn.T.tocsr().astype(np.float32),
                 indptr=A.indptr.astype(np.int64), indices=A.indices.astype(np.int32), weight=w, neurons=t, sign=sign,
                 neuropils=nps, Mn=Mn, comp_n=np.asarray(Mn.sum(0)).ravel().astype(np.float64), adj=adj,
                 adj_source="connection-based Jaccard 9-NN (jknn9)",
                 info=dict(connectome=ds, n_edges=int(A.nnz), total_synapses=int(A.data.astype(np.float64).sum()),
                           w_syn=float(w_syn), hist_sign=hist_sign, unknown_sign=unknown_sign,
                           demand_proxy="syn_total", alln_inhib=alln_inhib, alln_flipped=n_flip,
                           load_s=round(time.time() - t0, 1)))
    return m


def fafb_proxy(w_syn=None, alln_inhib=False):
    """FAFB princeton5 with the CNS-model substitutions (synapse-count demand + connection-based adjacency):
    isolates the effect of the proxies from the effect of the connectome."""
    m = fs.load_model("princeton5", w_syn=w_syn)
    if alln_inhib:
        am = alln_mask(m.neurons, "FAFB"); flip = am & (m.sign > 0)
        m.info["alln_flipped"] = int(flip.sum())
        m.sign = m.sign.copy(); m.sign[flip] = -1
        rowflip = np.repeat(flip, np.diff(m.indptr))
        m.weight = np.where(rowflip, -np.abs(m.weight), m.weight).astype(np.float32)
    A = sp.csr_matrix((np.abs(m.weight), m.indices, m.indptr), shape=(m.N, m.N))
    t = m.neurons.copy()
    t["syn_out"] = np.asarray(A.sum(1)).ravel(); t["syn_in"] = np.asarray(A.sum(0)).ravel()
    t["syn_total"] = t.syn_out + t.syn_in
    m.neurons = t
    m.adj = np.load(os.path.join(CNS_CACHE, "fafb_adjacency_jknn9.npy")); check_stability(m.adj)
    m.adj_source = "connection-based Jaccard 9-NN (proxy control)"
    m.info["connectome"] = "FAFBproxy"; m.info["demand_proxy"] = "syn_total"
    return m


def _sel(m, mask):
    return np.flatnonzero(np.asarray(mask, bool))


def battery(m, ds):
    """Probe battery analogous to the FAFB battery. Returns (probes, readouts, named, report).
    probes: name -> (idx, rate_Hz); readouts: dict of index arrays ('brain' = DN + brain MNs (FAFB-equivalent),
    'vnc' = VNC motor + efferent, 'all' = union); named: readout cell types."""
    t = m.neurons
    ty = t.primary_type.astype(str); side = t.side.astype(str); sc = t.super_class.astype(str)
    sub = t.sub_class.astype(str)
    reg = t.region.fillna("NA").astype(str) if "region" in t.columns else pd.Series(["CB"] * len(t))
    rep = {}
    P = {}
    if ds == "FAFB":
        P["sugar"] = (m.idx(fs.SHIU_SUGAR_GRN), 150.0)
        P["bitter"] = (m.select(None, **{"class": "gustatory", "sub_class": "bitter"}), 150.0)
        P["JO_aud"] = (m.select(None, **{"class": "mechanosensory", "sub_class": "auditory"}), 150.0)
        P["JO_wg"] = (m.select(None, **{"class": "mechanosensory", "sub_class": "wind_gravity"}), 150.0)
    else:
        P["sugar"] = (_sel(m, ty.str.match(r"^LB3[a-d]?$") & (side == "left")), 150.0)
        P["bitter"] = (_sel(m, ty.str.match(r"^LB1[a-e]$")), 150.0)
        if ds == "MCNS":
            P["JO_aud"] = (_sel(m, sub == "auditory"), 150.0)
            P["JO_wg"] = (_sel(m, sub == "wind_gravity"), 150.0)
        else:
            P["JO_aud"] = (_sel(m, sub.isin(["johnstons_organ_A_neuron", "johnstons_organ_B_neuron"])), 150.0)
            P["JO_wg"] = (_sel(m, sub.isin(["johnstons_organ_C_neuron", "johnstons_organ_D_neuron",
                                            "johnstons_organ_E_neuron"])), 150.0)
    P["LC4"] = (_sel(m, ty == "LC4"), 100.0)
    P["LC9"] = (_sel(m, ty == "LC9"), 100.0)
    P["turnL"] = (np.union1d(_sel(m, (ty == "LC9") & (side == "left")), _sel(m, (ty == "DNa02") & (side == "left"))), 180.0)
    for k, (idx, r) in P.items():
        rep[k] = dict(n=int(len(idx)), types=pd.Series(ty.values[idx]).value_counts().head(8).to_dict(), rate=r)
    if ds == "FAFB":
        mn9 = m.idx([fs.SHIU_MN9])
        brain = np.union1d(m.select(None, super_class=["descending", "motor"]), mn9)
        vnc = np.zeros(0, np.int64)
    else:
        mn9 = _sel(m, ty == "MN9")
        brain = np.union1d(_sel(m, sc == "descending"), _sel(m, (sc == "motor") & (reg != "VNC")))
        brain = np.union1d(brain, mn9)
        vnc = _sel(m, sc.isin(["motor", "efferent"]) & (reg == "VNC"))
    ro = dict(brain=brain, vnc=vnc, all=np.union1d(brain, vnc))
    named = {"MN9": mn9}
    for tyn in ["DNp01", "DNg100", "DNp09", "DNa02", "DNa01", "DNg97"]:
        for sd in ["left", "right"]:
            named[f"{tyn}_{sd[0].upper()}"] = _sel(m, (ty == tyn) & (side == sd))
    named["MDN"] = _sel(m, ty == "MDN")
    rep["readouts"] = {k: int(len(v)) for k, v in ro.items()}
    rep["named"] = {k: int(len(v)) for k, v in named.items()}
    return P, ro, named, rep
