# flysim_mech: copy of research/dynamics/flysim_dyn.py + optional excitotoxic reperfusion term (k_exc, r0_exc; default off -> identical numerics)
"""
flysim_dyn — Dynamics-track copy of flysim.py with additions: demand_scale & per-neuron atp_floor/injury_scale interventions, block_mask, coma_arrest, q10_injury, inh_gain/exc_gain.
flysim — whole-brain spiking + metabolic/ionic simulator on the FAFB v783 Drosophila connectome.

Fast layer  : exact re-implementation of the Shiu et al. 2024 (Nature 634:210) LIF model
              (philshiu/Drosophila_brain_model, model.py, brian2 'linear' integration), event-driven
              in numpy/scipy.  dt = 0.1 ms.
Slow layer  : optional phenomenological metabolic / ionic layer updated every `dt_slow_ms` (1 ms):
              per-neuron ATP, Na/K-pump failure -> resting-potential drift, per-neuropil extracellular
              K+ pools with diffusion, regenerative K+ release (spreading depolarization), depolarization
              block, injury integral -> stochastic permanent death, reversible synaptic failure.

Everything is pure numpy/scipy (no numba/brian2) so the step functions can be ported to torch/cupy.
See README.md in this folder for the API.
"""
from __future__ import annotations
import os, json, time, math
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
import scipy.sparse as sp

CACHE = os.environ.get("FLYSIM_CACHE", os.path.join(os.environ.get("RESEARCH_ROOT", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs")), "cache"))

# ----------------------------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------------------------
LIF_DEFAULTS = dict(            # Shiu et al. 2024 model.py default_params (verified against GitHub source)
    dt=0.1,                     # ms (brian2 default clock)
    v_rest=-52.0,               # mV  v_0
    v_reset=-52.0,              # mV  v_rst
    v_th=-45.0,                 # mV
    tau_m=20.0,                 # ms  t_mbr
    tau_syn=5.0,                # ms  tau
    t_ref=2.2,                  # ms  t_rfc (v and g frozen while refractory; Poisson targets have t_ref=0)
    t_delay=1.8,                # ms  t_dly
    w_syn=0.275,                # mV per synapse (added to g on spike arrival)
    f_poi=250.0,                # Poisson input weight = f_poi * w_syn = 68.75 mV, added directly to v
    active_tol=1e-5,            # mV; neurons within tol of rest (and g<tol) are snapped to rest and skipped (0 = exact)
)

METAB_DEFAULTS = dict(
    dt_slow_ms=5.0,             # slow-variable update interval (ms); 1.0 = finer but ~3x slower
    speedup=1.0,                # multiplies all slow rates (time compression for long protocols; 1 = real time)
    # --- ATP:  da/dt = [m*b_i*Qs*S*(1-a)/(K_adp+1-a) - d_i*a]/tau_atp ;  e_i = a_i/a0 (normalized energy)
    #     saturating (ADP-limited, finite-reserve) production; K_adp chosen so that a = a0 at rest, S=1
    tau_atp_s=20.0,             # s; anoxic e-fold decay time of a median neuron at rest
    a0=0.95,                    # baseline ATP fraction
    reserve=4.0,                # max production / resting demand (metabolic reserve)
    q10_supply=1.5,             # temperature factor on maximal production
    demand_mode="area",         # 'area' (Attwell-Laughlin membrane area), 'length', 'surface_to_volume', 'uniform'
    demand_beta=0.5,            # base_i = clip((x_i/median x)^beta, base_min, base_max)
    base_min=0.25, base_max=4.0,
    r_half_hz=100.0,            # spike demand: d_i = base_i*Q*(1 + r_i/r_half + k_K*(K_i-K0)/K0)
    k_K=0.2,                    # extra pump demand per unit relative K+ elevation
    tau_rate_s=0.1,             # rate estimator time constant
    T0=22.0, q10_demand=2.5, q10_pump=1.8,   # temperature (deg C) factors
    # --- pump failure -> resting potential drift
    e_ok=0.75,                  # dead-zone: no drift while e >= e_ok
    dV_atp_max=40.0, n_atp=1.5, # dV_atp = dV_max * ((e_ok - e)/e_ok)^n
    # --- extracellular K+ (per neuropil compartment), mM
    K0=5.0, K_max=80.0,
    g_K=0.6,                    # fraction of resting conductance that is K+: dV_K = g_K*RT/F*ln(K_i/K0)
    RT_F=25.6,                  # mV at 22 C (RT/F)
    k_spk=0.05,                 # mM/s per Hz of mean compartment firing rate
    k_rel=60.0,                 # mM/s regenerative release when all neurons of a compartment depolarized
    tau_rel_depl_s=8.0,         # release store depletes while releasing (terminates SD in healthy tissue)
    tau_rel_rec_s=30.0,         # store recovers (pump dependent)
    k_clr=1.0,                  # 1/s pump/glial K+ clearance at full energy
    glia_w=0.5,                 # clearance = k_clr*Qp*(glia_w*e_glia + (1-glia_w)*e_neur)
    k_leak=2.0,                 # mM/s passive K+ leak when compartment pump capacity is 0
    pump_ok=0.9,                # leak dead-zone: no passive leak while compartment pump capacity >= pump_ok
    k_diff=0.5,                 # 1/s inter-neuropil diffusion (scaled adjacency, row max = 1)
    # --- spreading depolarization / depolarization block
    V_sd=-30.0, k_sd=2.0, sd_jitter=2.0,     # dep_i = sigmoid((vr_eff - V_sd_i)/k_sd)
    V_block=-32.0, k_block=2.0, tau_h_ms=10.0, theta_block=100.0,  # v_th_eff = v_th + theta*(1-h)
    k_adapt=0.0,                # optional spike-frequency adaptation: v_rest shift -k_adapt*rate (mV/Hz); 0 = Shiu
    V_span=50.0,                # mV; total metabolic+K depolarization soft-capped (tanh) at v_rest+V_span (~ -2 mV plateau)
    # --- injury and death
    k_rep=0.02,                 # 1/s repair of injury integral when energy restored and not depolarized
    J50_s=60.0, J_hill=4.0, h_max=1/30.0,    # hazard(J) = h_max*J^m/(J^m+J50^m)  [1/s]
    # --- reversible synaptic failure (presynaptic release fraction vs energy)
    e_syn_lo=0.3, e_syn_hi=0.6,
    coma_arrest=1.0,            # demand multiplier in depolarization-blocked (comatose) neurons (1 = none; fly-like 0.03-0.2)
    q10_injury=1.0,             # temperature factor on injury accumulation (1 = temperature-independent)
    k_exc=0.0,                  # MECH variant: excitotoxic injury rate (J-units per s) = k_exc*max(r-r0,0)/r0 * Kex/(Kex+K0)
    r0_exc=50.0,                # Hz; firing-rate threshold of the excitotoxic term
)

# ----------------------------------------------------------------------------------------------
# Model loading
# ----------------------------------------------------------------------------------------------
@dataclass
class Model:
    N: int
    indptr: np.ndarray          # out-CSR (rows = presynaptic neuron)
    indices: np.ndarray         # int32 postsynaptic index
    weight: np.ndarray          # float32 signed synapse count (sign * syn_count)
    neurons: pd.DataFrame       # neuron table (idx order)
    sign: np.ndarray            # int8 presynaptic sign
    neuropils: list             # compartment names
    Mn: sp.csr_matrix           # N x P row-normalized neuron->neuropil synapse fractions
    comp_n: np.ndarray          # effective number of neurons per compartment (column sums of Mn)
    adj: np.ndarray             # P x P K+ diffusion adjacency (symmetric, row-max normalized)
    adj_source: str
    MnT: sp.csr_matrix = None   # P x N transpose (cached)
    w_syn: float = 0.275        # mV per synapse used by run() unless lif={'w_syn':...}
    topk: tuple = None          # (cidx (N,3) int64, cw (N,3) float32): top-3 neuropils per neuron
    info: dict = field(default_factory=dict)

    def idx(self, root_ids):
        """root_id(s) -> int indices (missing ids dropped, reported in model.info['last_missing'])."""
        rid = np.atleast_1d(np.asarray(root_ids, dtype=np.int64))
        ids = self.neurons.root_id.values
        p = np.minimum(np.searchsorted(ids, rid), self.N - 1)
        ok = ids[p] == rid
        self.info["last_missing"] = rid[~ok].tolist()
        return p[ok].astype(np.int64)

    def select(self, primary_type=None, side=None, **kw):
        """indices by annotation, e.g. select('LC4'), select('DNa02', side='left'), select(sub_class='sugar/water')."""
        m = np.ones(self.N, bool)
        t = self.neurons
        if primary_type is not None:
            m &= t.primary_type.isin(np.atleast_1d(primary_type)).values
        if side is not None:
            m &= (t.side == side).values
        for k, v in kw.items():
            m &= t[k].isin(np.atleast_1d(v)).values
        return np.flatnonzero(m)


def _sign_vector(nt, glut_inhibitory=True, unknown_sign=1, modulatory_sign=1):
    nt = pd.Series(nt).fillna("UNKNOWN").str.upper().values
    s = np.ones(len(nt), np.int8)
    s[nt == "GABA"] = -1
    if glut_inhibitory:
        s[nt == "GLUT"] = -1
    s[np.isin(nt, ["DA", "SER", "OCT"])] = modulatory_sign
    s[nt == "UNKNOWN"] = unknown_sign
    return s


def neuropil_adjacency_shared(Mraw, min_frac=0.0):
    """Fallback adjacency: number of neurons with synapses in both neuropils (BtB), symmetric, row-max=1."""
    B = (Mraw > 0).astype(np.float32)
    A = (B.T @ B).toarray().astype(np.float64)
    np.fill_diagonal(A, 0)
    A = A / A.max(axis=1, keepdims=True).clip(1e-12)
    A = 0.5 * (A + A.T)
    A[A < min_frac] = 0
    return A


def _topk_compartments(Mraw, k=3):
    D = Mraw.toarray().astype(np.float32)
    ci = np.argsort(-D, axis=1)[:, :k]
    w = np.take_along_axis(D, ci, 1)
    s = w.sum(1, keepdims=True)
    w = np.where(s > 0, w / np.where(s > 0, s, 1), 0).astype(np.float32)
    return ci.astype(np.int64), w


def load_model(connectome="princeton5", glut_inhibitory=True, unknown_sign=1, modulatory_sign=1,
               adjacency="auto", cache=CACHE, w_syn=None):
    """Load the FAFB v783 connectome.

    connectome : 'princeton5'  -> connections_princeton (2).csv.gz, FlyWire v783 Princeton synapses,
                                  pairs with >=5 synapses (3,732,460 edges); sign from neuron NT
                                  (neurons (2).csv.gz nt_type, fallback to edge-majority nt_type).
                 'shiu783'     -> Shiu et al. Connectivity_783.parquet (15,091,983 unthresholded pairs,
                                  Shiu's own 'Excitatory' sign), for exact replication of the paper.
    glut_inhibitory : Shiu convention (True). False -> glutamate excitatory (princeton5 only).
    unknown_sign / modulatory_sign : sign for neurons with no NT call / DA,SER,OCT (Shiu: +1).
    adjacency : 'auto' (geometric cache if present, else shared-neuron), 'shared', 'geometric'.
    w_syn : mV per synapse. Default: 0.275 for shiu783 (paper value); 0.185 for princeton5 — calibrated so the
            network is subcritical like Shiu's (with 0.275 the Princeton >=5 matrix reverberates indefinitely
            after a 150-ms sugar-GRN stimulus; at 0.185 activity stops within 50 ms of stimulus offset and
            MN9 = 107 Hz vs 97 Hz on shiu783 at 150 Hz). Pass w_syn=0.275 for literal Shiu parameters.
    """
    t0 = time.time()
    nt_tab = pd.read_csv(os.path.join(cache, "fafb783_neurons.csv.gz"))
    N = len(nt_tab)
    if connectome == "princeton5":
        A = sp.load_npz(os.path.join(cache, "fafb783_syncount_pre_post_csr.npz")).tocsr()
        sign = _sign_vector(nt_tab.nt.values, glut_inhibitory, unknown_sign, modulatory_sign)
        thr = ">=5 synapses per (pre,post) pair summed over neuropils (Princeton v783 thresholded release)"
    elif connectome == "shiu783":
        A = sp.load_npz(os.path.join(cache, "shiu783_syncount_pre_post_csr.npz")).tocsr()
        sign = np.load(os.path.join(cache, "shiu783_presyn_sign.npy")).astype(np.int8)
        sign[sign == 0] = 1
        thr = "none (Shiu Connectivity_783.parquet, all pairs, >=1 synapse; different synapse set)"
    else:
        raise ValueError(connectome)
    A.sort_indices()
    deg = np.diff(A.indptr)
    w = (A.data.astype(np.float32) * np.repeat(sign, deg).astype(np.float32))
    Mraw = sp.load_npz(os.path.join(cache, "fafb783_neuron_neuropil_syn.npz")).tocsr().astype(np.float32)
    rs = np.asarray(Mraw.sum(1)).ravel()
    Mn = sp.diags((1.0 / np.where(rs > 0, rs, 1)).astype(np.float32)) @ Mraw
    Mn = Mn.tocsr().astype(np.float32)
    nps = pd.read_csv(os.path.join(cache, "fafb783_neuropils.csv")).neuropil.tolist()
    comp_n = np.asarray(Mn.sum(0)).ravel().astype(np.float64)
    geo = os.path.join(cache, "fafb783_neuropil_adjacency_geometric.npy")
    if adjacency == "geometric" or (adjacency == "auto" and os.path.exists(geo)):
        adj = np.load(geo); src = "geometric (" + geo + ")"
    else:
        adj = neuropil_adjacency_shared(Mraw); src = "shared-neuron counts (B^T B, row-max normalized)"
    topk = _topk_compartments(Mraw, 3)
    if w_syn is None:
        w_syn = 0.185 if connectome == "princeton5" else 0.275
    m = Model(N=N, topk=topk, w_syn=float(w_syn), MnT=Mn.T.tocsr().astype(np.float32), indptr=A.indptr.astype(np.int64), indices=A.indices.astype(np.int32), weight=w,
              neurons=nt_tab, sign=sign, neuropils=nps, Mn=Mn, comp_n=comp_n, adj=adj, adj_source=src,
              info=dict(connectome=connectome, edge_threshold=thr, n_edges=int(A.nnz),
                        total_synapses=int(round(float(A.data.astype(np.float64).sum()))),
                        glut_inhibitory=glut_inhibitory, w_syn=float(w_syn), load_s=round(time.time() - t0, 1)))
    return m

# ----------------------------------------------------------------------------------------------
# Stimulus / schedule helpers
# ----------------------------------------------------------------------------------------------
def piecewise(points):
    """[(t_s, value), ...] -> callable t_s -> value (piecewise constant, value may be scalar or array)."""
    ts = [p[0] for p in points]; vs = [p[1] for p in points]
    def f(t):
        i = np.searchsorted(ts, t, side="right") - 1
        return vs[max(i, 0)]
    return f


def ramp(points):
    """[(t_s, value), ...] -> piecewise-linear callable (scalars only)."""
    ts = np.array([p[0] for p in points], float); vs = np.array([p[1] for p in points], float)
    return lambda t: float(np.interp(t, ts, vs))


def xenova_pulse(strength=180.0, hold_ms=650.0, decay_ms=220.0, duration_ms=2500.0, t_on_ms=0.0):
    """Rate envelope of the Xenova fruit-fly-simulation PulseBank ('turn' profile by default;
    'paint' profile: hold 80, decay 180, duration 1000). Returns callable t_ms -> Hz."""
    def f(t_ms):
        ms = t_ms - t_on_ms
        if ms < 0 or ms >= duration_ms:
            return 0.0
        r = strength * math.exp(-max(0.0, ms - hold_ms) / decay_ms)
        return r if r >= 1 else 0.0
    return f


def _as_rate_fn(r):
    return r if callable(r) else (lambda t_ms, _r=float(r): _r)

# ----------------------------------------------------------------------------------------------
# Simulation state / result
# ----------------------------------------------------------------------------------------------
@dataclass
class State:
    t_ms: float
    v: np.ndarray; g: np.ndarray
    ref_until: np.ndarray       # step index until which neuron is refractory
    ring: list                  # delayed spike lists
    step: int
    # slow
    a: np.ndarray = None; a0: np.ndarray = None; r: np.ndarray = None; J: np.ndarray = None
    h: np.ndarray = None; alive: np.ndarray = None; K: np.ndarray = None; R: np.ndarray = None
    a_glia: np.ndarray = None; V_sd: np.ndarray = None


class Result(dict):
    """dict with keys: counts (n_epochs x N int32), epochs, pop_ms (n_ms x n_groups int32), group_names,
    raster_t_ms / raster_idx, slow (dict of time series), state (final State), wall_s, params."""
    __getattr__ = dict.__getitem__


def init_state(model, metab=None, seed=0):
    N = model.N
    p = dict(LIF_DEFAULTS)
    D = int(round(p["t_delay"] / p["dt"]))
    st = State(t_ms=0.0, v=np.full(N, p["v_rest"], np.float32), g=np.zeros(N, np.float32),
               ref_until=np.zeros(N, np.int64), ring=[np.empty(0, np.int64) for _ in range(D)], step=0)
    if metab is not None:
        mp = dict(METAB_DEFAULTS); mp.update(metab)
        a0 = np.full(N, mp["a0"])
        rng = np.random.default_rng(seed + 7)
        a0 = a0.astype(np.float32)
        st.a = a0.copy(); st.a0 = a0; st.r = np.zeros(N, np.float32); st.J = np.zeros(N, np.float32); st.h = np.ones(N, np.float32)
        st.alive = np.ones(N, bool); st.K = np.full(len(model.neuropils), mp["K0"]); st.R = np.ones(len(model.neuropils))
        st.a_glia = np.full(len(model.neuropils), mp["a0"])
        st.V_sd = (mp["V_sd"] + mp["sd_jitter"] * rng.standard_normal(N)).astype(np.float32)
    return st


def demand_base(model, mp):
    t = model.neurons
    mode = mp["demand_mode"]
    if mode == "uniform":
        x = np.ones(model.N)
    elif mode == "area":
        x = t.area_nm.values.astype(float)
    elif mode == "length":
        x = t.length_nm.values.astype(float)
    elif mode == "surface_to_volume":
        x = t.area_nm.values.astype(float) / np.maximum(t.size_nm.values.astype(float), 1)
    elif mode in t.columns:
        x = t[mode].values.astype(float)
    else:
        raise ValueError(mode)
    med = np.nanmedian(x[x > 0])
    x = np.where(np.isfinite(x) & (x > 0), x, med)
    return np.clip((x / med) ** mp["demand_beta"], mp["base_min"], mp["base_max"])

# ----------------------------------------------------------------------------------------------
# Core engine
# ----------------------------------------------------------------------------------------------
def _gather(indptr, s):
    st = indptr[s]
    ln = indptr[s + 1] - st
    tot = int(ln.sum())
    if tot == 0:
        return None, ln
    cs = np.cumsum(ln)
    offs = np.arange(tot, dtype=np.int64) + np.repeat(st - (cs - ln), ln)
    return offs, ln


def _sig(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -50, 50)))


def run(model, duration_ms, stim=None, metab=None, supply=1.0, temperature=22.0, interventions=None,
        state=None, epochs=None, groups="super_class", record_idx=None, max_raster=2_000_000,
        rec_every_ms=10.0, seed=0, lif=None, verbose=False, inh_gain=1.0, exc_gain=1.0):
    """Simulate `duration_ms` of biological time.

    stim          : list of dicts {'idx': int array, 'rate': Hz or callable(t_ms)->Hz,
                    't_on': ms, 't_off': ms}  -> independent Poisson input (Shiu: +f_poi*w_syn to v).
    metab         : None (pure Shiu LIF) or dict overriding METAB_DEFAULTS (enables slow layer).
    supply        : O2/glucose availability S in [0,1]: scalar, (P,) per-neuropil array, or callable
                    t_s -> scalar | (P,) array (see piecewise / ramp).
    temperature   : deg C, scalar or callable t_s -> deg C.
    interventions : dict of optional knobs (scalars or callables of t_s):
                    'k_rel_scale' (SD-release blocker e.g. NMDA/ketamine analog), 'clr_scale' (K+ buffering),
                    'injury_scale' (neuroprotection), 'atp_floor' (min ATP fraction), 'v_offset' (mV, scalar
                    or (N,) array; hyperpolarization/depolarization), 'K_clamp' (mM; clamp all K_c),
                    'silence' (int array: neurons whose output is removed), 'revive' (bool: dead neurons stay
                    dead unless True -> resurrect with J=0 at run start — counterfactual only).
    state         : State to continue from (returned in result['state']); otherwise a fresh baseline.
    epochs        : list of (t0_ms, t1_ms) relative to run start; per-neuron spike counts per epoch.
    groups        : neuron-table column used for population rate per ms, or None.
    record_idx    : None | 'all' | int array — neurons whose spikes are stored in the raster (capped).
    """
    P = dict(LIF_DEFAULTS); P["w_syn"] = getattr(model, "w_syn", P["w_syn"]); P.update(lif or {})
    dt = P["dt"]; N = model.N
    nsteps = int(round(duration_ms / dt))
    D = int(round(P["t_delay"] / dt)); NREF = int(round(P["t_ref"] / dt))
    em = np.float32(math.exp(-dt / P["tau_m"])); eg = np.float32(math.exp(-dt / P["tau_syn"]))
    cpl = np.float32(P["tau_syn"] / (P["tau_m"] - P["tau_syn"]) * (math.exp(-dt / P["tau_m"]) - math.exp(-dt / P["tau_syn"])))
    w_scale = np.float32(P["w_syn"])
    poi_w = np.float32(P["f_poi"] * P["w_syn"])
    rng = np.random.default_rng(seed)
    iv = dict(interventions or {})
    ivf = lambda k, t, d=None: (iv[k](t) if callable(iv[k]) else iv[k]) if k in iv else d

    use_metab = metab is not None
    if state is None:
        state = init_state(model, metab, seed)
    elif use_metab and state.a is None:
        st2 = init_state(model, metab, seed); st2.v, st2.g, st2.ref_until, st2.ring, st2.step = state.v, state.g, state.ref_until, state.ring, state.step
        state = st2
    v, g, ref_until, ring = state.v, state.g, state.ref_until, state.ring
    step0 = state.step
    indptr, indices, weight = model.indptr, model.indices, model.weight * w_scale
    if inh_gain != 1.0 or exc_gain != 1.0:
        weight = np.where(weight < 0, weight * np.float32(inh_gain), weight * np.float32(exc_gain)).astype(np.float32)

    # stimulation
    stim = stim or []
    stim_idx = [np.asarray(s["idx"], np.int64) for s in stim]
    stim_rate = [_as_rate_fn(s.get("rate", 150.0)) for s in stim]
    stim_on = [s.get("t_on", 0.0) for s in stim]; stim_off = [s.get("t_off", 1e18) for s in stim]
    poisson_target = np.zeros(N, bool)
    for ix in stim_idx:
        poisson_target[ix] = True     # Shiu: refractory set to 0 for Poisson targets

    silence = np.zeros(N, bool)
    if "silence" in iv and iv["silence"] is not None:
        silence[np.asarray(iv["silence"], np.int64)] = True
    out_scale = None if not silence.any() else (~silence).astype(np.float32)

    # recording
    epochs = epochs or [(0.0, duration_ms)]
    ep_steps = [(int(round(a / dt)), int(round(b / dt))) for a, b in epochs]
    counts = np.zeros((len(epochs), N), np.int32)
    n_ms = int(math.ceil(duration_ms))
    if groups is not None:
        glab, gnames = pd.factorize(model.neurons[groups].fillna("NA"))
        glab = glab.astype(np.int64); ng = len(gnames)
    else:
        glab = np.zeros(N, np.int64); gnames = ["all"]; ng = 1
    pop = np.zeros((n_ms, ng), np.int32)
    rec_mask = None
    if record_idx is not None:
        rec_mask = np.ones(N, bool) if (isinstance(record_idx, str) and record_idx == "all") else np.zeros(N, bool)
        if not isinstance(record_idx, str):
            rec_mask[np.asarray(record_idx, np.int64)] = True
    rast_t, rast_i, n_rast = [], [], 0

    # metabolic setup
    if use_metab:
        mp = dict(METAB_DEFAULTS); mp.update(metab)
        Kslow = int(round(mp["dt_slow_ms"] / dt))
        dts = mp["dt_slow_ms"] * 1e-3 * mp["speedup"]           # slow-time step in (compressed) seconds
        base = demand_base(model, mp).astype(np.float32)
        Mn32 = model.Mn; MnT32 = model.MnT
        has_comp = np.asarray(Mn32.sum(1)).ravel() > 0
        comp_n = np.maximum(model.comp_n, 1e-9)
        cidx, cw = model.topk
        cwf = cw.ravel(); cif = cidx.ravel()
        def mix(xc, fill):
            y = (cw * xc[cidx]).sum(1, dtype=np.float32)
            if fill: y[~has_comp] = fill
            return y
        def cmean(x):
            return np.bincount(cif, weights=(cw * x[:, None]).ravel(), minlength=len(model.neuropils)) / comp_n
        dep = np.zeros(N, np.float32); hinf = np.ones(N, np.float32)
        Adj = model.adj; Lap = Adj - np.diag(Adj.sum(1))
        Pn = len(model.neuropils)
        a, a0, r, J, h, alive, K, R, ag, Vsd = (state.a, state.a0, state.r, state.J, state.h, state.alive,
                                                 state.K, state.R, state.a_glia, state.V_sd)
        inv_a0 = (1.0 / a0).astype(np.float32)
        if getattr(state, "Jx", None) is None:
            state.Jx = np.zeros(N, np.float32)
        Jx = state.Jx
        Kadp = np.float32((1 - mp["a0"]) * (mp["reserve"] / mp["a0"] - 1.0))
        if ivf("revive", 0.0, False):
            alive[:] = True; J[:] = 0
        vr = np.full(N, P["v_rest"], np.float32); vreset = np.full(N, P["v_reset"], np.float32)
        vth = np.full(N, P["v_th"], np.float32)
        frel = np.ones(N, np.float32)
        spk_acc = np.zeros(N, np.float32)
        slow_rec = {k: [] for k in ["t_s", "e_mean", "e_lt05", "K", "dep_frac_np", "e_np", "dep_frac", "block_frac",
                                     "dead_frac", "vr_mean", "syn_frac", "J_mean", "S_mean", "T", "Jx_mean"]}
        rec_every = max(1, int(round(rec_every_ms / mp["dt_slow_ms"])))
        nslow = 0
        a_floor_default = 0.0
        rng_death = np.random.default_rng(seed + 99991)

        def slow_update(t_ms_abs):
            nonlocal nslow
            ts = t_ms_abs * 1e-3 * mp["speedup"]
            S = supply(ts) if callable(supply) else supply
            S = np.broadcast_to(np.asarray(S, np.float32), (Pn,)).astype(np.float32) if np.ndim(S) else np.full(Pn, float(S), np.float32)
            T = temperature(ts) if callable(temperature) else temperature
            Qd = np.float32(mp["q10_demand"] ** ((T - mp["T0"]) / 10.0)); Qp = mp["q10_pump"] ** ((T - mp["T0"]) / 10.0)
            # rate estimate (Hz)
            r[:] += (spk_acc * np.float32(1.0 / (mp["dt_slow_ms"] * 1e-3)) - r) * np.float32(1 - math.exp(-mp["dt_slow_ms"] * 1e-3 / mp["tau_rate_s"]))
            spk_acc[:] = 0
            Kex = np.maximum(K - mp["K0"], 0).astype(np.float32)
            dVk_c = (mp["g_K"] * mp["RT_F"] * np.log(np.maximum(K, 0.1) / mp["K0"])).astype(np.float32)
            Z = Mn32 @ np.stack([S, Kex, dVk_c], 1).astype(np.float32)   # neuron-level mixtures over all neuropils
            S_i = Z[:, 0].copy(); S_i[~has_comp] = float(S.mean()); Kex_i = Z[:, 1]; dVk_i = Z[:, 2]
            d = base * Qd * (1.0 + r * np.float32(1.0 / mp["r_half_hz"]) + np.float32(mp["k_K"] / mp["K0"]) * Kex_i)
            dsc = ivf("demand_scale", ts, None)
            if dsc is not None:
                d = d * (np.float32(dsc) if np.ndim(dsc) == 0 else np.asarray(dsc, np.float32))
            if mp["coma_arrest"] != 1.0:
                d = d * (1.0 - np.float32(1.0 - mp["coma_arrest"]) * (1.0 - h))
            Qs = np.float32(mp["q10_supply"] ** ((T - mp["T0"]) / 10.0))
            prod = (np.float32(mp["reserve"]) * Qs) * base * S_i * (1 - a) / (Kadp + 1 - a)
            a[:] += np.float32(dts / mp["tau_atp_s"]) * (prod - d * a)
            af = ivf("atp_floor", ts, 0.0)
            np.clip(a, 0.0, 1.0, out=a)
            if af:
                np.maximum(a, (np.float32(af) if np.ndim(af) == 0 else np.asarray(af, np.float32)) * a0, out=a)
            a[~alive] = 0.0
            ag[:] += dts * (mp["reserve"] * Qs * S * (1 - ag) / (Kadp + 1 - ag) - Qd * ag) / mp["tau_atp_s"]
            e = a * inv_a0
            eg_ = ag / mp["a0"]
            # resting potential drift: pump failure (dead-zone e_ok) + extracellular K (per-compartment Nernst-like term mixed)
            dVa = np.clip((mp["e_ok"] - e) * np.float32(1.0 / mp["e_ok"]), 0, 1)
            dVa **= np.float32(mp["n_atp"]); dVa *= np.float32(mp["dV_atp_max"])
            vo = ivf("v_offset", ts, 0.0)
            span = np.float32(mp["V_span"])
            dtot = span * np.tanh((dVa + dVk_i) / span)
            if mp["k_adapt"]:
                dtot = dtot - np.float32(mp["k_adapt"]) * r          # optional spike-frequency adaptation (mV per Hz)        # soft cap: depolarization saturates at v_rest+V_span
            vr_eff = np.float32(P["v_rest"]) + dtot + (np.float32(vo) if np.ndim(vo) == 0 else np.asarray(vo, np.float32))
            # spreading depolarization indicator and depolarization block, computed only near threshold
            dep[:] = 0.0
            cand = np.flatnonzero(vr_eff > (Vsd - 8 * mp["k_sd"]))
            if cand.size:
                x = _sig((vr_eff[cand] - Vsd[cand]) / mp["k_sd"]).astype(np.float32)
                x[x < 1e-3] = 0.0
                dep[cand] = x * alive[cand]
            hinf.fill(1.0)
            cb = np.flatnonzero(vr_eff > (mp["V_block"] - 8 * mp["k_block"]))
            if cb.size:
                y = (1.0 - _sig((vr_eff[cb] - mp["V_block"]) / mp["k_block"])).astype(np.float32)
                y[y > 1 - 1e-3] = 1.0
                hinf[cb] = y
            if cb.size or (h < 1).any():
                h[:] += (hinf - h) * np.float32(1 - math.exp(-mp["dt_slow_ms"] / mp["tau_h_ms"]))
                h[h > 1 - 1e-4] = 1.0
            inact = 1.0 - h
            # compartment K+
            ksc = ivf("k_rel_scale", ts, 1.0); csc = ivf("clr_scale", ts, 1.0)
            XT = (MnT32 @ np.stack([np.clip(e, 0, 1.2) * alive, dep, r * alive], 1).astype(np.float32)) / comp_n[:, None]
            e_np, dep_np, r_np = XT[:, 0].astype(np.float64), XT[:, 1].astype(np.float64), XT[:, 2].astype(np.float64)
            pump = Qp * (mp["glia_w"] * np.clip(eg_, 0, 1.2) + (1 - mp["glia_w"]) * e_np)
            rel = ksc * mp["k_rel"] * dep_np * R * np.clip(1 - K / mp["K_max"], 0, 1)
            dK = (mp["k_spk"] * r_np + rel + mp["k_leak"] * np.clip((mp["pump_ok"] - pump) / mp["pump_ok"], 0, 1)
                  - csc * mp["k_clr"] * pump * (K - mp["K0"]) + mp["k_diff"] * (Lap @ K))
            K[:] = np.clip(K + dts * dK, 0.5, mp["K_max"])
            R[:] += dts * (-dep_np * R / mp["tau_rel_depl_s"] + np.clip(pump, 0, 1) * (1 - R) * (1 - dep_np) / mp["tau_rel_rec_s"])
            np.clip(R, 0, 1, out=R)
            kc = ivf("K_clamp", ts, None)
            if kc is not None:
                K[:] = kc
            # injury integral (Ca2+-load proxy) and stochastic death
            isc = ivf("injury_scale", ts, 1.0)
            kx = float(mp.get("k_exc", 0.0))
            exc = None
            if kx > 0:
                r0x = np.float32(mp["r0_exc"])
                exc = (np.maximum(r - r0x, 0) * (np.float32(1.0) / r0x) * (Kex_i / (Kex_i + np.float32(mp["K0"])))).astype(np.float32)
                exc[~alive] = 0.0
            if cand.size or (J > 0).any() or (exc is not None and (exc > 0).any()):
                cond = (dep > 0) | (J > 0)
                if exc is not None:
                    cond |= exc > 0
                jj = np.flatnonzero(cond)
                ej = np.clip(e[jj], 0, 1)
                qi = np.float32(mp["q10_injury"] ** ((T - mp["T0"]) / 10.0))
                iscj = np.float32(isc) if np.ndim(isc) == 0 else np.asarray(isc, np.float32)[jj]
                add = qi * iscj * dep[jj] * (1 - ej)
                if exc is not None:
                    xadd = (qi * iscj * np.float32(kx)) * exc[jj]
                    add = add + xadd
                    Jx[jj] += np.float32(dts) * xadd
                J[jj] += np.float32(dts) * (add - np.float32(mp["k_rep"]) * ej * (1 - dep[jj]) * J[jj])
                J[jj] = np.maximum(J[jj], 0)
                Jm = J[jj].astype(np.float64) ** mp["J_hill"]
                haz = mp["h_max"] * Jm / (Jm + mp["J50_s"] ** mp["J_hill"])
                if haz.size and haz.max() > 0:
                    die = jj[alive[jj] & (rng_death.random(jj.size) < -np.expm1(-haz * dts))]
                    if die.size:
                        alive[die] = False
            # push to fast layer
            vr[:] = vr_eff; vreset[:] = vr_eff + np.float32(P["v_reset"] - P["v_rest"])
            vth[:] = np.float32(P["v_th"]) + np.float32(mp["theta_block"]) * inact
            vth[~alive] = 1e9
            bm = ivf("block_mask", ts, None)
            if bm is not None:
                vth[bm] = 1e9; g[bm] = 0.0
            if ivf("block_all", ts, False):          # TTX / isoelectric-coma analog: no spikes, synaptic drive cleared
                vth[:] = 1e9; g[:] = 0.0
            np.clip((e - mp["e_syn_lo"]) * np.float32(1.0 / (mp["e_syn_hi"] - mp["e_syn_lo"])), 0, 1, out=frel)
            frel[~alive] = 0.0
            if nslow % rec_every == 0:
                sr = slow_rec
                sr["t_s"].append(ts); sr["e_mean"].append(float(e[alive].mean()) if alive.any() else 0.0)
                sr["e_lt05"].append(float((e < 0.5).mean()))
                sr["K"].append(K.astype(np.float32).copy()); sr["dep_frac_np"].append(dep_np.astype(np.float32))
                sr["e_np"].append(e_np.astype(np.float32))
                sr["dep_frac"].append(float((dep > 0.5).mean())); sr["block_frac"].append(float((inact > 0.5).mean()))
                sr["dead_frac"].append(float(1 - alive.mean())); sr["vr_mean"].append(float(vr_eff.mean()))
                sr["syn_frac"].append(float(frel.mean())); sr["J_mean"].append(float(J.mean()))
                sr["S_mean"].append(float(S.mean())); sr["T"].append(float(T)); sr["Jx_mean"].append(float(Jx.mean()))
            nslow += 1

        slow_update(state.t_ms)
    else:
        vr = np.float32(P["v_rest"]); vreset = np.float32(P["v_reset"]); vth = np.float32(P["v_th"]); frel = None

    tmp = np.empty(N, np.float32)
    tol = float(P.get("active_tol", 1e-5))
    act = np.zeros(N, bool)                    # active set: neurons not exactly at rest
    vr_arr = vr if use_metab else None
    v_r0 = (np.abs(v - vr) > 0) | (g != 0)
    act |= v_r0
    if use_metab:
        act |= vr > (vth - 1.0)
        vr_prev = vr.copy()
    act_idx = np.flatnonzero(act)
    pending = []
    ref_ring = [np.empty(0, np.int64) for _ in range(NREF)]
    dense_frac = 0.3
    act_sum = 0
    t_wall = time.time()
    for k in range(nsteps):
        s_abs = step0 + k
        t_ms = k * dt
        if pending:
            act_idx = np.concatenate([act_idx] + pending); pending = []
        dense = act_idx.size > dense_frac * N
        act_sum += act_idx.size
        # 1) state update (exact linear integration), refractory neurons frozen
        ref_idx = np.concatenate(ref_ring) if NREF else np.empty(0, np.int64)
        if ref_idx.size:
            v_keep = v[ref_idx]; g_keep = g[ref_idx]
        if dense:
            np.subtract(v, vr, out=v); v *= em; v += vr
            np.multiply(g, cpl, out=tmp); v += tmp; g *= eg
        elif act_idx.size:
            vv = v[act_idx]; gg = g[act_idx]
            vrr = vr_arr[act_idx] if use_metab else vr
            vv = vrr + (vv - vrr) * em + gg * cpl
            v[act_idx] = vv; g[act_idx] = gg * eg
        if ref_idx.size:
            v[ref_idx] = v_keep; g[ref_idx] = g_keep
        # 2) threshold (only active neurons can cross)
        if dense:
            spk = np.flatnonzero(v > vth)
        elif act_idx.size:
            vt = vth[act_idx] if use_metab else vth
            spk = act_idx[v[act_idx] > vt]
        else:
            spk = act_idx[:0]
        if spk.size:
            spk = spk[ref_until[spk] <= s_abs]
        # 3) synaptic delivery of spikes emitted D steps ago
        slot = s_abs % D
        arr = ring[slot]
        if arr.size:
            offs, ln = _gather(indptr, arr)
            if offs is not None:
                post = indices[offs]; w = weight[offs]
                if frel is not None or out_scale is not None:
                    f = np.ones(arr.size, np.float32)
                    if frel is not None: f = f * frel[arr]
                    if out_scale is not None: f = f * out_scale[arr]
                    w = w * np.repeat(f, ln)
                if offs.size < 4000:
                    np.add.at(g, post, w)
                    up = np.unique(post)
                else:
                    gi = np.bincount(post, weights=w, minlength=N)
                    up = np.flatnonzero(gi); g[up] += gi[up].astype(np.float32)
                newa = up[~act[up]]
                if newa.size:
                    act[newa] = True; pending.append(newa)
        ring[slot] = spk
        # Poisson input (to v)
        for j, ix in enumerate(stim_idx):
            if stim_on[j] <= t_ms < stim_off[j]:
                rate = stim_rate[j](t_ms)
                if rate > 0:
                    hit = ix[rng.random(ix.size) < rate * dt * 1e-3]
                    if hit.size:
                        v[hit] += poi_w
                        newa = hit[~act[hit]]
                        if newa.size:
                            act[newa] = True; pending.append(newa)
        # 4) reset + refractoriness
        if spk.size:
            if use_metab:
                v[spk] = vreset[spk]
                spk_acc[spk] += 1
            else:
                v[spk] = vreset
            g[spk] = 0.0
            nr = spk[~poisson_target[spk]]
            ref_until[nr] = s_abs + NREF
            ref_ring[s_abs % NREF if NREF else 0] = nr
            for e_i, (a_, b_) in enumerate(ep_steps):
                if a_ <= k < b_:
                    counts[e_i, spk] += 1
            np.add.at(pop[int(t_ms)], glab[spk], 1)
            if rec_mask is not None and n_rast < max_raster:
                sr_ = spk[rec_mask[spk]]
                if sr_.size:
                    rast_t.append(np.full(sr_.size, t_ms, np.float32)); rast_i.append(sr_.astype(np.int32)); n_rast += sr_.size
        elif NREF:
            ref_ring[s_abs % NREF] = np.empty(0, np.int64)
        # 5) slow layer
        if use_metab and (k + 1) % Kslow == 0:
            slow_update(state.t_ms + (k + 1) * dt)
            ch = np.flatnonzero((np.abs(vr - vr_prev) > tol) | (vr > vth - 1.0))
            vr_prev[:] = vr
            newa = ch[~act[ch]]
            if newa.size:
                act[newa] = True; pending.append(newa)
        # prune active set every 10 steps: snap neurons within tol of rest back to exact rest
        if k % 10 == 9 and not dense and act_idx.size:
            if pending:
                act_idx = np.concatenate([act_idx] + pending); pending = []
            vrr = vr_arr[act_idx] if use_metab else vr
            vt = vth[act_idx] if use_metab else vth
            quiet = (np.abs(v[act_idx] - vrr) < tol) & (np.abs(g[act_idx]) < tol) & (ref_until[act_idx] <= s_abs) & (vrr < vt - 1.0)
            if quiet.any():
                q = act_idx[quiet]
                v[q] = vrr[quiet] if use_metab else vr
                g[q] = 0.0; act[q] = False
                act_idx = act_idx[~quiet]
        elif dense and k % 10 == 9:
            act_idx = np.flatnonzero(act)
        if verbose and k % 1000 == 0 and k:
            print(f"t={t_ms:.0f} ms  active={act_idx.size} wall={time.time()-t_wall:.1f}s", flush=True)
    wall = time.time() - t_wall
    state.step = step0 + nsteps; state.t_ms = state.t_ms + duration_ms
    res = Result(counts=counts, epochs=epochs, pop_ms=pop, group_names=list(gnames),
                 raster_t_ms=np.concatenate(rast_t) if rast_t else np.empty(0, np.float32),
                 raster_idx=np.concatenate(rast_i) if rast_i else np.empty(0, np.int32),
                 state=state, wall_s=wall, duration_ms=duration_ms, lif=P,
                 metab=(mp if use_metab else None), mean_active=act_sum / max(nsteps, 1))
    if use_metab:
        sl = {k: (np.array(v_) if k not in ("K", "dep_frac_np", "e_np") else np.vstack(v_)) for k, v_ in slow_rec.items()}
        res["slow"] = sl
        res["alive"] = alive.copy(); res["e"] = (a / a0).astype(np.float32); res["J"] = J.astype(np.float32)
    return res

# ----------------------------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------------------------
def rates(res, idx=None, epoch=0):
    """Mean firing rate (Hz) per neuron in an epoch."""
    a, b = res["epochs"][epoch]
    r = res["counts"][epoch] / ((b - a) * 1e-3)
    return r if idx is None else r[np.asarray(idx)]


def run_trials(model, n_trials, duration_ms, seed0=0, **kw):
    """Repeat `run` with different seeds; returns mean rates (Hz) per neuron (epoch 0) and list of results' walls."""
    acc = None; walls = []
    for i in range(n_trials):
        r_ = run(model, duration_ms, seed=seed0 + i, **kw)
        x = rates(r_); acc = x if acc is None else acc + x; walls.append(r_["wall_s"])
    return acc / n_trials, walls


def functional_identity_index(rate_ref, rate_test, min_rate=1.0, eps=1.0):
    """Functional Identity Index (FII) between a reference (pre-insult) and test (post-resuscitation)
    evoked-response profile over the same stimulus.
    Neurons considered: responsive (>= min_rate Hz) in either condition.
      corr    = Pearson r of log(rate+eps) over considered neurons (pattern similarity)
      recall  = fraction of reference-responsive neurons still responsive (pathway persistence)
      precision = fraction of test-responsive neurons that were responsive in reference (no novel paths)
      FII = corr_clipped * sqrt(recall * precision)  in [0,1]
    """
    rr, rt = np.asarray(rate_ref, float), np.asarray(rate_test, float)
    A, B = rr >= min_rate, rt >= min_rate
    U = A | B
    if U.sum() < 3 or A.sum() == 0:
        return dict(FII=float("nan"), corr=float("nan"), recall=float("nan"), precision=float("nan"), n=int(U.sum()))
    x, y = np.log(rr[U] + eps), np.log(rt[U] + eps)
    c = float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 and y.std() > 0 else 0.0
    rec = float((A & B).sum() / A.sum()); prec = float((A & B).sum() / max(B.sum(), 1))
    return dict(FII=max(c, 0.0) * math.sqrt(rec * prec), corr=c, recall=rec, precision=prec, n=int(U.sum()))


def population_rate(res, group=None, bin_ms=10):
    """Population rate (Hz per neuron) time series for a group name (None = all) with bins of bin_ms."""
    pop = res["pop_ms"]
    x = pop.sum(1) if group is None else pop[:, res["group_names"].index(group)]
    nb = len(x) // bin_ms
    x = x[: nb * bin_ms].reshape(nb, bin_ms).sum(1)
    return x


# Shiu et al. 2024 reference neurons (root ids, FAFB v783)
SHIU_SUGAR_GRN = [720575940624963786, 720575940630233916, 720575940637568838, 720575940638202345,
                  720575940617000768, 720575940630797113, 720575940632889389, 720575940621754367,
                  720575940621502051, 720575940640649691, 720575940639332736, 720575940616885538,
                  720575940639198653, 720575940620900446, 720575940617937543, 720575940632425919,
                  720575940633143833, 720575940612670570, 720575940628853239, 720575940629176663,
                  720575940611875570]
SHIU_MN9 = 720575940660219265
XENOVA_READOUTS = dict(walk=["DNp09", "DNg100", "DNg97"], turn=["DNa02", "DNa11", "DNg13"], reverse=["MDN"], escape=["DNp01"])
