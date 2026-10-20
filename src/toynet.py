"""Small synthetic connectome for tests and smoke runs (no FlyWire data needed).

make_toy_model(fs) returns a `fs.Model` (fs = flysim, flysim_dyn or flysim_mech) with the same fields that
`load_model` fills for FAFB, so `fs.run` exercises the identical code path: signed out-CSR, a neuron table,
P neuropil compartments with row-normalised synapse fractions, a shared-neuron K+ adjacency and top-3 compartments.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp


def make_toy_model(fs, N=300, P=4, p_conn=0.04, frac_inh=0.3, syn_max=6, w_syn=0.185, seed=0):
    rng = np.random.default_rng(seed)
    A = sp.random(N, N, density=p_conn, format="csr", random_state=rng,
                  data_rvs=lambda k: rng.integers(1, syn_max + 1, k)).astype(np.float32)
    A.setdiag(0); A.eliminate_zeros(); A.sort_indices()
    nt = np.where(rng.random(N) < frac_inh, "GABA", "ACH")
    sign = fs._sign_vector(nt)
    deg = np.diff(A.indptr)
    w = (A.data * np.repeat(sign, deg)).astype(np.float32)
    comp = rng.integers(0, P, N)
    Mraw = sp.csr_matrix((rng.integers(5, 50, N).astype(np.float32), (np.arange(N), comp)), shape=(N, P))
    Mraw = (Mraw + sp.random(N, P, density=0.3, format="csr", random_state=rng, data_rvs=lambda k: rng.integers(1, 5, k)).astype(np.float32)).tocsr()
    rs = np.asarray(Mraw.sum(1)).ravel()
    Mn = (sp.diags((1.0 / np.where(rs > 0, rs, 1)).astype(np.float32)) @ Mraw).tocsr().astype(np.float32)
    neurons = pd.DataFrame(dict(idx=np.arange(N), root_id=np.arange(N, dtype=np.int64) + 10**17, nt=nt,
                                super_class=np.where(np.arange(N) < 20, "sensory", "central"),
                                primary_type=[f"T{i % 10}" for i in range(N)], side="left",
                                area_nm=rng.lognormal(0, 0.5, N) * 1e8, length_nm=rng.lognormal(0, 0.5, N) * 1e6,
                                size_nm=rng.lognormal(0, 0.5, N) * 1e9))
    return fs.Model(N=N, indptr=A.indptr.astype(np.int64), indices=A.indices.astype(np.int32), weight=w,
                    neurons=neurons, sign=sign, neuropils=[f"NP{i}" for i in range(P)], Mn=Mn,
                    comp_n=np.asarray(Mn.sum(0)).ravel().astype(np.float64),
                    adj=fs.neuropil_adjacency_shared(Mraw), adj_source="toy shared-neuron",
                    MnT=Mn.T.tocsr().astype(np.float32), w_syn=float(w_syn), topk=fs._topk_compartments(Mraw, 3),
                    info=dict(connectome="toy", seed=seed))


SENSORY = np.arange(20)
