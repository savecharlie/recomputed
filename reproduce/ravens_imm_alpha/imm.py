#!/usr/bin/env python3
"""The individual mobility model of arXiv:2609.01858 Sec. IV A, rebuilt.

Two rules, taken verbatim from the paper:
    Pnew(t) = q * n(t)^(-beta)                                          (1)
    Pi(t)   = mi(t)^alpha / sum_j mj(t)^alpha       (preferential return)(2)

The known answers this instrument must reproduce BEFORE it is allowed to say
anything new (Sec. IV A):
    <n(t)>  = [(1+beta) q t]^(1/(1+beta))            for q t >> 1        (5)
    sigma_n = sqrt(<n>/(1+2 beta))
    <mi(t)> ~ t   at alpha = 1
Written from the equations, not from the authors' code.

One reading decision, stated because it changes the answer: at a step the walker
either discovers a new site (which then has m=1) or returns to a known one. I
count BOTH as a step, so sum_j mj(t) = t exactly -- the paper says so explicitly
("the total number of visits sum_j mj(t) equals the total number of steps t").
"""
import numpy as np


def run(T, alpha, beta, q, rng, track_ranks=3):
    """One trajectory. Returns n(t) and m_i(t) for the top `track_ranks` sites,
    with ranks re-evaluated at each recorded time (the paper ranks by visitation
    intensity, i=1 the most visited)."""
    m = np.zeros(T + 2)          # visit counts, site k in slot k
    n = 1                        # site 0 exists at t=0 (P(n,0)=delta_{n,1})
    m[0] = 1.0
    nt = np.empty(T + 1, dtype=np.int64)
    mt = np.empty((T + 1, track_ranks))
    nt[0] = n
    mt[0] = np.sort(m[:n])[::-1][:track_ranks].tolist() + [0.0] * max(0, track_ranks - n)
    w = m[:1] ** alpha           # kept incrementally below for speed
    for t in range(1, T + 1):
        if rng.random() < q * n ** (-beta):
            m[n] = 1.0
            n += 1
        else:
            wk = m[:n] ** alpha
            j = rng.choice(n, p=wk / wk.sum())
            m[j] += 1.0
        nt[t] = n
        top = np.sort(m[:n])[::-1]
        mt[t, :] = [top[k] if k < n else 0.0 for k in range(track_ranks)]
    return nt, mt


def ensemble(R, T, alpha, beta, q, seed=0, track_ranks=3):
    rng = np.random.default_rng(seed)
    N = np.empty((R, T + 1)); M = np.empty((R, T + 1, track_ranks))
    for r in range(R):
        N[r], M[r] = run(T, alpha, beta, q, rng, track_ranks)
    return N, M
