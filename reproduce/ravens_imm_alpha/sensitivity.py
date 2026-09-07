"""Two things must be true before the alpha interval means anything.

(1) SUSPECT THE RULER. n(200 d) was read off a log-log figure by eye -- 57 for
    FT, 190 for BN. q is fitted to that, and q shifts how many sites compete for
    the returns, which is exactly what sets p1. So: how far does the inferred
    alpha move if my reading of the figure is 25% wrong in either direction?
    And if their beta were 0.1 off?

(2) WHY THE EYE-GUIDE COULD NOT SEE IT. The paper's alpha diagnostic is the
    log-log slope of m_i(t) over the plotted window, compared to a dashed line
    of slope 1. Measure that slope as a function of alpha at the ravens' actual
    run length. If it stays near 1 over a wide band of alpha, the diagnostic
    was under-powered and p1 is simply the better statistic -- which is the
    constructive half of this note.
"""
import numpy as np
from multiprocessing import Pool
from alpha_scan import run, fit_q, SPECIES, DAYS

R = 40


def p1_of(alpha, beta, q, T, seed):
    rng = np.random.default_rng(seed)
    return float(np.mean([run(T, alpha, beta, q, rng, ranks=1)[1][0] for _ in range(R)]))


def invert(beta, q, T, target, seed=77):
    """alpha such that simulated p1 == target, by bisection on a monotone curve."""
    lo, hi = 0.4, 1.8
    for _ in range(14):
        a = 0.5 * (lo + hi)
        lo, hi = (a, hi) if p1_of(a, beta, q, T, seed) < target else (lo, a)
    return 0.5 * (lo + hi)


def cell(args):
    sp, dn, db = args
    d = SPECIES[sp]
    T = int(round(DAYS * d["lam"]))
    beta = d["beta"] + db
    n200 = d["n200"] * (1 + dn)
    q = fit_q(beta, T, n200)
    return sp, dn, db, q, invert(beta, q, T, d["p1"])


if __name__ == "__main__":
    tasks = [(sp, dn, db) for sp in SPECIES for dn in (-0.25, 0.0, 0.25) for db in (-0.10, 0.0, 0.10)]
    with Pool(5) as pool:
        res = pool.map(cell, tasks)
    print("alpha inferred from the paper's own p1, under perturbations of my figure reading")
    print("  sp   n(200d)   beta     q      alpha")
    for sp, dn, db, q, a in res:
        d = SPECIES[sp]
        print(f"  {sp}   {d['n200']*(1+dn):6.1f}   {d['beta']+db:4.2f}  {q:5.3f}   {a:5.3f}")
    for sp in SPECIES:
        al = [a for s, _, _, _, a in res if s == sp]
        print(f"  -> {sp}: alpha in {min(al):.2f} .. {max(al):.2f} across all nine perturbations")

    print("\nlog-log slope of m_i(t) over the plotted window (days 10-200), vs alpha")
    print("  the paper compares this to a dashed line of slope 1")
    for sp, d in SPECIES.items():
        T = int(round(DAYS * d["lam"]))
        q = fit_q(d["beta"], T, d["n200"])
        lo = int(round(10 * d["lam"]))
        cps = sorted(set(np.unique(np.round(np.logspace(np.log10(lo), np.log10(T), 24))).astype(int)))
        print(f"  {sp} (T={T}):")
        for alpha in (0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4):
            rng = np.random.default_rng(5)
            acc = np.zeros((len(cps), 3))
            for _ in range(R):
                _, _, out = run(T, alpha, d["beta"], q, rng, ranks=3, checks=cps)
                acc += np.array([out[c][:3] * c for c in cps])
            acc /= R
            sl = [np.polyfit(np.log(cps), np.log(np.maximum(acc[:, i], 1e-9)), 1)[0] for i in range(3)]
            print(f"    alpha {alpha:4.2f}   slopes  i=1 {sl[0]:.3f}   i=2 {sl[1]:.3f}   i=3 {sl[2]:.3f}")
