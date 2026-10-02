"""How much of the 72-state result depends on corners of the parameter space?

The minimiser went to d = 2 + 1e-15 (the paper requires d > 2) and to
q = -ln(eps)/ln(L) = 0.374, which is OUTSIDE the paper's window 1/2 < q < 1.
Both are legitimate as Markov chains -- at a fixed L nothing in the
construction refers to q, and d > 2 is satisfied -- but a result that lives
only on a boundary should be reported as living on a boundary.  So: the same
minimisation under progressively stricter constraints.
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "2")
import sys, numpy as np
sys.path.insert(0, ".")
import model as M
from scipy.optimize import minimize

RNG = np.random.default_rng(20261002)

def feasible(L, eps, d, kmul, K, dmin, qlo, qhi):
    if not (0.0 < eps < 1.0) or d < dmin or kmul < 1.0 or K <= 0.0:
        return None
    if eps <= 2.0 * np.tan(np.pi / L):
        return None
    q = -np.log(eps) / np.log(L)
    if not (qlo <= q <= qhi):
        return None
    P = M.rates_eps(L, eps, d=d, kappa=None, K=K)
    P = M.rates_eps(L, eps, d=d, kappa=kmul * P["kappa_min"], K=K)
    if P["p"] >= 1.0 or not M.valid(P)[0]:
        return None
    return P

def best(L, dmin=2.0, qlo=-np.inf, qhi=np.inf, dfix=None, n=50):
    def neg(x):
        d = dfix if dfix is not None else x[1]
        P = feasible(L, x[0], d, x[2], x[3], dmin, qlo, qhi)
        if P is None: return np.inf
        r = M.osb(P)
        if not r["precondition"] or not np.isfinite(r["ratio"]): return np.inf
        return -r["ratio"]
    lo = 2.0 * np.tan(np.pi / L)
    if lo >= 1.0: return None
    bst = (np.inf, None)
    for _ in range(n):
        x0 = [lo + (1 - lo) * RNG.random() ** 2,
              max(dmin, 2.0) + 8 * RNG.random(), 1 + 2 * RNG.random(),
              10 ** (3 * RNG.random())]
        r = minimize(neg, x0, method="Nelder-Mead",
                     options=dict(maxiter=800, xatol=1e-11, fatol=1e-13))
        if r.fun < bst[0]: bst = (r.fun, r.x)
    if bst[1] is None or not np.isfinite(bst[0]): return None
    x = bst[1]
    d = dfix if dfix is not None else x[1]
    P = feasible(L, x[0], d, x[2], x[3], dmin, qlo, qhi)
    return (-bst[0], P) if P is not None else None

def smallest(label, **kw):
    for L in range(7, 400):
        b = best(L, **kw)
        if b and b[0] > 1.0:
            ratio, P = b
            r = M.osb(P)
            print(f"{label:<44} L={L:4d}  {3*L:5d} states  ratio={ratio:7.4f}"
                  f"  eps={P['eps']:.5f}  d={P['d']:7.4f}  q={P['q']:.4f}"
                  f"  sigma={r['sigma']:.5f}  rhs={r['rhs']:.5f}")
            return L, P
    print(f"{label:<44} none below L=400")
    return None

print(__doc__)
print("Smallest violating system under each constraint set:\n")
smallest("free (d>2, eps free)",                  dmin=2.0)
smallest("d >= 2.1",                              dmin=2.1)
smallest("d >= 2.5",                              dmin=2.5)
smallest("d = 3 exactly (the paper's choice)",    dfix=3.0, dmin=2.0)
smallest("paper's window 1/2 < q < 1",            dmin=2.0, qlo=0.5, qhi=1.0)
smallest("paper's window AND d = 3",              dfix=3.0, dmin=2.0, qlo=0.5, qhi=1.0)
print()
print("Reading: the conjecture is false on objects small enough to have been")
print("searched exhaustively.  The very smallest instances sit against d -> 2")
print("and outside the paper's q window; both qualifications are real and both")
print("leave the headline intact, because every row above is a genuine Markov")
print("generator with strictly positive rates.")
