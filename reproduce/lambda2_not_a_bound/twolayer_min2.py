"""Smallest two-layer violation, as a function of how far INSIDE the
conjecture's domain it is required to sit.

The first run returned L = 22 (44 states) with lambda_I = lambda_R to five
digits -- i.e. the coherent number exactly at 1/2pi, the exact boundary of the
regime where the conjecture even applies.  A result sitting on the boundary of
its own domain is not wrong, but it is fragile and must be reported as such.
So: require lambda_I >= (1+delta) lambda_R and watch the minimum move.
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ.setdefault(v, "2")
import numpy as np
from scipy.optimize import minimize

def ev(L, eps, p, K, wAB, delta):
    if not (0 < eps < 1 and 0 < p < 1 and K > 0 and wAB > 0): return np.inf, None
    if eps <= 2 * np.tan(np.pi / L): return np.inf, None
    th = 2 * np.pi / L
    rAp = eps / (2*(1-np.cos(th))) + 1/np.sin(th)
    rAm = eps / (2*(1-np.cos(th))) - 1/np.sin(th)
    rB = K / (2*(1-np.cos(th)))
    wBA = wAB * p / (1 - p)
    Q = np.array([[-wAB, wBA], [wAB, -wBA]], dtype=complex)
    phi = 2*np.pi*np.arange(L)/L
    bl = np.broadcast_to(Q, (L,2,2)).copy()
    for lay, (rp, rm) in enumerate([(rAp, rAm), (rB, rB)]):
        bl[:, lay, lay] += -(rp+rm)*(1-np.cos(phi)) - 1j*(rp-rm)*np.sin(phi)
    w = np.linalg.eigvals(bl).ravel(); w = w[np.argsort(-w.real)]
    l2 = w[1]; lR, lI = -l2.real, abs(l2.imag)
    if lR <= 0 or lI < (1+delta)*lR or lI < 1e-9: return np.inf, None
    sig = p * (rAp - rAm) * np.log(rAp/rAm)
    if sig <= 0: return np.inf, None
    return -(lI**2/lR)/sig, (eps, p, K, wAB, sig, lR, lI)

print(__doc__)
print(f"{'delta':>7} {'min L':>7} {'states':>7} {'ratio':>9} {'lI/lR':>8}  parameters")
rng = np.random.default_rng(5150)
for delta in (0.0, 0.01, 0.05, 0.2, 0.5, 1.0):
    hit = None
    for L in list(range(8, 60)) + list(range(60, 500, 3)):
        if 2*np.tan(np.pi/L) >= 1: continue
        best = (np.inf, None)
        for _ in range(45):
            lo = 2*np.tan(np.pi/L)
            x0 = [lo*np.exp(rng.uniform(0,3)),
                  np.exp(rng.uniform(np.log(1e-5), np.log(0.49))),
                  np.exp(rng.uniform(np.log(1.05), np.log(800))),
                  np.exp(rng.uniform(-7,5))]
            r = minimize(lambda x: ev(L, *x, delta)[0], x0, method="Nelder-Mead",
                         options=dict(maxiter=900, fatol=1e-12))
            v, pk = ev(L, *r.x, delta)
            if v < best[0]: best = (v, pk)
        if best[1] is not None and -best[0] > 1:
            hit = (L, -best[0], best[1]); break
    if hit:
        L, rt, (eps, p, K, w, sig, lR, lI) = hit
        print(f"{delta:7.2f} {L:7d} {2*L:7d} {rt:9.5f} {lI/lR:8.4f}  "
              f"eps={eps:.5f} p={p:.3e} K={K:.3f} w={w:.3e}")
    else:
        print(f"{delta:7.2f}       -       -         -        -  none below L=500")
