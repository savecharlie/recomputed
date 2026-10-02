"""Does the counterexample need three layers at FINITE size?

The paper's reason for three is asymptotic: a cycle A->B->C->A needs three
states, so a two-layer system can carry no inter-layer driving at all, and then
the Q spectrum is 0 and -(1-eps)/(1-p) = -1 + eps - O(eps^2), which sits ABOVE
the k=1 mode at -1 + O(eps^2).  lambda_2 of R is then a REAL eigenvalue,
lambda_I = 0, and the conjecture's precondition is never met.

At the 72-state instance eps = 0.305, nowhere near asymptotic, so the ordering
could in principle come out differently.  Brute force over the two-layer box.

(The first version of this file diagonalised the full 2L x 2L matrix and was
still running after six minutes.  The two-layer model is shift invariant too,
so it block-diagonalises into L 2x2 blocks and the whole sweep takes seconds.
Same lesson as every other instrument here: use the structure.)
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ.setdefault(v, "2")
import sys, numpy as np
sys.path.insert(0, ".")

def spectrum(L, eps, p, K, wAB):
    th = 2 * np.pi / L
    rAp = eps / (2 * (1 - np.cos(th))) + 1 / np.sin(th)
    rAm = eps / (2 * (1 - np.cos(th))) - 1 / np.sin(th)
    if rAm <= 0: return None
    rB = K / (2 * (1 - np.cos(th)))
    wBA = wAB * p / (1 - p)                 # so that p^ss = (p, 1-p) per site
    Q = np.array([[-wAB, wBA], [wAB, -wBA]], dtype=complex)
    phi = 2 * np.pi * np.arange(L) / L
    blocks = np.broadcast_to(Q, (L, 2, 2)).copy()
    for lay, (rp, rm) in enumerate([(rAp, rAm), (rB, rB)]):
        blocks[:, lay, lay] += (-(rp + rm) * (1 - np.cos(phi))
                                - 1j * (rp - rm) * np.sin(phi))
    w = np.linalg.eigvals(blocks).ravel()
    sigma = L * (p / L) * (rAp - rAm) * np.log(rAp / rAm)   # only the driven ring
    return w, sigma

rng = np.random.default_rng(77)
tried = cohr = viol = 0
best_lI = 0.0
examples = []
for _ in range(120000):
    L = int(np.exp(rng.uniform(np.log(8), np.log(800))))
    lo = 2 * np.tan(np.pi / L)
    if lo >= 1: continue
    eps = lo * np.exp(rng.uniform(0, 3.0))
    if eps >= 1: continue
    p = np.exp(rng.uniform(np.log(1e-6), np.log(0.49)))
    K = np.exp(rng.uniform(np.log(1.05), np.log(800)))
    wAB = np.exp(rng.uniform(-7, 5))
    out = spectrum(L, eps, p, K, wAB)
    if out is None: continue
    w, sigma = out
    tried += 1
    w = w[np.argsort(-w.real)]
    l2 = w[1]; lR, lI = -l2.real, abs(l2.imag)
    best_lI = max(best_lI, lI)
    if lI < lR or lI < 1e-9: continue
    cohr += 1
    if sigma < lI ** 2 / lR:
        viol += 1
        if len(examples) < 5:
            examples.append((L, eps, p, K, wAB, sigma, lR, lI))

print(__doc__)
print(f"  two-layer parameter points inside the model : {tried}")
print(f"  with lambda_I >= lambda_R (precondition met) : {cohr}")
print(f"  violating OSB                               : {viol}")
print(f"  largest |Im lambda_2| seen anywhere          : {best_lI:.3e}")
print()
if cohr == 0:
    print("  => the second eigenvalue of EVERY two-layer system in this sweep is")
    print("     real.  The oscillating k=1 mode is always buried under the real")
    print("     inter-layer relaxation mode, exactly as the paper argues -- and")
    print("     the argument is not merely asymptotic: it holds at eps of order 1")
    print("     too.  Three layers are needed at finite size as well, which is a")
    print("     stronger statement than the paper makes.")
else:
    print("  => the three-layer requirement is NOT absolute at finite size:")
    for e in examples:
        print(f"     L={e[0]} eps={e[1]:.4f} p={e[2]:.3e} K={e[3]:.3f} w={e[4]:.3e}"
              f" sigma={e[5]:.4f} lR={e[6]:.4f} lI={e[7]:.4f}")
