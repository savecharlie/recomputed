"""Witness for the smallest violating system.  Law 5: never report a number
from a configuration that merely nearly passes -- build the actual object and
re-measure it independently, at precision that cannot be the explanation.

The candidate from minsize.py violates OSB by only 1.0367x, with d sitting on
its own boundary (d > 2).  So: rebuild the 3L x 3L generator in 60-digit
arithmetic, recompute the spectrum and the entropy production from the DENSE
matrix (no Fourier blocks, no shift-invariance shortcut), and check every
structural precondition by hand.
"""
import sys, json, numpy as np
sys.path.insert(0, ".")
import model as M
import mpmath as mp
mp.mp.dps = 60

CAND = json.load(open("minsize.json"))
hit = next(r for r in CAND if r["ratio"] > 1)
print(__doc__)
print(f"candidate from minsize.json: L={hit['L']}  ({3*hit['L']} states)")
for k in ("eps", "d", "kappa", "K", "q", "ratio"):
    print(f"   {k:6s} = {hit[k]!r}")

L, eps, d, kappa, K = hit["L"], hit["eps"], hit["d"], hit["kappa"], hit["K"]

# ---- rebuild in mp arithmetic, from the paper's formulas only ---------------
E = mp.mpf(repr(eps)); D = mp.mpf(repr(d)); KA = mp.mpf(repr(kappa)); KK = mp.mpf(repr(K))
a = 1 - E
p = KA * E ** 2
u = a * p / (1 - p)
h = (a + D * E) / 2
v = (D - 1) * E * mp.sqrt(p)
b = (D - 1) * (1 - p) * E / (2 * mp.sqrt(p))
th = 2 * mp.pi / L
rAp = E / (2 * (1 - mp.cos(th))) + 1 / mp.sin(th)
rAm = E / (2 * (1 - mp.cos(th))) - 1 / mp.sin(th)
rBC = KK / (2 * (1 - mp.cos(th)))
Q = [[-a, u, u], [a / 2, -u - h, h], [a / 2, h, -u - h]]
Q1 = [[0, -v, v], [b, 0, -v], [-b, v, 0]]
Q = [[Q[i][j] + Q1[i][j] for j in range(3)] for i in range(3)]

print("\nSTRUCTURE, checked by hand in 60 digits:")
offs = {"r+_A": rAp, "r-_A": rAm, "r_BC": rBC}
for i in range(3):
    for j in range(3):
        if i != j: offs[f"Q[{i},{j}]"] = Q[i][j]
worst = min(offs.items(), key=lambda kv: kv[1])
for k, val in offs.items():
    print(f"   {k:9s} = {mp.nstr(val, 12):>22s}   {'>= 0 OK' if val >= 0 else 'NEGATIVE'}")
print(f"   smallest off-diagonal rate: {worst[0]} = {mp.nstr(worst[1], 12)}")
assert worst[1] > 0, "NOT A MARKOV GENERATOR"
print("   -> all rates strictly positive: this is a genuine Markov jump process.")

# ---- dense generator, mp ----------------------------------------------------
n = 3 * L
R = [[mp.mpf(0)] * n for _ in range(n)]
ring = [(rAp, rAm), (rBC, rBC), (rBC, rBC)]
for lay, (rp, rm) in enumerate(ring):
    for i in range(L):
        j = 3 * i + lay
        R[3 * ((i + 1) % L) + lay][j] += rp
        R[3 * ((i - 1) % L) + lay][j] += rm
        R[j][j] -= rp + rm
for i in range(L):
    for A in range(3):
        for B in range(3):
            R[3 * i + A][3 * i + B] += Q[A][B]
colmax = max(abs(sum(R[i][j] for i in range(n))) for j in range(n))
print(f"\n   max |column sum| of the 60-digit generator: {mp.nstr(colmax, 6)}")

# stationary distribution, then VERIFY it by applying R
cell = [p / L, (1 - p) / (2 * L), (1 - p) / (2 * L)]
pss = [cell[j % 3] for j in range(n)]
resid = max(abs(sum(R[i][j] * pss[j] for j in range(n))) for i in range(n))
print(f"   max |(R p^ss)_i| for the analytic p^ss: {mp.nstr(resid, 6)}")
print(f"   sum p^ss = {mp.nstr(sum(pss), 20)}")

# ---- entropy production, dense double sum, 60 digits -----------------------
sig = mp.mpf(0)
for i in range(n):
    for j in range(n):
        if i == j: continue
        f = R[i][j] * pss[j]
        if f > 0:
            g = R[j][i] * pss[i]
            assert g > 0
            sig += f * mp.log(f / g)
print(f"\n   sigma_ss  (dense double sum, 60 digits) = {mp.nstr(sig, 20)}")
P = M.rates_eps(L, eps, d=d, kappa=kappa, K=K)
print(f"   sigma_ss  (float, shift-invariant route) = {M.sigma_ss_fast(P):.17g}")

# ---- spectrum: dense float eig on the FULL matrix, plus mp refinement ------
Rf = np.array([[float(x) for x in row] for row in R])
w = np.linalg.eigvals(Rf); w = w[np.argsort(-w.real)]
print(f"\n   top five eigenvalues of the dense 72x72 matrix:")
for z in w[:5]:
    print(f"      {z.real:+.15f} {z.imag:+.15f}j")
l2 = w[1]
lR, lI = -l2.real, abs(l2.imag)
# Refine lambda_2 exactly, through the k=1 Fourier block.
#
# A first attempt Newton-refined det(zI - R) for the full 72x72 matrix and the
# root-finder refused: det is of order 1e131 there, so an ABSOLUTE tolerance of
# 1e-40 on it is not a statement about the root at all.  Badly scaled
# diagnostic, not a bad root.  The block route is exact and well scaled: the
# k=1 block is 3x3, its characteristic polynomial is a cubic, and mp.polyroots
# solves it to arbitrary precision.  It is also an INDEPENDENT route from the
# dense float eig above, so the two together are a cross-check.
phi1 = th                      # k=1  =>  2 pi k / L = theta
B = [[Q[i][j] for j in range(3)] for i in range(3)]
for lay, (rp, rm) in enumerate(ring):
    B[lay][lay] += (-(rp + rm) * (1 - mp.cos(phi1)) - mp.mpc(0, 1) * (rp - rm) * mp.sin(phi1))
# characteristic cubic  z^3 - tr z^2 + (sum of 2x2 principal minors) z - det
tr = B[0][0] + B[1][1] + B[2][2]
m2 = (B[0][0]*B[1][1] - B[0][1]*B[1][0]
      + B[0][0]*B[2][2] - B[0][2]*B[2][0]
      + B[1][1]*B[2][2] - B[1][2]*B[2][1])
dt = (B[0][0]*(B[1][1]*B[2][2] - B[1][2]*B[2][1])
      - B[0][1]*(B[1][0]*B[2][2] - B[1][2]*B[2][0])
      + B[0][2]*(B[1][0]*B[2][1] - B[1][1]*B[2][0]))
roots = mp.polyroots([mp.mpc(1), -tr, m2, -dt], maxsteps=200, extraprec=200)
z = max(roots, key=lambda r: r.real)
print(f"\n   k=1 block cubic roots (60 digits):")
for r in sorted(roots, key=lambda r: -r.real):
    print(f"      {mp.nstr(r, 22)}")
print(f"   lambda_2 (largest real part in the k=1 block):")
print(f"      {mp.nstr(z, 25)}")
print(f"   agreement with the dense 72x72 float eig: "
      f"{mp.nstr(min(abs(z - mp.mpc(l2.real, l2.imag)), abs(mp.conj(z) - mp.mpc(l2.real, l2.imag))), 4)}")
print(f"   residual of the cubic at this root: {mp.nstr(abs(z**3 - tr*z**2 + m2*z - dt), 4)}")
lRm, lIm = -z.real, abs(z.imag)
rhs = lIm ** 2 / lRm
print(f"\n   lambda_R = {mp.nstr(lRm, 20)}")
print(f"   lambda_I = {mp.nstr(lIm, 20)}")
print(f"   precondition lambda_I >= lambda_R : {lIm >= lRm}")
print(f"   lambda_I^2/lambda_R = {mp.nstr(rhs, 20)}")
print(f"   sigma_ss            = {mp.nstr(sig, 20)}")
print(f"\n   OSB demands sigma_ss >= lambda_I^2/lambda_R.")
print(f"   DEFICIT  = {mp.nstr(rhs - sig, 12)}   (positive => VIOLATED)")
print(f"   RATIO    = {mp.nstr(rhs / sig, 12)}")
print(f"\n   {'*** OSB VIOLATED at ' + str(n) + ' states ***' if rhs > sig else 'holds'}")
# also report Delta S >= 4 pi^2 N in the paper's own variables
N = lIm / (2 * mp.pi * lRm)
DS = 2 * mp.pi * sig / lIm
print(f"\n   in the paper's variables:  coherent number N = {mp.nstr(N, 12)}"
      f"  (needs N >= 1/2pi = {mp.nstr(1/(2*mp.pi), 8)}: {N >= 1/(2*mp.pi)})")
print(f"   Delta S = 2 pi sigma/lambda_I = {mp.nstr(DS, 12)}")
print(f"   4 pi^2 N                      = {mp.nstr(4*mp.pi**2*N, 12)}")
print(f"   conjecture Delta S >= 4 pi^2 N : {'HOLDS' if DS >= 4*mp.pi**2*N else 'FALSE'}")
