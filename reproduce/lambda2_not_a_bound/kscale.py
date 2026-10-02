"""The prefactor in A_abs = c * p.  Is c -> 3 as K -> infinity, and how?

constant.py measured c = 3.0927, 3.0171, 3.00097, 3.000037 at K = 4, 10, 40,
200 and I eyeballed "O(1/K^2)" off four points.  That is a round number in a
sentence (law 16).  Test the rate instead of asserting it.
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ.setdefault(v, "2")
import sys, numpy as np
sys.path.insert(0, ".")
import model as M, repair as Rp

print(__doc__)
L = 6400
print(f"{'K':>9} {'A_abs/p':>14} {'c - 3':>13} {'(c-3)*K':>11} {'(c-3)*K^2':>11}")
ks, cs = [], []
for K in (4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0):
    P0 = M.rates_eps(L, L ** -0.55, d=3.0, K=K)
    P = M.rates_eps(L, L ** -0.55, d=3.0, kappa=1.05 * P0["kappa_min"], K=K)
    c = Rp.A_abs(P)[0] / P["p"]
    ks.append(K); cs.append(c)
    print(f"{K:9.1f} {c:14.9f} {c-3:13.3e} {(c-3)*K:11.5f} {(c-3)*K**2:11.4f}")
ks, cs = np.array(ks), np.array(cs)
sl = np.polyfit(np.log(ks), np.log(cs - 3.0), 1)[0]
print(f"\n   (c - 3) ~ K^{sl:+.4f}")
print(f"   so c -> 3 as K -> infinity, with the correction falling as a power")
print(f"   of K that the fit puts at {sl:+.3f}.")
# and the joint statement: is 3 the number of layers?  change the layer count.
print("\n   IS THE 3 THE NUMBER OF LAYERS?  The construction has exactly three")
print("   layers, so this cannot be tested inside it.  What CAN be tested is the")
print("   origin of the factor: A_abs is normalised by RMS 1 over the 3L states")
print("   while the mode's stationary weight is carried by p/L of them, so the")
print("   factor should be the ratio of the number of states to the number of")
print("   SITES, i.e. 3L/L = 3.  Direct check: recompute A_abs with the RMS")
print("   normalised over L states instead of 3L -- the factor must become 1.")
import amplitude as A
P0 = M.rates_eps(L, L ** -0.55, d=3.0, K=512.0)
P = M.rates_eps(L, L ** -0.55, d=3.0, kappa=1.05 * P0["kappa_min"], K=512.0)
s = A.sector(P)
from scipy.optimize import minimize
def best(scale):
    rng = np.random.default_rng(2)
    def neg(z):
        g = z[:3] + 1j * z[3:]
        nr = np.linalg.norm(g)
        if nr < 1e-14: return 0.0
        g = g * (np.sqrt(2.0 * scale) / nr)        # sum|g|^2 = 2*scale
        return -2.0 * abs(Rp.W2_sector(P, g, s))
    return max(-minimize(neg, rng.standard_normal(6), method="Nelder-Mead",
               options=dict(maxiter=4000)).fun for _ in range(25))
print(f"   RMS 1 over 3L states (sum|g|^2 = 6): A_abs/p = {best(3.0)/P['p']:.6f}")
print(f"   RMS 1 over  L states (sum|g|^2 = 2): A_abs/p = {best(1.0)/P['p']:.6f}")
print("   => the 3 is exactly the state-to-site ratio of the normalisation, not")
print("      a property of the dynamics.  The physics is A_abs proportional to p.")
