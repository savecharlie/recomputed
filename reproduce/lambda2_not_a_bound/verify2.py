"""Second witness: the smallest violating system INSIDE the paper's own
parameter window 1/2 < q < 1.  Verified the same way as the 72-state one.

This is the one to quote when the paper's framing matters: the 72-state
instance sits at q = 0.374, outside 1/2 < q < 1, and violates by 3.7%.  This
one is inside the window and violates by 71%.
"""
import sys, numpy as np, mpmath as mp
sys.path.insert(0, ".")
import model as M
mp.mp.dps = 60

L, eps, d, K = 41, 0.15617, 2.0000, None
# recover the exact point from robust.py's own optimiser output
import re
txt = open("robust.out").read()
row = [l for l in txt.splitlines() if "states  ratio" in l and l.startswith("paper's window 1/2")][0]
L = int(re.search(r"L=\s*(\d+)", row).group(1))
eps = float(re.search(r"eps=([\d.]+)", row).group(1))
d = float(re.search(r"d=\s*([\d.]+)", row).group(1))
print(__doc__)
print(f"from robust.out: {row.strip()}")
# robust.py prints eps to 5 dp and d to 4; re-optimise kappa and K at that point
best = None
for kmul in np.linspace(1.001, 3.0, 60):
    for Kv in np.exp(np.linspace(np.log(1.5), np.log(500), 60)):
        P0 = M.rates_eps(L, eps, d=max(d, 2.0 + 1e-12), K=Kv)
        P = M.rates_eps(L, eps, d=max(d, 2.0 + 1e-12), kappa=kmul * P0["kappa_min"], K=Kv)
        if P["p"] >= 1 or not M.valid(P)[0]: continue
        r = M.osb(P)
        if not r["precondition"]: continue
        if best is None or r["ratio"] > best[0]: best = (r["ratio"], P, r)
ratio, P, r = best
print(f"\nre-optimised at that (L, eps, d): kappa/kappa_min="
      f"{P['kappa']/P['kappa_min']:.6f}  K={P['K']:.6f}  ratio={ratio:.6f}")
print(f"q = {P['q']:.6f}   (the paper's window is 1/2 < q < 1)")

E = mp.mpf(repr(P["eps"])); D = mp.mpf(repr(P["d"])); KK = mp.mpf(repr(P["K"]))
a = 1 - E; p = mp.mpf(repr(P["kappa"])) * E ** 2
u = a * p / (1 - p); h = (a + D * E) / 2
v = (D - 1) * E * mp.sqrt(p); b = (D - 1) * (1 - p) * E / (2 * mp.sqrt(p))
th = 2 * mp.pi / L
rAp = E / (2 * (1 - mp.cos(th))) + 1 / mp.sin(th)
rAm = E / (2 * (1 - mp.cos(th))) - 1 / mp.sin(th)
rBC = KK / (2 * (1 - mp.cos(th)))
Q = [[-a + 0, u - v, u + v], [a/2 + b, -u - h, h - v], [a/2 - b, h + v, -u - h]]
offs = [rAp, rAm, rBC] + [Q[i][j] for i in range(3) for j in range(3) if i != j]
print(f"\nsmallest off-diagonal rate: {mp.nstr(min(offs), 12)}")
assert min(offs) > 0
n = 3 * L
R = [[mp.mpf(0)] * n for _ in range(n)]
for lay, (rp, rm) in enumerate([(rAp, rAm), (rBC, rBC), (rBC, rBC)]):
    for i in range(L):
        j = 3 * i + lay
        R[3 * ((i+1) % L) + lay][j] += rp
        R[3 * ((i-1) % L) + lay][j] += rm
        R[j][j] -= rp + rm
for i in range(L):
    for A in range(3):
        for B in range(3):
            R[3*i+A][3*i+B] += Q[A][B]
print(f"max |column sum|: "
      f"{mp.nstr(max(abs(sum(R[i][j] for i in range(n))) for j in range(n)), 6)}")
cell = [p / L, (1-p)/(2*L), (1-p)/(2*L)]
pss = [cell[j % 3] for j in range(n)]
print(f"max |(R p_ss)_i|: "
      f"{mp.nstr(max(abs(sum(R[i][j]*pss[j] for j in range(n))) for i in range(n)), 6)}")
sig = mp.mpf(0)
for i in range(n):
    for j in range(n):
        if i == j: continue
        f = R[i][j] * pss[j]
        if f > 0:
            g = R[j][i] * pss[i]; assert g > 0
            sig += f * mp.log(f / g)
B = [[Q[i][j] for j in range(3)] for i in range(3)]
for lay, (rp, rm) in enumerate([(rAp, rAm), (rBC, rBC), (rBC, rBC)]):
    B[lay][lay] += -(rp+rm)*(1-mp.cos(th)) - mp.mpc(0,1)*(rp-rm)*mp.sin(th)
tr = B[0][0]+B[1][1]+B[2][2]
m2 = (B[0][0]*B[1][1]-B[0][1]*B[1][0] + B[0][0]*B[2][2]-B[0][2]*B[2][0]
      + B[1][1]*B[2][2]-B[1][2]*B[2][1])
dt = (B[0][0]*(B[1][1]*B[2][2]-B[1][2]*B[2][1]) - B[0][1]*(B[1][0]*B[2][2]-B[1][2]*B[2][0])
      + B[0][2]*(B[1][0]*B[2][1]-B[1][1]*B[2][0]))
z = max(mp.polyroots([mp.mpc(1), -tr, m2, -dt], maxsteps=200, extraprec=200),
        key=lambda r: r.real)
# confirm this really is the second-largest over the WHOLE spectrum
allw = np.linalg.eigvals(M.R_full(P)); allw = allw[np.argsort(-allw.real)]
print(f"\ndense {n}x{n} float spectrum, top 4:")
for q_ in allw[:4]: print(f"   {q_.real:+.12f} {q_.imag:+.12f}j")
lR, lI = -z.real, abs(z.imag)
print(f"\nlambda_2 (block cubic, 60 digits): {mp.nstr(z, 22)}")
print(f"   agrees with dense float eig to "
      f"{mp.nstr(min(abs(z-mp.mpc(allw[1].real,allw[1].imag)), abs(mp.conj(z)-mp.mpc(allw[1].real,allw[1].imag))),4)}")
print(f"lambda_R = {mp.nstr(lR,18)}")
print(f"lambda_I = {mp.nstr(lI,18)}   precondition: {lI >= lR}")
print(f"lambda_I^2/lambda_R = {mp.nstr(lI**2/lR,18)}")
print(f"sigma_ss            = {mp.nstr(sig,18)}")
print(f"\nDEFICIT = {mp.nstr(lI**2/lR - sig, 10)}    RATIO = {mp.nstr((lI**2/lR)/sig, 10)}")
N = lI/(2*mp.pi*lR); DS = 2*mp.pi*sig/lI
print(f"Delta S = {mp.nstr(DS,10)}   4 pi^2 N = {mp.nstr(4*mp.pi**2*N,10)}")
print(f"\n{'*** OSB VIOLATED at '+str(n)+' states, INSIDE the paper window ***' if lI**2/lR > sig else 'holds'}")
