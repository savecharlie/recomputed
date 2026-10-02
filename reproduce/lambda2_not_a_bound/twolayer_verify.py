"""Witness for a TWO-layer violation.  Dense, independent of the block route,
60-digit arithmetic, every structural claim checked by hand.

If this holds, two layers suffice at finite size -- which the paper's
"this is why our counterexample consists of three layers, not two" does not
claim but is easy to read as claiming.  The paper's argument is a Gershgorin
bound and Gershgorin is loose; at eps ~ 0.07 the real ordering can differ.
"""
import sys, numpy as np, mpmath as mp
mp.mp.dps = 60

L, eps, p, K, wAB = 216, 0.0747, 6.851e-02, 4.521, 9.119e-01

E, P_, KK, W = (mp.mpf(repr(x)) for x in (eps, p, K, wAB))
th = 2 * mp.pi / L
rAp = E / (2 * (1 - mp.cos(th))) + 1 / mp.sin(th)
rAm = E / (2 * (1 - mp.cos(th))) - 1 / mp.sin(th)
rB = KK / (2 * (1 - mp.cos(th)))
wBA = W * P_ / (1 - P_)
print(__doc__)
print(f"L={L}  ->  {2*L} states")
print(f"  r+_A = {mp.nstr(rAp,12)}   r-_A = {mp.nstr(rAm,12)}   r_B = {mp.nstr(rB,12)}")
print(f"  A->B = {mp.nstr(W,12)}     B->A = {mp.nstr(wBA,12)}")
assert min(rAp, rAm, rB, W, wBA) > 0
print("  all rates strictly positive.")

n = 2 * L
R = [[mp.mpf(0)] * n for _ in range(n)]
for lay, (rp, rm) in enumerate([(rAp, rAm), (rB, rB)]):
    for i in range(L):
        j = 2 * i + lay
        R[2 * ((i + 1) % L) + lay][j] += rp
        R[2 * ((i - 1) % L) + lay][j] += rm
        R[j][j] -= rp + rm
for i in range(L):
    R[2*i+1][2*i+0] += W;    R[2*i+0][2*i+0] -= W
    R[2*i+0][2*i+1] += wBA;  R[2*i+1][2*i+1] -= wBA
print(f"  max |column sum| = "
      f"{mp.nstr(max(abs(sum(R[i][j] for i in range(n))) for j in range(n)), 6)}")

pss = [(P_ / L) if j % 2 == 0 else ((1 - P_) / L) for j in range(n)]
print(f"  sum p^ss       = {mp.nstr(sum(pss), 20)}")
print(f"  max |(R p)_i|  = "
      f"{mp.nstr(max(abs(sum(R[i][j]*pss[j] for j in range(n))) for i in range(n)), 6)}")

sig = mp.mpf(0)
for i in range(n):
    for j in range(n):
        if i == j: continue
        f = R[i][j] * pss[j]
        if f > 0:
            g = R[j][i] * pss[i]; assert g > 0
            sig += f * mp.log(f / g)
print(f"\n  sigma_ss (dense double sum, 60 digits) = {mp.nstr(sig, 20)}")
print(f"  sigma from the ring formula            = "
      f"{mp.nstr(L * (P_/L) * (rAp - rAm) * mp.log(rAp/rAm), 20)}")

Rf = np.array([[float(x) for x in row] for row in R])
w = np.linalg.eigvals(Rf); w = w[np.argsort(-w.real)]
print(f"\n  dense {n}x{n} float spectrum, top 5:")
for z in w[:5]: print(f"     {z.real:+.14f} {z.imag:+.14f}j")

phi1 = th
B = [[-W, wBA], [W, -wBA]]
for lay, (rp, rm) in enumerate([(rAp, rAm), (rB, rB)]):
    B[lay][lay] += -(rp+rm)*(1-mp.cos(phi1)) - mp.mpc(0,1)*(rp-rm)*mp.sin(phi1)
tr = B[0][0] + B[1][1]
dt = B[0][0]*B[1][1] - B[0][1]*B[1][0]
roots = mp.polyroots([mp.mpc(1), -tr, dt], maxsteps=200, extraprec=200)
z = max(roots, key=lambda r: r.real)
lR, lI = -z.real, abs(z.imag)
print(f"\n  lambda_2 from the k=1 block quadratic (60 digits):")
print(f"     {mp.nstr(z, 22)}")
print(f"     agrees with the dense float eig to "
      f"{mp.nstr(min(abs(z-mp.mpc(w[1].real,w[1].imag)), abs(mp.conj(z)-mp.mpc(w[1].real,w[1].imag))),4)}")
rhs = lI**2/lR
print(f"\n  lambda_R = {mp.nstr(lR,18)}")
print(f"  lambda_I = {mp.nstr(lI,18)}    precondition lambda_I>=lambda_R: {lI>=lR}")
print(f"  lambda_I^2/lambda_R = {mp.nstr(rhs,18)}")
print(f"  sigma_ss            = {mp.nstr(sig,18)}")
print(f"\n  DEFICIT = {mp.nstr(rhs-sig,10)}   RATIO = {mp.nstr(rhs/sig,10)}")
N = lI/(2*mp.pi*lR); DS = 2*mp.pi*sig/lI
print(f"  Delta S = {mp.nstr(DS,10)}  <  4 pi^2 N = {mp.nstr(4*mp.pi**2*N,10)} ?"
      f"  {DS < 4*mp.pi**2*N}")
print(f"\n  {'*** OSB VIOLATED with TWO layers, '+str(n)+' states ***' if rhs > sig else 'holds'}")
