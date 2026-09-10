"""The estimator returns 0.945 on a process whose answer is exact.  Look at WHERE."""
import numpy as np, math
tau, dtr, sigma = 20.0, 0.1, 0.015811393924397392
NT, NENS = 12000, 512
rng = np.random.default_rng(11)
a = math.exp(-dtr/tau); s = sigma*math.sqrt(1-a*a)
R = np.empty((NT, NENS)); R[0] = sigma*rng.standard_normal(NENS)
for i in range(1, NT):
    R[i] = a*R[i-1] + s*rng.standard_normal(NENS)
gd = R - R.mean()
gvar = float((gd*gd).mean())
print("var/th          %.4f" % (gvar/sigma**2))
lags = [0,1,10,50,100,200,400,800,1200,1600]
print(" lag t    est/var     exp(-t/tau)   ratio")
for L in lags:
    c = gvar if L==0 else float((gd[:-L]*gd[L:]).mean())
    ex = math.exp(-L*dtr/tau)
    print("%6.1f   %9.5f   %9.5f   %7.4f" % (L*dtr, c/gvar, ex, (c/gvar)/ex))
# full integral, exact truth = sigma^2 * tau
Lc = 1600
cov = np.empty(Lc+1); cov[0]=gvar
for L in range(1,Lc+1): cov[L]=float((gd[:-L]*gd[L:]).mean())
zc = np.nonzero(cov<=0)[0]; cut = int(zc[0]) if len(zc) else len(cov)
print("\nfirst zero crossing at t = %.1f (%.2f tau)" % (cut*dtr, cut*dtr/tau))
for name, c in (("to first zero", cut), ("to 8 tau", Lc+1)):
    D = float(np.trapz(cov[:c], dx=dtr))
    print("  %-14s D/th = %.4f" % (name, D/(sigma**2*tau)))
D_analytic_disc = dtr*(gvar*(1+a)/(1-a) - gvar)/1.0
print("  exact discrete sum  D/th = %.4f  (dtr*sum_k a^k form)"
      % ((dtr*gvar*(0.5 + a/(1-a)))/(sigma**2*tau)))
