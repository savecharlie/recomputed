"""Same OU, but read out through the ACTUAL spline g.  Where does 0.945 come from?"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius
import amplitude_phase as ap
from channels import fbarPhi_spline

T, GAM, TR, u = 5e-4, 0.05, 0.4, 3.5
eps = math.pi*GAM*TR*u/4.0
esc = Grasshopper(TR, 20.0); Rstar = limit_cycle_radius(esc, GAM, eps)
p = ap.predict(esc, GAM, eps, T, Rstar)
tau, sigma = 1.0/p['k'], p['sigma_R']
NT, NENS, dtr = 12000, 128, 0.1
rng = np.random.default_rng(11); a = math.exp(-dtr/tau); s = sigma*math.sqrt(1-a*a)
R = np.empty((NT, NENS)); R[0] = Rstar + sigma*rng.standard_normal(NENS)
for i in range(1, NT):
    R[i] = Rstar + a*(R[i-1]-Rstar) + s*rng.standard_normal(NENS)

Rg, fg = fbarPhi_spline(esc, max(1e-4, Rstar-8*sigma), Rstar+8*sigma)
print("spline grid: n=%d  span=%.4g  spacing/sigma=%.4g" % (len(Rg), Rg[-1]-Rg[0], (Rg[1]-Rg[0])/sigma))
d2 = np.diff(fg, 2)
print("fg: range %.6g   monotone=%s   max|2nd diff|/|1st diff| = %.3g"
      % (fg[-1]-fg[0], bool(np.all(np.diff(fg) > 0) or np.all(np.diff(fg) < 0)),
         np.abs(d2).max()/np.abs(np.diff(fg)).mean()))
gvar_th = (eps*p['fbarPhi_prime']*sigma)**2
for label, g in (("spline g", eps*np.interp(R, Rg, fg)),
                 ("linear  g", eps*p['fbarPhi_prime']*(R-Rstar))):
    gd = g - g.mean(); gvar = float((gd*gd).mean())
    print("\n%s   gvar/th = %.4f" % (label, gvar/gvar_th))
    print("   lag t     est/var    exp(-t/tau)   ratio")
    for L in (1, 10, 50, 100, 200, 400, 800):
        c = float((gd[:-L]*gd[L:]).mean()); ex = math.exp(-L*dtr/tau)
        print("   %6.1f   %9.5f  %9.5f   %7.4f" % (L*dtr, c/gvar, ex, (c/gvar)/ex))
    Lc = 1600; cov = np.empty(Lc+1); cov[0] = gvar
    for L in range(1, Lc+1): cov[L] = float((gd[:-L]*gd[L:]).mean())
    zc = np.nonzero(cov <= 0)[0]; cut = int(zc[0]) if len(zc) else len(cov)
    print("   cut %.1f   D/D_cor_th = %.4f   D/(gvar*tau) = %.4f"
          % (cut*dtr, float(np.trapz(cov[:cut], dx=dtr))/p['D_cor'],
             float(np.trapz(cov[:cut], dx=dtr))/(gvar*tau)))
