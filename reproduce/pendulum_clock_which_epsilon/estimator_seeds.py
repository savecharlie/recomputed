"""Was 0.945 a bias or an unlucky draw?  Same code path as estimator_check, many seeds, FFT."""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius
import amplitude_phase as ap
from channels import fbarPhi_spline, boxcar
from estimator_scatter import autocov

T, GAM, TR, u = 5e-4, 0.05, 0.4, 3.5
eps = math.pi*GAM*TR*u/4.0
esc = Grasshopper(TR, 20.0); Rstar = limit_cycle_radius(esc, GAM, eps)
p = ap.predict(esc, GAM, eps, T, Rstar)
tau, sigma = 1.0/p['k'], p['sigma_R']
NT, NENS, dtr = 12000, 512, 0.1
LC = int(round(8*tau/dtr))
Rg, fg = fbarPhi_spline(esc, max(1e-4, Rstar-8*sigma), Rstar+8*sigma)
w = max(1, int(round(2.0*np.pi/dtr)))

def est(g):
    gd = g - g.mean()
    cov = autocov(gd, min(LC, gd.shape[0]-2))
    zc = np.nonzero(cov <= 0)[0]
    cut = int(zc[0]) if len(zc) else len(cov)
    return float(np.trapz(cov[:cut], dx=dtr))/p['D_cor'], cut*dtr

rows=[]
for seed in range(1, 25):
    rng = np.random.default_rng(seed)
    a = math.exp(-dtr/tau); s = sigma*math.sqrt(1-a*a)
    R = np.empty((NT, NENS)); R[0] = Rstar + sigma*rng.standard_normal(NENS)
    for i in range(1, NT):
        R[i] = Rstar + a*(R[i-1]-Rstar) + s*rng.standard_normal(NENS)
    r_raw, c_raw = est(eps*np.interp(R, Rg, fg))
    r_box, c_box = est(eps*np.interp(boxcar(R, w), Rg, fg))
    rows.append((seed, r_raw, c_raw, r_box, c_box))
    print("seed %2d   raw %.4f (cut %5.1f)   boxcar %.4f (cut %5.1f)" % rows[-1], flush=True)
a_ = np.array([r[1] for r in rows]); b_ = np.array([r[3] for r in rows])
for nm, v in (("raw", a_), ("boxcar", b_)):
    print("%-7s mean %.4f  s.d. %.4f (%.1f%%)  s.e.m %.4f  min %.3f max %.3f"
          % (nm, v.mean(), v.std(ddof=1), 100*v.std(ddof=1)/v.mean(),
             v.std(ddof=1)/math.sqrt(len(v)), v.min(), v.max()))
json.dump([list(map(float,r)) for r in rows], open('estimator_seeds.json','w'), indent=1)
