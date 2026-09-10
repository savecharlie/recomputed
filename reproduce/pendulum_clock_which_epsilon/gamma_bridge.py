"""gamma = 0.10 with the no-filter ruler.  Where does the linear-response result stop?

sigma_nofilter.py: gamma=0.05 -> 0.9806 +- 0.0011 (3 seeds) against 0.9813 predicted, and
gamma=0.20 -> ~0.958 against 0.9315 predicted.  So the isostable prediction is exact at
gamma=0.05 and over-predicts the deficit by 2.7% at gamma=0.20, where gamma/omega_0 is no
longer small and Var(Rbar) = A^2 Var(psi) drops terms it can no longer afford to drop.

PREDICTION, written before running: 0.9634.  If it lands, the result holds through
gamma = 0.1 and fails between 0.1 and 0.2.  If it comes back near 0.975 the failure starts
earlier and the agreement at 0.05 needs re-examining.

Iris (Opus 5), fire 265.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from sigma_nofilter import run

print("gamma = 0.10, u = 3.5.  PREDICTED 0.9634 (isostable), 0.9681 (boxcar, 1 seed)")
rs = []
for sd in (5, 6, 7):
    r = run(3.5, 0.10, seed=sd, t_total=3200.0)
    rs.append(r['ratio'])
    print("   seed=%d   sigma/th = %.4f   tau/th = %.4f" % (sd, r['ratio'], r['tau_ratio']), flush=True)
a = np.array(rs)
print("   MEAN %.4f  s.d. %.4f  s.e.m %.4f" % (a.mean(), a.std(ddof=1), a.std(ddof=1)/math.sqrt(len(a))))
json.dump(dict(gamma=0.10, u=3.5, ratios=rs, mean=float(a.mean()),
               sem=float(a.std(ddof=1)/math.sqrt(len(a))), predicted=0.9634),
          open('gamma_bridge.json', 'w'), indent=1)
