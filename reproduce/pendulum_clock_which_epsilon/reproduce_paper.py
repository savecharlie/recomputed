"""Step 1: reproduce the paper's own claim, gamma == eps, Q_inf ~= eps^2/2.

If this does not come out, nothing later in this directory means anything.
"""
import sys, time, json
sys.path.insert(0, '.')
import numpy as np
from clock import VanDerPol, Grasshopper, Graham, run

T = 0.005
rows = []
for esc in (VanDerPol(2.0), Grasshopper(0.4, 20.0), Graham(0.1, 20.0)):
    for eps in (0.02, 0.05, 0.1, 0.2):
        t0 = time.time()
        r = run(esc, gamma=eps, eps=eps, T=T, n_traj=2000, t_total=1200.0,
                dt=1e-3, seed=hash((esc.name, eps)) % 2**31)
        pred_Q = eps ** 2 / 2
        pred_D = eps * T / (2 * r['Rstar'] ** 2)
        pred_s = eps * r['Rstar'] ** 2 / (2 * T)
        rows.append(dict(model=esc.name, eps=eps, Rstar=r['Rstar'],
                         D_phi=r['D_phi'], D_pred=pred_D,
                         sigma=r['sigma_irr'], sigma_pred=pred_s,
                         sigma_pump=r['sigma_pump'],
                         Omega=r['Omega'], Q=r['Q'], Q_pred=pred_Q))
        print("%-22s eps=%.3f R*=%.4f | D_phi %.4e (pred %.4e) | sigma %.4g (pred %.4g)"
              " | Q %.4e (pred %.4e)  ratio %.3f   [%.0fs]"
              % (esc.name, eps, r['Rstar'], r['D_phi'], pred_D, r['sigma_irr'],
                 pred_s, r['Q'], pred_Q, r['Q'] / pred_Q, time.time() - t0), flush=True)
json.dump(rows, open('reproduce_paper.json', 'w'), indent=1)
