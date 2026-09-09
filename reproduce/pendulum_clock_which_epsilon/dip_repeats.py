"""How big is the single-point scatter?  Repeats at two u values, new seeds.

The dip sweep has one point (u=3.5) sitting 11% below prediction while its
neighbours sit within 5%.  With one run per point I cannot separate run-to-run
scatter from a systematic over-read of the amplitude-to-phase term where it grows
large -- and the fire-260 comparison points the same way (0.88-0.91 at u=12.7).
So: three fresh seeds at u=3.5, and two at u=1.2 as a control at the other end.

Iris (Opus 5), Sep 9 2026.
"""
import sys, time, json, math
sys.path.insert(0, '.')
import numpy as np
import amplitude_phase as ap
from clock import Grasshopper, run

T, GAM, TR = 5e-4, 0.05, 0.4
esc = Grasshopper(TR, 20.0)
rows = []
for u, seeds in ((3.5, (101, 102, 103)), (1.2, (201, 202))):
    eps = math.pi * GAM * TR * u / 4.0
    for sd in seeds:
        t0 = time.time()
        r = run(esc, gamma=GAM, eps=eps, T=T, n_traj=4000, t_total=1200.0, dt=2e-3,
                t_burn=400.0, seed=sd)
        p = ap.predict(esc, GAM, eps, T, r['Rstar'])
        h = r['D_phi_int'] / p['D_dir']
        rows.append(dict(u=u, seed=sd, h_meas=h, h_pred=p['h'],
                         D_meas=r['D_phi_int'], D_pred=p['D_phi']))
        print("u=%.2f seed=%d  h_meas=%.4f  h_pred=%.4f  ratio=%.4f  [%.0fs]"
              % (u, sd, h, p['h'], h / p['h'], time.time() - t0), flush=True)
        json.dump(rows, open('dip_repeats.json', 'w'), indent=1)
for u in (3.5, 1.2):
    hs = [r['h_meas'] for r in rows if r['u'] == u]
    print("u=%.1f  h = %s   mean %.4f  spread %.4f (%.1f%%)"
          % (u, " ".join("%.4f" % x for x in hs), np.mean(hs), np.std(hs, ddof=1),
             100 * np.std(hs, ddof=1) / np.mean(hs)), flush=True)
