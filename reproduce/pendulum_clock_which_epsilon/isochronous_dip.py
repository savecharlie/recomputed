"""THE SHARP TEST of the amplitude-to-phase theory: a non-monotonic dip, and the
prediction that its minimum sits exactly at R* = sqrt(2) theta_r.

amplitude_phase.py derives, with no fitted parameters,

    D_Phi = (gamma T/2)[1/R*^2 + (eps fbar_Phi'(R*)/k)^2]
    h = D_Phi / D_dir = 1 + (2-u)^2/(4(u-1))   for the sharp grasshopper,  u=(R*/theta_r)^2

h is NOT monotonic.  It equals ONE exactly at u = 2, because that is where the
cycle-averaged frequency is stationary in amplitude -- the escapement is isochronous
at its own running amplitude -- and rises on both sides.  Finding #6 measured h
rising from 1.0 to 2.7 and read it as monotonic in eps/gamma; it only looked
monotonic because every point sat at u >= 1.6.

Nothing in the paper, and nothing in finding #6, predicts a minimum.  If the
measured h dips to 1 at u=2 and comes back up below it, the mechanism is settled.

T is dropped to 5e-4 (from 5e-3) so the thermal amplitude width sigma_R = sqrt(T/2)
is small against the distance to the escapement's critical angle -- at 5e-3 the
Gaussian readout smears the low-u side flat.  D_Phi scales with T, so the relative
precision of the estimator is unchanged; only the prediction gets sharper.

Iris (Opus 5), Sep 9 2026.
"""
import sys, time, json, math
sys.path.insert(0, '.')
import numpy as np
import amplitude_phase as ap
from clock import Grasshopper, run

T = 5e-4
GAM = 0.05
TR = 0.4
esc = Grasshopper(TR, 20.0)
US = (1.2, 1.5, 2.0, 2.6, 3.5, 5.0)
TTOT, NTRAJ, DT = 1200.0, 4000, 2e-3

rows = []
print("u -> eps = pi gamma theta_r u / 4;  gamma=%.3f theta_r=%.2f T=%.1e "
      "n_traj=%d t_total=%.0f" % (GAM, TR, T, NTRAJ, TTOT), flush=True)
print("%6s %8s %7s | %10s %10s %10s | %7s %7s %7s"
      % ("u", "eps", "R*", "D_dir", "D_meas", "D_pred", "h_meas", "h_pred", "h_sharp"), flush=True)
for u in US:
    eps = math.pi * GAM * TR * u / 4.0
    t0 = time.time()
    r = run(esc, gamma=GAM, eps=eps, T=T, n_traj=NTRAJ, t_total=TTOT, dt=DT,
            t_burn=400.0, seed=abs(hash(("dip", u))) % 2 ** 31)
    p = ap.predict(esc, GAM, eps, T, r['Rstar'])
    Dm = r['D_phi_int']
    row = dict(u=u, eps=eps, Rstar=r['Rstar'], D_dir=p['D_dir'], D_cor=p['D_cor'],
               D_pred=p['D_phi'], D_meas=Dm, D_meas_origin=r['D_phi'],
               h_meas=Dm / p['D_dir'], h_pred=p['h'], h_sharp=float(ap.gh_h(u)),
               Omega=r['Omega'], sigma_irr=r['sigma_irr'], Q=r['Q'],
               sigma_R=p['sigma_R'], k=p['k'], secs=time.time() - t0)
    rows.append(row)
    print("%6.2f %8.4f %7.3f | %.4e %.4e %.4e | %7.3f %7.3f %7.3f  [%.0fs]"
          % (u, eps, r['Rstar'], p['D_dir'], Dm, p['D_phi'],
             row['h_meas'], row['h_pred'], row['h_sharp'], row['secs']), flush=True)
    json.dump(rows, open('isochronous_dip.json', 'w'), indent=1)

hm = [r['h_meas'] for r in rows]
i = int(np.argmin(hm))
print("\nmeasured minimum of h at u = %.2f  (theory: exactly 2)" % rows[i]['u'], flush=True)
print("h_meas: " + "  ".join("%.3f" % x for x in hm), flush=True)
