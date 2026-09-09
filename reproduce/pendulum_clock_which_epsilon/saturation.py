"""The other face of the same formula: D_Phi STOPS FALLING with amplitude.

Paper's Eq. 18 (and every textbook phase-reduction) gives D_Phi = D/(2 R*^2): swing
harder, diffuse less.  Adding the amplitude-to-phase channel, for the grasshopper,

    D_Phi = (gamma T/2)[1/R*^2 + (u/4)/R*^2 * ...]  ->  gamma T / (8 theta_r^2)

a FLOOR set only by the damping, the temperature and the escapement's critical
angle -- the amplitude cancels out of it.  Past the isochronous amplitude a bigger
swing buys nothing, and the uncertainty product Q_inf -> gamma^2 R*^2/(8 theta_r^2)
gets strictly WORSE.

Lever: hold eps and DROP gamma.  u = 4 eps/(pi gamma theta_r) grows, R*^2 = theta_r^2 u
grows, and the frequency shift eps fbar_Phi ~ 2 eps/(pi R*) SHRINKS, so the run stays
inside weak nonlinearity as the test gets stronger.  Naive theory says D_Phi/gamma^2
is flat; the floor says D_Phi/gamma is flat.  Over this range those differ by 8x.

Measurement window must be long against the amplitude correlation time 1/gamma or the
slope has not reached its asymptote -- t_total and t_burn are both scaled as 1/gamma.

Iris (Opus 5), Sep 9 2026.
"""
import sys, time, json
sys.path.insert(0, '.')
import numpy as np
import amplitude_phase as ap
from clock import Grasshopper, run

T = 5e-4
EPS = 0.05
TR = 0.4
esc = Grasshopper(TR, 20.0)
GAMS = (0.05, 0.025, 0.0125, 0.00625)

rows = []
print("eps=%.3f theta_r=%.2f T=%.1e ; floor = gamma T/(8 theta_r^2) = gamma * %.4e"
      % (EPS, TR, T, T / (8 * TR * TR)), flush=True)
print("%8s %6s %7s | %10s %10s %10s | %8s %8s"
      % ("gamma", "u", "R*", "D_dir", "D_meas", "D_pred", "D/gam", "D/gam^2"), flush=True)
for g in GAMS:
    tt = max(1200.0, 20.0 / g)
    t0 = time.time()
    r = run(esc, gamma=g, eps=EPS, T=T, n_traj=4000, t_total=tt, dt=2e-3,
            t_burn=max(400.0, 20.0 / g), seed=abs(hash(("sat", g))) % 2 ** 31)
    p = ap.predict(esc, g, EPS, T, r['Rstar'])
    Dm = r['D_phi_int']
    row = dict(gamma=g, u=float(ap.gh_u(g, EPS, TR)), Rstar=r['Rstar'],
               D_dir=p['D_dir'], D_cor=p['D_cor'], D_pred=p['D_phi'], D_meas=Dm,
               D_meas_origin=r['D_phi'], Omega=r['Omega'], Q=r['Q'],
               sigma_irr=r['sigma_irr'], t_total=tt, secs=time.time() - t0)
    rows.append(row)
    print("%8.5f %6.2f %7.3f | %.4e %.4e %.4e | %.3e %.3e  [%.0fs]"
          % (g, row['u'], r['Rstar'], p['D_dir'], Dm, p['D_phi'],
             Dm / g, Dm / g ** 2, row['secs']), flush=True)
    json.dump(rows, open('saturation.json', 'w'), indent=1)

print("\nD_meas/gamma over the sweep (flat == the floor is real):", flush=True)
print("  " + "  ".join("%.3e" % (r['D_meas'] / r['gamma']) for r in rows), flush=True)
print("D_meas/gamma^2 (flat == the paper's 1/R*^2 only):", flush=True)
print("  " + "  ".join("%.3e" % (r['D_meas'] / r['gamma'] ** 2) for r in rows), flush=True)
