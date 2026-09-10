"""Is the 2.07% sigma_R shortfall the integrator?  Halve dt and look.

sigma_check validated the sigma_R code path on a known sigma (1.002 with a wobble 25%
bigger than the slow fluctuation).  So the shortfall is in the trajectories, and the
remaining suspect is the symplectic-Euler step itself: run_channels uses dt = 2e-3 with
om += (-th - gamma om + eps f) dt + sqrt(2 gamma T dt) z ; th += om dt.
Theory sigma_R^2 = T/2 is a continuum result.

Same scheme, same escapement, no phase bookkeeping -- just the amplitude.

Iris (Opus 5), fire 264.
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius
import amplitude_phase as ap
from channels import boxcar

T, GAM, TR = 5e-4, 0.05, 0.4
esc = Grasshopper(TR, 20.0)


def run(u, dt, n_traj=512, t_total=1200.0, t_burn=400.0, seed=5, T=T):
    eps = math.pi*GAM*TR*u/4.0
    Rstar = limit_cycle_radius(esc, GAM, eps)
    p = ap.predict(esc, GAM, eps, T, Rstar)
    tau = 1.0/p['k']
    rng = np.random.default_rng(seed)
    sig = math.sqrt(2.0*GAM*T*dt)
    th = np.full(n_traj, Rstar); om = np.zeros(n_traj)
    for _ in range(int(round(t_burn/dt))):
        om += (-th - GAM*om + eps*esc.f(th, om))*dt + sig*rng.standard_normal(n_traj)
        th += om*dt
    rec_stride = max(1, int(round(0.1/dt)))
    dtr = rec_stride*dt
    n = int(round(t_total/dt))
    R = np.empty((n//rec_stride, n_traj)); j = 0
    for i in range(n):
        om += (-th - GAM*om + eps*esc.f(th, om))*dt + sig*rng.standard_normal(n_traj)
        th += om*dt
        if (i+1) % rec_stride == 0 and j < R.shape[0]:
            R[j] = np.sqrt(th*th + om*om); j += 1
    R = R[:j]
    w = max(1, int(round(2.0*np.pi/dtr)))
    Rb = boxcar(R, w)
    bf = math.sqrt(2*tau**2/(w*dtr)**2*(w*dtr/tau - 1 + math.exp(-w*dtr/tau)))
    return float(Rb.std(ddof=1))/bf/p['sigma_R'], float(R.std(ddof=1))/p['sigma_R'], dtr


if __name__ == '__main__':
    print("u     dt        dtr     sigma_R corrected / theory   (raw, uncorrected)")
    for u in (3.5, 6.37):
        for dt in (2e-3, 1e-3, 5e-4):
            c, raw, dtr = run(u, dt)
            print("%-5.2f %-9.0e %-7.3f  %.4f                     (%.4f)"
                  % (u, dt, dtr, c, raw), flush=True)
