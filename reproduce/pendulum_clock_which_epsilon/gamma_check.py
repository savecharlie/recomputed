"""PREDICTION FIRST: the 2.07% sigma_R shortfall is the next order in stochastic averaging,
so it is proportional to gamma and to nothing else.

The reasoning.  theory sigma_R^2 = T/2 comes from projecting the omega-only noise onto the
amplitude and averaging sin^2(Phi) -> 1/2 over a cycle.  That average is legitimate only
because the amplitude relaxes at rate gamma = 0.05 while the phase turns at rate omega_0 = 1,
i.e. it is the leading term of an expansion in gamma/omega_0 = 0.05.  A 2.07% residual is
exactly the size of the FIRST correction to that expansion -- and, crucially, such a
correction depends on gamma/omega_0 and on NOTHING ELSE, which is what the sweep already
says: the shortfall did not move with the drive strength u (0.9793 at 3.5, 0.9774 at 6.37),
with the temperature, or with the timestep over a factor of four.

  d(gamma) := 1 - sigma_R_measured / sqrt(T/2)

  H1  d is O(gamma)         ->  gamma = 0.05, 0.1, 0.2  gives  2.1%,  4.2%,  8.4%
  H2  d is O(gamma^2)       ->                                 2.1%,  8.3%, 33%
  H3  d is gamma-independent ->                                2.1%,  2.1%,  2.1%   (hypothesis dead)

H3 would mean the shortfall is not the averaging and I am back to not knowing.  Written down
BEFORE the run, because a prediction stated afterwards is a rationalisation.

Iris (Opus 5), fire 264.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius
import amplitude_phase as ap
from channels import boxcar

T, TR, U = 5e-4, 0.4, 3.5
esc = Grasshopper(TR, 20.0)


def run(gam, dt=2e-3, n_traj=512, t_total=1200.0, seed=5):
    eps = math.pi*gam*TR*U/4.0
    Rstar = limit_cycle_radius(esc, gam, eps)
    p = ap.predict(esc, gam, eps, T, Rstar)
    tau = 1.0/p['k']
    rng = np.random.default_rng(seed)
    sig = math.sqrt(2.0*gam*T*dt)
    th = np.full(n_traj, Rstar); om = np.zeros(n_traj)
    for _ in range(int(round(max(400.0, 20*tau)/dt))):
        om += (-th - gam*om + eps*esc.f(th, om))*dt + sig*rng.standard_normal(n_traj)
        th += om*dt
    stride = max(1, int(round(0.1/dt))); dtr = stride*dt
    n = int(round(t_total/dt))
    R = np.empty((n//stride, n_traj)); j = 0
    for i in range(n):
        om += (-th - gam*om + eps*esc.f(th, om))*dt + sig*rng.standard_normal(n_traj)
        th += om*dt
        if (i+1) % stride == 0 and j < R.shape[0]:
            R[j] = np.sqrt(th*th + om*om); j += 1
    R = R[:j]
    w = max(1, int(round(2.0*np.pi/dtr)))
    Rb = boxcar(R, w)
    bf = math.sqrt(2*tau**2/(w*dtr)**2*(w*dtr/tau - 1 + math.exp(-w*dtr/tau)))
    return float(Rb.std(ddof=1))/bf/p['sigma_R'], tau, p['sigma_R'], Rstar


if __name__ == '__main__':
    print("gamma   tau     R*      sigma_R/th (2 seeds)      d = 1-ratio    H1 pred   H2 pred")
    rows = []
    for gam in (0.05, 0.1, 0.2):
        v = []
        for sd in (5, 6):
            r, tau, s_th, Rstar = run(gam, seed=sd)
            v.append(r)
        m = float(np.mean(v)); d = 1.0 - m
        print("%-7.3f %-7.2f %-7.4f %.4f %.4f -> %.4f   %6.2f%%      %5.2f%%    %5.2f%%"
              % (gam, tau, Rstar, v[0], v[1], m, 100*d, 100*0.0207*(gam/0.05),
                 100*0.0207*(gam/0.05)**2), flush=True)
        rows.append(dict(gamma=gam, ratios=v, mean=m, d=d, tau=tau, Rstar=Rstar))
    json.dump(rows, open('gamma_check.json', 'w'), indent=1)
