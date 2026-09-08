"""Step 3: is the residual eps-dependence amplitude-to-phase feedthrough?

Sweep A and B agree that Q_inf/(gamma^2/2) is a function of eps/gamma alone, that
sigma_irr matches gamma R*^2/(2T) at every point (1-8%), and that the whole excess
sits in D_Phi.  That fingers ONE mechanism: the paper replaces the fluctuating
radius R by its deterministic R* inside

    Phi_dot = -1 + eps fbar_Phi(R) + noise,

and drops the resulting slow frequency wander.  If that is the cause, the excess
must vanish for an escapement whose fbar_Phi is identically zero -- because then the
phase velocity does not know the amplitude at all.

Van der Pol is exactly that escapement:

    fbar_Phi(R) = (1/2pi) int cos(Phi) (a - R^2 cos^2 Phi)(R sin Phi)/R dPhi = 0

by parity, every term odd in sin(Phi)cos(Phi).  And the noiseless run agrees:
Omega = -0.999935 for van der Pol against -1.0364 for the grasshopper.

Prediction: with van der Pol, Q_inf = gamma^2/2 flat across the SAME eps/gamma range
(0.75 to 4) over which the grasshopper climbed to 2.7x.
"""
import sys, time, json
sys.path.insert(0, '.')
import numpy as np
from clock import VanDerPol, run

T = 0.005
esc = VanDerPol(2.0)          # R*^2 = 8(1 - gamma/(2 eps))
gamma = 0.05
rows = []
for eps in (0.0375, 0.05, 0.075, 0.1, 0.15, 0.2):
    t0 = time.time()
    r = run(esc, gamma=gamma, eps=eps, T=T, n_traj=2000, t_total=800.0, dt=2e-3,
            t_burn=max(10.0/gamma, 150.0), seed=abs(hash(("vdp", eps))) % 2**31)
    Dp = gamma*T/(2*r['Rstar']**2); sp = gamma*r['Rstar']**2/(2*T)
    rows.append(dict(gamma=gamma, eps=eps, Rstar=r['Rstar'], D_phi=r['D_phi'],
                     D_pred=Dp, sigma=r['sigma_irr'], sigma_pred=sp,
                     Omega=r['Omega'], Q=r['Q'], Q_gamma=gamma**2/2))
    print("vdp  eps/gamma=%4.1f  R*=%.3f | D/Dpred %.2f | sigma/pred %.2f | Omega %.4f"
          " | Q/(gamma^2/2) = %.3f   [%.0fs]"
          % (eps/gamma, r['Rstar'], r['D_phi']/Dp, r['sigma_irr']/sp, r['Omega'],
             r['Q']/(gamma**2/2), time.time()-t0), flush=True)
json.dump(rows, open('vdp_control.json','w'), indent=1)
