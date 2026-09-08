"""Step 2: which parameter is the law actually about?

In arXiv:2609.04957 a single eps multiplies h = -theta_dot + f, so eps IS the
linear damping coefficient AND the noise strength (D = eps T by FDT) AND the
bookkeeping prefactor on the escapement f.  The paper reads the answer
Q_inf ~= eps^2/2 as "depends solely on the degree of nonlinearity".

Its own Eqs. 18 and 21 say otherwise:

    D_Phi   ~= eps T / (2 R*^2)      <- the eps here is D = eps T, i.e. FRICTION
    sigma   ~= eps R*^2 / (2 T)      <- the eps here is the friction dissipating
    Q       =  2 D_Phi sigma / Omega^2  ->  eps^2 / 2,  R*^2 CANCELS

R* is the only place the escapement f enters, and it cancels exactly.  So the
prediction for the two-parameter model

    omega_dot = -theta - gamma omega + eps f(theta, omega) + sqrt(2 gamma T) xi

is  Q_inf ~= gamma^2 / 2,  independent of eps  --  i.e. Q_inf = 1/(2 Q_qual^2)
with Q_qual = omega_0/gamma the ordinary mechanical quality factor.

The test: hold gamma fixed and sweep eps over a decade (which moves R* a lot),
then hold eps fixed and sweep gamma.
"""
import sys, time, json
sys.path.insert(0, '.')
import numpy as np
from clock import Grasshopper, VanDerPol, run

T = 0.005
esc = Grasshopper(0.4, 20.0)     # R*^2 = 4 eps theta_r / (pi gamma): a clean eps lever
rows = []

def one(gamma, eps, tag):
    t0 = time.time()
    r = run(esc, gamma=gamma, eps=eps, T=T, n_traj=2000, t_total=800.0,
            dt=2e-3, t_burn=max(10.0 / gamma, 150.0),
            seed=abs(hash((tag, gamma, eps))) % 2**31)
    row = dict(tag=tag, gamma=gamma, eps=eps, Rstar=r['Rstar'], D_phi=r['D_phi'],
               sigma=r['sigma_irr'], sigma_pump=r['sigma_pump'], Omega=r['Omega'],
               Q=r['Q'], Q_gamma=gamma**2/2, Q_eps=eps**2/2)
    rows.append(row)
    print("%-9s gamma=%.4f eps=%.4f  R*=%.3f | Q=%.4e | gamma^2/2=%.4e (x%.2f)"
          "  eps^2/2=%.4e (x%.2f)   [%.0fs]"
          % (tag, gamma, eps, r['Rstar'], r['Q'], gamma**2/2, r['Q']/(gamma**2/2),
             eps**2/2, r['Q']/(eps**2/2), time.time()-t0), flush=True)

print("--- A: gamma fixed at 0.05, eps swept over a factor of 16 ---")
for eps in (0.0125, 0.025, 0.05, 0.1, 0.2):
    one(0.05, eps, "sweep-eps")

print("--- B: eps fixed at 0.05, gamma swept ---")
for gamma in (0.0125, 0.025, 0.05, 0.1, 0.2):
    one(gamma, 0.05, "sweep-gam")

json.dump(rows, open('separate_gamma_eps.json', 'w'), indent=1)
