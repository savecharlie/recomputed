"""What the noiseless limit cycle says about the two constants the averaging assumes.

The averaged amplitude equation is  dx = -k x dt + sqrt(gamma T) dW, whose stationary
variance is gamma T/(2k).  Both constants come from a phase average taken as if Phi
advanced uniformly at rate 1:

    diffusion   D_R = 2 gamma T <sin^2 Phi>      assumed <sin^2 Phi> = 1/2
    drift rate  k   = gamma/2 - eps d(fbar_R)/dR  with fbar_R = <sin Phi f>_Phi

On the ACTUAL orbit Phi does not advance uniformly:
    dPhi/dt = -1 - (gamma/2) sin 2Phi + eps f cos Phi / R
so a TIME average and a PHASE average are different things.  This measures both
on the noiseless cycle, exactly, by integrating the real orbit.

PREDICTION, written before running (fire 265):  <sin^2 Phi>_time - 1/2 vanishes at
first order in gamma and in eps (the correcting weight is odd against sin^2), so this
should come back within ~1e-3 of 0.5 and NOT explain a 2% sigma_R deficit.  If it
comes back near 0.48 the deficit is solved.

Iris (Opus 5), fire 265.
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius
import amplitude_phase as ap

THETA_R, N = 0.4, 20.0
GAMMA, T = 0.05, 5e-4


def orbit(esc, gamma, eps, Rstar, dt=1e-5, n_cyc=40):
    """Integrate the noiseless system onto its cycle, return one clean period."""
    th, om = Rstar, 0.0
    steps = int(round(n_cyc * 2 * math.pi / dt))
    # burn onto the cycle
    for _ in range(steps):
        om += (-th - gamma * om + eps * esc.f(th, om)) * dt
        th += om * dt
    # now record until theta crosses zero upward twice
    ths, oms = [], []
    prev = th
    crossings = 0
    for _ in range(steps):
        om += (-th - gamma * om + eps * esc.f(th, om)) * dt
        th += om * dt
        if prev < 0.0 <= th and om > 0:
            crossings += 1
            if crossings == 1:
                ths, oms = [], []
            elif crossings == 2:
                break
        prev = th
        ths.append(th); oms.append(om)
    return np.array(ths), np.array(oms), dt


for u in (3.5, 6.37):
    eps = math.pi * GAMMA * THETA_R * u / 4.0
    esc = Grasshopper(THETA_R, N)
    Rstar = limit_cycle_radius(esc, GAMMA, eps)
    th, om, dt = orbit(esc, GAMMA, eps, Rstar)
    R = np.sqrt(th * th + om * om)
    Phi = np.arctan2(om, th)
    P = len(th) * dt
    s2 = float((np.sin(Phi) ** 2).mean())          # TIME average over one period
    # the same average weighted by R^2, which is what enters d(R^2)/dt
    print("u = %.2f  eps = %.5f  R* = %.5f  period = %.6f (2pi = %.6f)"
          % (u, eps, Rstar, P, 2 * math.pi))
    print("   <sin^2 Phi>_time      = %.6f   (assumed 0.5, ratio %.5f)" % (s2, s2 / 0.5))
    print("   <R>_time/R*           = %.6f" % (R.mean() / Rstar))
    print("   R wobble sd           = %.6f   (sigma_R th = %.6f)"
          % (R.std(), math.sqrt(T / 2)))
    # k, both ways
    p = ap.predict(esc, GAMMA, eps, T, Rstar)
    print("   k/gamma (phase avg)   = %.6f" % (p['k'] / GAMMA))
    print("   sigma_R th            = %.8f" % p['sigma_R'])
