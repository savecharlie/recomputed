"""The last honest gap in the amplitude-to-phase prediction, and where it goes.

Against the fire-260 sweep, the full prediction (Hermite readout + the known
through-origin fit bias) lands at 0.98 +/- 0.08 on eight of ten points.  The two
misses are both SUB-CRITICAL: R* = 0.318 < theta_r = 0.4, where the theory says
there is no amplitude-to-phase channel at all (fbar_Phi == 0 there, exactly) and
D_Phi should be pure direct diffusion.  Measured, they sit 10-14% high.

The explanation is Jensen, and the reason it is LARGE there rather than negligible
is a structural fact worth writing down: below the critical angle the escapement's
phase-averaged force is the CONSTANT 2/pi, so its slope is zero and the amplitude
relaxation rate is k = gamma/2 rather than k = gamma.  The amplitude is held half
as tightly, sigma_R^2 = gamma T/(2k) = T, and sigma_R/R* reaches 0.22.

D_dir is (gamma T/2) <1/R^2>, and <1/R^2> exceeds 1/R*^2 by about 3 sigma^2/R*^2.

    gamma=0.05  eps=0.0125  R*=0.318   k/gamma=0.53   sigma_R/R*=0.216
      D_dir at R*        1.2337e-03      measured/pred  1.144   (10-14% high)
      D_dir with <1/R^2> 1.4796e-03      measured/pred  0.970

So the excess is real, has the right sign and roughly the right size.  The Gaussian
estimate OVER-corrects (predicting ~20% where 10-14% is observed) because the true
stationary amplitude distribution is repelled from the origin and is not Gaussian
out in the tail that dominates <1/R^2>.  Called explained in direction and
magnitude, not closed.

Iris (Opus 5), Sep 9 2026.
"""
import sys, json
sys.path.insert(0, '.')
import numpy as np
import amplitude_phase as ap
from clock import Grasshopper

esc, T, TR = Grasshopper(0.4, 20.0), 0.005, 0.4
x, w = np.polynomial.hermite_e.hermegauss(120)
w = w / w.sum()
print("%8s %8s %7s %8s %8s | %11s %11s %11s %7s"
      % ("gamma", "eps", "R*", "k/gamma", "sig/R*", "D_dir(R*)", "D_dir<1/R2>",
         "D_meas", "ratio"))
for row in json.load(open('separate_gamma_eps.json')):
    g, e, R = row['gamma'], row['eps'], row['Rstar']
    if R > TR:
        continue
    p = ap.predict(esc, g, e, T, R, hermite=False)
    sig = p['sigma_R']
    Rs = R + sig * x
    inv = float(np.sum(w * np.where(Rs > 1e-3, 1.0 / np.maximum(Rs, 1e-3) ** 2, 0.0)))
    D0, D1 = g * T / (2 * R * R), g * T / 2.0 * inv
    print("%8.4f %8.4f %7.3f %8.3f %8.3f | %.5e %.5e %.5e %7.3f"
          % (g, e, R, p['k'] / g, sig / R, D0, D1, row['D_phi'], row['D_phi'] / D1))
print("\nanalytic 3 sigma^2/R*^2 = %.1f%%; Gauss-Hermite <1/R^2> is %.1f%% above 1/R*^2"
      % (300 * sig * sig / (R * R), 100 * (inv * R * R - 1)))
