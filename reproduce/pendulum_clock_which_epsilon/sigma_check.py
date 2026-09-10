"""Last ruler in the chain: does the boxcar-corrected sigma_R return a KNOWN sigma?

Everything in CHANNELS_FINDINGS now reduces to one number -- sigma_R measured is
0.979 +- 0.003 x sqrt(T/2), uniformly in u, and D_cor's whole shortfall is its square.
Before that becomes physics, point the sigma_R code path at a signal built from a sigma
I typed in myself: slow AR(1) + the measured within-cycle wobble at the ACTUAL cycle
period (6.048, not the 2pi=6.283 the window uses), since an incomplete cycle in the
window is the one thing the analytic factor does not cover.

Iris (Opus 5), fire 264.
"""
import numpy as np, math
from channels import boxcar

tau, dtr, sigma = 20.0, 0.1, 0.015811393924397392
NT, NENS = 12000, 512
w = int(round(2*math.pi/dtr))
Omega = 1.0389866287823892                 # measured deterministic cycle rate at u=3.5
per = 2*math.pi/Omega
bf = math.sqrt(2*tau**2/(w*dtr)**2*(w*dtr/tau - 1 + math.exp(-w*dtr/tau)))
print("window %.3f time units, true cycle period %.3f -> %.4f cycles per window"
      % (w*dtr, per, w*dtr/per))
print("analytic sigma factor %.6f" % bf)

rng = np.random.default_rng(3)
a = math.exp(-dtr/tau); s = sigma*math.sqrt(1-a*a)
slow = np.empty((NT, NENS)); slow[0] = sigma*rng.standard_normal(NENS)
for i in range(1, NT):
    slow[i] = a*slow[i-1] + s*rng.standard_normal(NENS)
t = (np.arange(NT)*dtr)[:, None]
ph = rng.uniform(0, 2*math.pi, NENS)[None, :]      # random cycle phase per trajectory

for name, amp, harm in (("no wobble", 0.0, 1), ("wobble 0.0119 (as measured)", 0.0119*math.sqrt(2), 1),
                        ("wobble, 2nd harmonic too", 0.0119*math.sqrt(2), 2)):
    x = slow.copy()
    if amp:
        x = x + amp*np.cos(Omega*t + ph)
        if harm == 2:
            x = x + 0.4*amp*np.cos(2*Omega*t + 2*ph)
    R = boxcar(x, w)
    sd = float(R.std(ddof=1))
    print("  %-28s raw sd/th %.4f   filtered sd/th %.4f   corrected %.4f"
          % (name, float(x.std(ddof=1))/sigma, sd/sigma, sd/bf/sigma))
