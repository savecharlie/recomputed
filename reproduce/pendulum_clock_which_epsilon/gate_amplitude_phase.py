"""GATE the amplitude-to-phase theory against things already known, before it
is allowed to say anything new.  Six checks; any red one stops the fire.

Iris (Opus 5), Sep 9 2026.
"""
import sys, math
sys.path.insert(0, '.')
import numpy as np
import amplitude_phase as ap
from clock import Grasshopper, VanDerPol, Graham, limit_cycle_radius

fails = []
def chk(name, got, want, tol, rel=True):
    err = abs(got - want) / (abs(want) if (rel and want) else 1.0)
    ok = err <= tol
    print(("  %-62s %12.6g vs %12.6g   %s (%.2g)"
           % (name, got, want, "PASS" if ok else "FAIL", err)))
    if not ok:
        fails.append(name)

# --- 1. sharp-sgn quadrature vs my closed forms -----------------------------
print("1. closed-form phase averages vs quadrature on the SHARP escapement")
class SharpGH:
    def __init__(s, tr): s.theta_r = tr
    def f(s, th, om): return -np.sign(th - s.theta_r * np.sign(om))
for tr in (0.4, 0.1):
    for R in (0.25 * tr, 0.9 * tr, 1.01 * tr, 1.5 * tr, 4.0 * tr, 20.0 * tr):
        q = ap.phase_averages(SharpGH(tr), R, P=2000001)
        chk("fbar_R  theta_r=%.2f R/theta_r=%5.2f" % (tr, R / tr), q[0],
            float(ap.gh_fbar_R(R, tr)), 2e-3)
        chk("fbar_Phi theta_r=%.2f R/theta_r=%5.2f" % (tr, R / tr), q[1],
            float(ap.gh_fbar_Phi(R, tr)), 3e-3, rel=(R > 1.02 * tr))

# --- 2. van der Pol is isochronous at every amplitude (parity) --------------
print("2. van der Pol fbar_Phi == 0 for all R  (why finding #6's control was flat)")
for R in (0.3, 1.0, 2.0, 5.0):
    chk("vdP fbar_Phi(R=%.1f)" % R, ap.phase_averages(VanDerPol(2.0), R)[1], 0.0,
        1e-12, rel=False)
print("   ... and its fbar_R matches the corrected Table 1 entry  R - R^3/8:")
for R in (0.5, 2.0, 3.0):
    chk("vdP fbar_R(R=%.1f)" % R, ap.phase_averages(VanDerPol(2.0), R)[0],
        R - R ** 3 / 8.0, 1e-6)

# --- 3. the limit cycle and u agree with the independent bisection ----------
print("3. R*^2 = theta_r^2 u  against clock.py's bisection on the smoothed force")
for gamma, eps in ((0.05, 0.05), (0.05, 0.2), (0.0125, 0.05), (0.00625, 0.05)):
    R = limit_cycle_radius(Grasshopper(0.4, 20.0), gamma, eps)
    chk("R* gamma=%.5f eps=%.3f" % (gamma, eps), R,
        0.4 * math.sqrt(ap.gh_u(gamma, eps, 0.4)), 6e-3)

# --- 4. k = gamma exactly for the grasshopper -------------------------------
print("4. amplitude relaxation rate k == gamma exactly (fbar_R goes as 1/R)")
for gamma, eps in ((0.05, 0.05), (0.05, 0.2), (0.0125, 0.05)):
    R = limit_cycle_radius(Grasshopper(0.4, 20.0), gamma, eps)
    p = ap.predict(Grasshopper(0.4, 20.0), gamma, eps, 0.005, R, hermite=False)
    chk("k/gamma gamma=%.4f eps=%.3f" % (gamma, eps), p['k'] / gamma, 1.0, 0.03)

# --- 5. Hermite sum collapses onto the b_1 form as sigma_R -> 0 -------------
print("5. Hermite readout -> linear-response form as T -> 0 (same mechanism, two maths)")
esc = Grasshopper(0.4, 20.0)
R = limit_cycle_radius(esc, 0.05, 0.05)
for T in (5e-3, 5e-5):
    a = ap.predict(esc, 0.05, 0.05, T, R, hermite=True)
    b = ap.predict(esc, 0.05, 0.05, T, R, hermite=False)
    print("   T=%-8g  D_cor hermite=%.4e  linear=%.4e  ratio=%.4f  sigma_R/R*=%.3f"
          % (T, a['D_cor'], b['D_cor'], a['D_cor'] / b['D_cor'], a['sigma_R'] / R))
chk("hermite/linear at T=5e-5 (sigma_R tiny)",
    ap.predict(esc, 0.05, 0.05, 5e-5, R)['D_cor']
    / ap.predict(esc, 0.05, 0.05, 5e-5, R, hermite=False)['D_cor'], 1.0, 0.05)

# --- 6. the fit-bias formula, against a synthetic Var(J) it did not see -----
print("6. fit_bias_fraction vs a synthetic through-origin fit (instrument check)")
for tau, tmax in ((20.0, 800.0), (80.0, 800.0), (5.0, 400.0)):
    D_dir, D_cor = 1.0, 3.0
    ts = np.linspace(tmax / 40, tmax, 40)
    var = 2 * D_dir * ts + 2 * D_cor * (ts - tau * (1 - np.exp(-ts / tau)))
    m = ts > tmax * 0.25
    slope = float(np.sum(ts[m] * var[m]) / np.sum(ts[m] ** 2) / 2.0)
    pred = D_dir + D_cor * (1.0 - ap.fit_bias_fraction(tau, tmax) / tau * 1.0)
    # fit_bias_fraction returns tau*<t>/<t^2>; deficit on D_cor is exactly that
    pred = D_dir + D_cor * (1.0 - ap.fit_bias_fraction(tau, tmax))
    chk("through-origin slope tau=%.0f tmax=%.0f" % (tau, tmax), slope, pred, 0.02)
    # and an intercept fit should recover the truth
    A = np.vstack([ts[m], np.ones(m.sum())]).T
    sl, ic = np.linalg.lstsq(A, var[m], rcond=None)[0]
    chk("  intercept fit recovers D_dir+D_cor", sl / 2.0, D_dir + D_cor, 0.01)

print()
if fails:
    print("RED: %d gate(s) failed:" % len(fails))
    for f in fails: print("   -", f)
    sys.exit(1)
print("ALL GATES GREEN — the theory module may now say something new.")
