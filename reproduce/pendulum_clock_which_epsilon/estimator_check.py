"""Is the 5% D_cor deficit physics, or is it my estimator?  Ask a process whose answer I know.

CHANNELS_FINDINGS (fire 263) measured D_cor from its definition and got a uniform
0.93-0.96 x theory at every u.  A uniform systematic across every condition is the
signature of a ruler, so before that becomes a finding the estimator gets pointed at a
process with a KNOWN D_cor.

The construction is exact, not approximate.  The Hermite closed form
    D_cor = eps^2/k * sum_{n>=1} b_n^2 n!/n,   b_n = E[fbar(R*+sigma Z) He_n(Z)]/n!
is derived for R a Gaussian OU process of variance sigma^2 and correlation time tau=1/k
read out through the nonlinear map g = eps*fbar_Phi(R).  So: simulate exactly that OU
process, with the SAME sigma, tau, dtr, record length and ensemble size as the pendulum
run, push it through the SAME spline g, and run the SAME estimator.  Any departure from
1.000 is the estimator, because the physics is now a closed form I typed in myself.

Three arms, to localise it:
  raw        OU straight into the estimator      -> truncation + finite sample only
  boxcar     OU smoothed over 2pi first          -> adds what the cycle-average does
  boxcar-lin same, with g linearised at R*       -> boxcar x truncation, no nonlinearity

Iris (Opus 5), fire 264, Sep 10 2026.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius
import amplitude_phase as ap
from channels import boxcar, fbarPhi_spline

T, GAM, TR = 5e-4, 0.05, 0.4
DTR, NT, NENS, TTOT = 0.1, 12000, 512, 1200.0


def d_cor_estimate(g, dtr, tau_th, cut_rule='zero'):
    """The estimator exactly as channels.analyse uses it."""
    gd = g - g.mean()
    gvar = float((gd * gd).mean())
    Lc = min(int(round(8.0 * tau_th / dtr)), gd.shape[0] - 2)
    cov = [gvar]
    for L in range(1, Lc + 1):
        cov.append(float((gd[:-L] * gd[L:]).mean()))
    cov = np.array(cov)
    zc = np.nonzero(cov <= 0)[0]
    cut = int(zc[0]) if len(zc) else len(cov)
    return float(np.trapz(cov[:cut], dx=dtr)), gvar, cut * dtr


def one(u, seed=7):
    eps = math.pi * GAM * TR * u / 4.0
    esc = Grasshopper(TR, 20.0)
    Rstar = limit_cycle_radius(esc, GAM, eps)
    p = ap.predict(esc, GAM, eps, T, Rstar)
    tau, sigma = 1.0 / p['k'], p['sigma_R']

    # exact discrete OU: R_{n+1} = R* + a (R_n - R*) + sqrt(sigma^2 (1-a^2)) xi
    rng = np.random.default_rng(seed)
    a = math.exp(-DTR / tau)
    s = sigma * math.sqrt(1.0 - a * a)
    R = np.empty((NT, NENS))
    R[0] = Rstar + sigma * rng.standard_normal(NENS)
    for i in range(1, NT):
        R[i] = Rstar + a * (R[i - 1] - Rstar) + s * rng.standard_normal(NENS)

    Rg, fg = fbarPhi_spline(esc, max(1e-4, Rstar - 8 * sigma), Rstar + 8 * sigma)
    w = max(1, int(round(2.0 * np.pi / DTR)))
    Rb = boxcar(R, w)
    bf = math.sqrt(2.0 * tau ** 2 / (w * DTR) ** 2 *
                   (w * DTR / tau - 1.0 + math.exp(-w * DTR / tau)))

    rows = {}
    for name, Ruse in (('raw', R), ('boxcar', Rb)):
        g = eps * np.interp(Ruse, Rg, fg)
        D, gv, cut = d_cor_estimate(g, DTR, tau)
        rows[name] = dict(D=D, ratio=D / p['D_cor'], gvar=gv, cut=cut,
                          tau_g=D / gv, sigma=float(Ruse.std(ddof=1)))
    g = eps * p['fbarPhi_prime'] * (Rb - Rb.mean())          # linearised readout
    D, gv, cut = d_cor_estimate(g, DTR, tau)
    D_cor_lin = (eps * p['fbarPhi_prime']) ** 2 * sigma ** 2 * tau
    rows['boxcar-lin'] = dict(D=D, ratio=D / D_cor_lin, gvar=gv, cut=cut,
                              tau_g=D / gv, sigma=float(Rb.std(ddof=1)))
    rows['meta'] = dict(u=u, eps=eps, Rstar=Rstar, tau=tau, sigma_th=sigma,
                        D_cor_th=p['D_cor'], D_cor_lin_th=D_cor_lin,
                        hermite_over_linear=p['D_cor'] / D_cor_lin,
                        boxcar_sigma_factor=bf, seed=seed)
    return rows


if __name__ == '__main__':
    out = []
    for u in [float(x) for x in sys.argv[1:]] or [3.5, 6.37]:
        r = one(u)
        m = r['meta']
        print("u = %.2f   tau = %.3f  sigma_R = %.5f   D_cor_th = %.5e "
              "(hermite/linear = %.4f)" % (u, m['tau'], m['sigma_th'],
                                           m['D_cor_th'], m['hermite_over_linear']))
        for k in ('raw', 'boxcar', 'boxcar-lin'):
            d = r[k]
            print("   %-11s ratio %7.4f   sigma_meas/th %6.4f   tau_g %6.2f   cut %6.1f"
                  % (k, d['ratio'], d['sigma'] / m['sigma_th'], d['tau_g'], d['cut']))
        print("   boxcar sigma factor predicted %.4f" % m['boxcar_sigma_factor'])
        out.append(r)
        sys.stdout.flush()
    json.dump(out, open('estimator_check.json', 'w'), indent=1)
