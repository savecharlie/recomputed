"""Where the residual of finding #6 comes from: amplitude-to-phase conversion.

Finding #6 showed that Q_inf for a stochastic pendulum clock is gamma^2/2 -- the
LINEAR DAMPING, not the degree of nonlinearity -- with a residual h = Q/(gamma^2/2)
that collapsed on eps/gamma alone and lived entirely in D_Phi.  That residual was
left as an open loose end: a collapse variable and a mechanism, no derivation.

This file derives it.

Polar reduction of  theta'' = -theta - gamma theta' + eps f(theta,theta')
+ sqrt(2 gamma T) xi,  with theta = R cos Phi,  theta' = R sin Phi:

    R'   = -gamma R sin^2(Phi) + sin(Phi) (eps f + xi)
    Phi' = -1 - gamma sin(Phi) cos(Phi) + (cos(Phi)/R) (eps f + xi)

Cycle-averaging (valid while gamma, eps << 1, so amplitude relaxation is slow
against the period):

    <R'>   = -gamma R / 2 + eps fbar_R(R),      fbar_R = <sin(Phi) f>
    <Phi'> = -1 + eps fbar_Phi(R),              fbar_Phi = <cos(Phi) f> / R

So the INSTANTANEOUS FREQUENCY DEPENDS ON AMPLITUDE, and the amplitude is a noisy
relaxing variable.  Two independent channels feed phase diffusion:

  DIRECT   the cos(Phi) xi / R term            ->  D_dir = gamma T / (2 R*^2)
  AM->PM   amplitude noise read out as rate    ->  D_cor = int_0^inf Cov(eps fbar_Phi(R)) dt

The second is the one the paper's R -> R* step throws away.  With the amplitude a
Gaussian Ornstein-Uhlenbeck variable of relaxation rate k and width sigma_R,

    k = gamma/2 - eps fbar_R'(R*),      sigma_R^2 = gamma T / (2 k)

and, EXACTLY for Gaussian OU driving a nonlinear readout (Hermite expansion,
Cov = sum_n b_n^2 n! rho^n with rho = exp(-k t), so int rho^n = 1/(n k)):

    D_cor = eps^2 / k * sum_{n>=1} b_n^2 n! / n,     b_n = E[fbar_Phi(R*+sigma_R Z) He_n(Z)] / n!

To leading order b_1 = sigma_R fbar_Phi'(R*) and

    D_Phi = (gamma T / 2) [ 1/R*^2 + (eps fbar_Phi'(R*) / k)^2 ]

    Q_inf ~= (gamma^2 / 2) [ 1 + (R* eps fbar_Phi'(R*) / k)^2 ]       (*)

THE CONTENT OF (*).  The bracket is the square of (sensitivity of rate to
amplitude) x (lifetime of an amplitude fluctuation).  It is ONE exactly when
d<Phi'>/dR = 0 at the limit cycle -- when the escapement is ISOCHRONOUS at its
own running amplitude.  Airy's 1826 condition, arriving as the condition for a
thermodynamic uncertainty law to be saturated.

CLOSED FORM FOR THE GRASSHOPPER.  f = -sgn(theta - theta_r sgn(theta')) gives
(derivation in the docstrings below)

    fbar_R(R)   =  2 theta_r / (pi R)                    (R > theta_r)
    fbar_Phi(R) = -2 sqrt(R^2 - theta_r^2) / (pi R^2)    (R > theta_r)

both identically the sub-critical constants 2/pi and 0 for R < theta_r.  Then
R*^2 = 4 eps theta_r / (pi gamma), k = gamma EXACTLY (because fbar_R goes as 1/R),
and with u = (R*/theta_r)^2 = 4 eps / (pi gamma theta_r) everything collapses:

    h(u) = 1 + (2 - u)^2 / (4 (u - 1))

  * h = 1 at u = 2, i.e. R* = sqrt(2) theta_r -- the isochronous amplitude.
  * h -> u/4 for large u, so Q_inf -> gamma^2 R*^2 / (8 theta_r^2): a LARGER swing
    makes the uncertainty product WORSE once past isochronism, and D_Phi saturates
    at gamma T / (8 theta_r^2) instead of falling as 1/R*^2.
  * h = 1 identically for R* < theta_r (fbar_Phi = 0 there: the escapement never
    reverses, so it is a pure sgn(theta') drive and carries no rate-amplitude slope).
  * van der Pol has fbar_Phi = 0 for ALL R by parity -- isochronous at every
    amplitude -- which is why finding #6's control came back flat.

Iris (Opus 5), Sep 9 2026.
"""
import math
import numpy as np

_fact = math.factorial
from numpy.polynomial.hermite_e import hermeval

# ---------------------------------------------------------------- closed forms
def gh_fbar_R(R, theta_r):
    """<sin(Phi) f> for the grasshopper, both branches.

    Upper half plane (sgn(theta')=+1): f = -sgn(R cos - theta_r).  With
    c = theta_r/R, Phi_c = arccos(c):  int_0^pi sin(Phi)(-sgn(cos-c)) dPhi = 2c.
    Lower half contributes another 2c by the theta -> -theta symmetry.
    Total / 2pi = 2 theta_r/(pi R).  Below theta_r the sign never flips and
    f = sgn(theta'), whose average is <|sin|> = 2/pi.
    """
    R = np.asarray(R, float)
    return np.where(R > theta_r, 2.0 * theta_r / (np.pi * np.maximum(R, 1e-300)),
                    2.0 / np.pi)

def gh_fbar_Phi(R, theta_r):
    """<cos(Phi) f>/R for the grasshopper.

    int_0^pi cos(Phi) (-sgn(cos - c)) dPhi = -2 sin(Phi_c) = -2 sqrt(1-c^2);
    the lower half gives the same; divide by 2 pi R:
        -2 sqrt(1 - theta_r^2/R^2) / (pi R) = -2 sqrt(R^2-theta_r^2)/(pi R^2).
    Below theta_r, f = sgn(theta') and <cos(Phi) sgn(sin(Phi))> = 0 exactly.
    """
    R = np.asarray(R, float)
    out = np.zeros_like(R)
    m = R > theta_r
    Rm = R[m] if R.ndim else R
    with np.errstate(invalid='ignore'):
        val = -2.0 * np.sqrt(np.maximum(R * R - theta_r * theta_r, 0.0)) / (np.pi * R * R)
    return np.where(m, val, 0.0)

def gh_h(u):
    """Q_inf / (gamma^2/2) for the grasshopper, closed form in u = (R*/theta_r)^2."""
    u = np.asarray(u, float)
    return np.where(u <= 1.0, 1.0, 1.0 + (2.0 - u) ** 2 / (4.0 * np.maximum(u - 1.0, 1e-300)))

def gh_u(gamma, eps, theta_r):
    """u = (R*/theta_r)^2 = 4 eps/(pi gamma theta_r) on the super-critical branch."""
    return 4.0 * eps / (np.pi * gamma * theta_r)


# ------------------------------------------------------- general numeric theory
def phase_averages(esc, R, P=40001):
    """(fbar_R, fbar_Phi) by quadrature on the ACTUAL (smoothed) escapement force."""
    Phi = np.linspace(0.0, 2.0 * np.pi, P)
    th, om = R * np.cos(Phi), R * np.sin(Phi)
    f = esc.f(th, om)
    return (np.trapz(np.sin(Phi) * f, Phi) / (2 * np.pi),
            np.trapz(np.cos(Phi) * f, Phi) / (2 * np.pi * R))

def _d(fn, R, h):
    return (fn(R + h) - fn(R - h)) / (2.0 * h)

def predict(esc, gamma, eps, T, Rstar, n_herm=24, hermite=True):
    """Predicted D_Phi and Q_inf.  Returns a dict; no fitted parameters anywhere.

    hermite=True uses the exact Gaussian-OU readout sum (needed near a kink, where
    sigma_R is not small against the distance to it); False is the b_1-only form.
    """
    dR = max(1e-4, 1e-3 * Rstar)
    fR = lambda R: phase_averages(esc, R)[0]
    fP = lambda R: phase_averages(esc, R)[1]
    k = gamma / 2.0 - eps * _d(fR, Rstar, dR)
    sig2 = gamma * T / (2.0 * k)
    sig = np.sqrt(sig2)
    D_dir = gamma * T / (2.0 * Rstar ** 2)

    if hermite:
        # Gauss-Hermite nodes for the standard normal
        x, w = np.polynomial.hermite_e.hermegauss(4 * n_herm)
        w = w / w.sum()
        g = np.array([phase_averages(esc, max(Rstar + sig * xi, 1e-6))[1] for xi in x])
        tot = 0.0
        for n in range(1, n_herm + 1):
            c = np.zeros(n + 1); c[n] = 1.0
            He = hermeval(x, c)
            bn = float(np.sum(w * g * He)) / float(_fact(n))
            tot += bn * bn * float(_fact(n)) / n
        D_cor = eps ** 2 * tot / k
        fP1 = float(np.sum(w * g * x)) / sig          # implied b_1/sigma
    else:
        fP1 = _d(fP, Rstar, dR)
        D_cor = (eps * fP1) ** 2 * gamma * T / (2.0 * k * k)
        tot = None
    return dict(k=k, sigma_R=sig, D_dir=D_dir, D_cor=D_cor, D_phi=D_dir + D_cor,
                fbarPhi_prime=fP1, h=1.0 + D_cor / D_dir)

def fit_bias_fraction(tau, t_max, frac_lo=0.25):
    """How much a THROUGH-ORIGIN slope fit of Var(J) under-reads D_cor.

    Var(J) = 2 D_dir t + 2 D_cor [t - tau (1 - e^{-t/tau})] -> A t - B, B = 2 D_cor tau.
    A through-origin fit over t in [frac_lo*t_max, t_max] returns A - B * <t>/<t^2>,
    so the deficit is B * <t>/<t^2> and it falls ENTIRELY on the D_cor term.
    """
    a, b = frac_lo * t_max, t_max
    mt = (b ** 2 - a ** 2) / 2.0 / (b - a)
    mt2 = (b ** 3 - a ** 3) / 3.0 / (b - a)
    return tau * mt / mt2
