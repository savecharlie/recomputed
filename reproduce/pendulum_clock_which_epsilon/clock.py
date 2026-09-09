"""
Stochastic pendulum clock: underdamped Langevin, phase current, entropy production.

Model (generalised from Izumida, arXiv:2609.04957, Eqs. 2-3):

    theta_dot = omega
    omega_dot = -theta - gamma*omega + eps*f(theta, omega) + sqrt(2*gamma*T) xi

The paper sets gamma == eps identically (its h = -theta_dot + f is multiplied by a
single eps).  Here the LINEAR DAMPING gamma and the ESCAPEMENT STRENGTH eps are
separate knobs, which is the whole point of this file.  Noise obeys FDT with
gamma, so D = gamma*T (the paper's D = eps*T is the gamma == eps case).

Current: J_Phi(t) = unwrapped Phi(t) - Phi(0), with Phi = atan2(omega, theta).
Stratonovich obeys the ordinary chain rule, so the integral of grad(Phi) o x_dot
IS the unwrapped geometric phase.  Nothing subtler is needed.

Iris (Opus 5), Sep 8 2026.
"""
import numpy as np

# ---------------------------------------------------------------- escapements
class VanDerPol:
    """f = (a - theta^2) * omega.   a = 2 recovers the paper's Eq. 24."""
    name = "van der Pol"
    def __init__(self, a=2.0):
        self.a = a
        self.params = {"a": a}
    def f(self, th, om):
        return (self.a - th * th) * om
    def dfdom(self, th, om):
        return self.a - th * th
    def fbar_R(self, R):
        # (1/2pi) int sin(Phi) f(R cos, R sin) dPhi  =  a R / 2 - R^3 / 8
        return self.a * R / 2.0 - R**3 / 8.0

class Grasshopper:
    """f = -sgn(theta - theta_r sgn(omega)), sgn smoothed as tanh(n x)."""
    name = "grasshopper escapement"
    def __init__(self, theta_r=0.4, n=20.0):
        self.theta_r, self.n = theta_r, n
        self.params = {"theta_r": theta_r, "n": n}
    def f(self, th, om):
        return -np.tanh(self.n * (th - self.theta_r * np.tanh(self.n * om)))
    def dfdom(self, th, om):
        u = th - self.theta_r * np.tanh(self.n * om)
        sech2_out = 1.0 - np.tanh(self.n * u) ** 2
        sech2_in = 1.0 - np.tanh(self.n * om) ** 2
        return sech2_out * self.n * self.theta_r * self.n * sech2_in
    def fbar_R(self, R):
        """PIECEWISE.  Table 1 prints 2 theta_r/(pi R), which is only the R > theta_r
        branch.  If the swing never reaches the critical angle the escapement never
        reverses, f = sgn(omega) throughout, and the phase average is the constant
        2/pi.  The paper never needs the other branch (its R* = sqrt(4 theta_r/pi)
        = 0.714 > theta_r = 0.4 for every epsilon, because epsilon cancels there).
        The two-parameter sweep DOES reach it: R*^2 = 4 eps theta_r/(pi gamma) drops
        below theta_r^2 once eps < pi gamma theta_r/4.  Using the printed branch
        below that gives R* about 12% too large.
        """
        R = np.asarray(R, dtype=float)
        return np.where(R > self.theta_r,
                        2.0 * self.theta_r / (np.pi * np.maximum(R, 1e-12)),
                        2.0 / np.pi)

class Graham:
    """f = -sgn(sin(psi0) theta - cos(psi0) omega), smoothed."""
    name = "Graham escapement"
    def __init__(self, psi0=0.1, n=20.0):
        self.psi0, self.n = psi0, n
        self.params = {"psi0": psi0, "n": n}
    def f(self, th, om):
        return -np.tanh(self.n * (np.sin(self.psi0) * th - np.cos(self.psi0) * om))
    def dfdom(self, th, om):
        z = np.sin(self.psi0) * th - np.cos(self.psi0) * om
        return (1.0 - np.tanh(self.n * z) ** 2) * self.n * np.cos(self.psi0)
    def fbar_R(self, R):
        return 2.0 * np.cos(self.psi0) / np.pi


def limit_cycle_radius(esc, gamma, eps, Rmax=60.0):
    """Solve -gamma R/2 + eps fbar_R(R) = 0 for the largest stable root.

    Paper's Eq. 19 is the gamma == eps case.  Bisection on a dense bracket so we
    do not depend on a good initial guess.
    """
    R = np.linspace(1e-4, Rmax, 400001)
    g = -gamma * R / 2.0 + eps * esc.fbar_R(R)
    s = np.sign(g)
    idx = np.nonzero((s[:-1] > 0) & (s[1:] <= 0))[0]   # downward crossing = stable
    if len(idx) == 0:
        return np.nan
    i = idx[-1]
    lo, hi = R[i], R[i + 1]
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if -gamma * mid / 2.0 + eps * esc.fbar_R(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------- integrator
def run(esc, gamma, eps, T, n_traj=2000, t_total=1000.0, dt=1e-3,
        t_burn=None, seed=0, chunks=40):
    """Euler-Maruyama on the ensemble.  Returns the measured observables.

    Everything is measured AFTER a burn-in so the ensemble sits on the
    stationary ring before the phase variance clock starts.
    """
    rng = np.random.default_rng(seed)
    Rstar = limit_cycle_radius(esc, gamma, eps)
    if t_burn is None:
        t_burn = max(15.0 / max(gamma, 1e-6), 150.0)   # ~7 amplitude relaxation times

    # start every trajectory on the deterministic cycle at the same phase
    th = np.full(n_traj, Rstar, dtype=np.float64)
    om = np.zeros(n_traj, dtype=np.float64)
    sig = np.sqrt(2.0 * gamma * T * dt)

    def step(th, om, k):
        z = rng.standard_normal((k, n_traj))
        for j in range(k):
            acc = -th - gamma * om + eps * esc.f(th, om)
            om = om + acc * dt + sig * z[j]
            th = th + om * dt
        return th, om

    n_burn = int(round(t_burn / dt))
    B = 4000
    done = 0
    while done < n_burn:
        k = min(B, n_burn - done)
        th, om = step(th, om, k)
        done += k

    # --- measurement phase -------------------------------------------------
    phi = np.arctan2(om, th)
    J = np.zeros(n_traj)                 # accumulated (unwrapped) phase
    n_steps = int(round(t_total / dt))
    per = n_steps // chunks
    ts, varJ, meanJ = [], [], []
    acc_om2 = acc_R2 = acc_dfdom = 0.0
    nacc = 0
    t = 0.0
    for c in range(chunks):
        for _sub in range(per // B + (1 if per % B else 0)):
            k = min(B, per - (_sub * B))
            if k <= 0:
                break
            z = rng.standard_normal((k, n_traj))
            for j in range(k):
                acc = -th - gamma * om + eps * esc.f(th, om)
                om = om + acc * dt + sig * z[j]
                th = th + om * dt
                new = np.arctan2(om, th)
                d = new - phi
                d -= 2.0 * np.pi * np.round(d / (2.0 * np.pi))   # unwrap
                J += d
                phi = new
                if j % 20 == 0:                                   # observables
                    R2 = th * th + om * om
                    acc_om2 += (om * om).mean()
                    acc_R2 += R2.mean()
                    acc_dfdom += esc.dfdom(th, om).mean()
                    nacc += 1
            t += k * dt
        ts.append(t); varJ.append(J.var(ddof=1)); meanJ.append(J.mean())

    ts = np.array(ts); varJ = np.array(varJ); meanJ = np.array(meanJ)
    # D_Phi from a straight-line fit through the origin on the last 3/4
    m = ts > ts[-1] * 0.25
    D_phi = float(np.sum(ts[m] * varJ[m]) / np.sum(ts[m] ** 2) / 2.0)
    Omega = float(np.sum(ts[m] * meanJ[m]) / np.sum(ts[m] ** 2))
    # Same slope, fitted WITH an intercept.  A through-origin fit systematically
    # UNDER-reads D_Phi whenever the phase noise has a finite correlation time:
    # Var(J) -> A t - B with B = 2 D_cor tau, so forcing the line through zero
    # tilts the slope down by B <t>/<t^2>.  Added Sep 9 2026 (fire 261) when the
    # amplitude-to-phase theory predicted exactly that bias; the intercept
    # absorbs B and the slope is unbiased for ANY finite-tau frequency noise.
    A_ = np.vstack([ts[m], np.ones(int(m.sum()))]).T
    sl_, ic_ = np.linalg.lstsq(A_, varJ[m], rcond=None)[0]
    D_phi_int = float(sl_ / 2.0)
    varJ_intercept = float(ic_)

    om2 = acc_om2 / nacc
    R2 = acc_R2 / nacc
    dfdom = acc_dfdom / nacc
    Teff = R2 / 2.0
    sigma_heat = gamma * (om2 - T) / T          # Sekimoto heat / T
    sigma_pump = eps * dfdom
    sigma_irr = sigma_heat - sigma_pump
    Q = 2.0 * D_phi / Omega ** 2 * sigma_irr

    return dict(Rstar=Rstar, D_phi=D_phi, D_phi_int=D_phi_int,
                varJ_intercept=varJ_intercept, Omega=Omega, om2=om2, R2=R2,
                Teff=Teff, sigma_heat=sigma_heat, sigma_pump=sigma_pump,
                sigma_irr=sigma_irr, Q=Q, ts=ts, varJ=varJ, meanJ=meanJ,
                gamma=gamma, eps=eps, T=T, model=esc.name)
