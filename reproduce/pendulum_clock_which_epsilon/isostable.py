"""Where the 2.1% went: the slow variable is not R, it is the isostable coordinate.

k is exact (floquet.py: Floquet exponent = gamma to 5 decimals for u >= 2).  So the
whole 4% shortfall in sigma_R^2 = gamma T/(2k) has to be in the DIFFUSION, and the
diffusion is where averaging the equation and solving the equation part company.

EXACT STATEMENT.  Near a stable limit cycle x*(t) of period P the only slow variable
is the isostable coordinate psi = z(t).dx, where z is the P-periodic ADJOINT Floquet
eigenvector,   zdot = -J^T z - k z,   and perturbations relax along the P-periodic
right eigenvector v,   vdot = J v + k v,   normalised z.v = 1.  With noise sqrt(2 gamma T)
on omega only,

    dpsi = -k psi dt + z_omega(t) sqrt(2 gamma T) dW      ->   Var(psi) = gamma T <z_om^2>_t / k

and the deviation of the CYCLE-AVERAGED radius that this produces is

    dRbar = psi * A,      A = < v(t) . xhat*(t) >_t

so the exact prediction for what channels.py measures is

    Var(Rbar) / (gamma T / 2k)  =  2 A^2 <z_omega^2>_t .

For a weakly perturbed harmonic oscillator v is radial and z_omega = sin(Phi), giving
2 * 1 * (1/2) = 1 -- the averaged theory.  Every departure of the true Floquet vectors
from radial/sinusoidal is a correction the averaging threw away.

MEASURED TARGET: sigma_R^meas / sqrt(T/2) = 0.9793 +- 0.0030, i.e. this ratio = 0.959.

Iris (Opus 5), fire 265.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius

THETA_R, N, GAMMA, T = 0.4, 20.0, 0.05, 5e-4


def f_theta(esc, th, om, h=1e-7):
    return (esc.f(th + h, om) - esc.f(th - h, om)) / (2 * h)


def jac(esc, gamma, eps, th, om):
    return np.array([[0.0, 1.0],
                     [-1.0 + eps * f_theta(esc, th, om), -gamma + eps * esc.dfdom(th, om)]])


def rk4_orbit(esc, gamma, eps, Rstar, dt, n_cyc=80):
    """Land on the cycle, then return one full period sampled at dt (RK4)."""
    def F(s):
        return np.array([s[1], -s[0] - gamma * s[1] + eps * esc.f(s[0], s[1])])
    def step(s):
        k1 = F(s); k2 = F(s + dt/2*k1); k3 = F(s + dt/2*k2); k4 = F(s + dt*k3)
        return s + dt/6*(k1 + 2*k2 + 2*k3 + k4)
    s = np.array([Rstar, 0.0])
    for _ in range(int(round(n_cyc * 2*math.pi/dt))):
        s = step(s)
    while not (s[0] < 0 <= step(s)[0] and s[1] > 0):
        s = step(s)
    s = step(s)
    traj = [s.copy()]
    prev = s[0]
    while True:
        s = step(s)
        traj.append(s.copy())
        if len(traj) > 20 and prev < 0 <= s[0] and s[1] > 0:
            break
        prev = s[0]
    return np.array(traj)


def analyse(u, dt=2e-4, verbose=True, lam=1.0):
    """lam rescales the right Floquet vector v -> lam*v (and so z -> z/lam).

    It exists to make a point checkable rather than merely argued: |v(0)| = 1 is a
    CONVENTION, so A and <z_om^2> are separately meaningless and only their product is
    physical.  See gauge.py.  Default 1.0 leaves every existing result untouched.
    """
    eps = math.pi * GAMMA * THETA_R * u / 4.0
    esc = Grasshopper(THETA_R, N)
    Rstar = limit_cycle_radius(esc, GAMMA, eps)
    X = rk4_orbit(esc, GAMMA, eps, Rstar, dt)
    n = len(X) - 1
    P = n * dt
    Js = np.array([jac(esc, GAMMA, eps, X[i, 0], X[i, 1]) for i in range(n + 1)])

    # --- monodromy, and k
    M = np.eye(2)
    for i in range(n):
        Ja, Jb = Js[i], Js[i + 1]
        Jm = 0.5 * (Ja + Jb)
        m1 = Ja @ M
        m2 = Jm @ (M + dt/2*m1)
        m3 = Jm @ (M + dt/2*m2)
        m4 = Jb @ (M + dt*m3)
        M = M + dt/6*(m1 + 2*m2 + 2*m3 + m4)
    ev, R_ = np.linalg.eig(M)
    j = int(np.argmax(np.abs(ev - 1.0)))          # the non-trivial multiplier
    mu = float(np.real(ev[j]))
    k = -math.log(abs(mu)) / P
    v0 = np.real(R_[:, j]); v0 /= np.linalg.norm(v0); v0 *= lam
    L = np.linalg.inv(R_)
    z0 = np.real(L[j, :])
    z0 = z0 / float(z0 @ v0)                       # z.v = 1

    # --- propagate the two P-periodic Floquet vectors
    def prop(vec, sign):
        """sign=+1: vdot = J v + k v.   sign=-1: zdot = -J^T z - k z."""
        out = np.empty((n + 1, 2)); out[0] = vec
        s = vec.copy()
        for i in range(n):
            Ja, Jb = Js[i], Js[i + 1]; Jm = 0.5*(Ja + Jb)
            if sign > 0:
                g = lambda A, y: A @ y + k * y
            else:
                g = lambda A, y: -A.T @ y - k * y
            k1 = g(Ja, s); k2 = g(Jm, s + dt/2*k1); k3 = g(Jm, s + dt/2*k2); k4 = g(Jb, s + dt*k3)
            s = s + dt/6*(k1 + 2*k2 + 2*k3 + k4)
            out[i + 1] = s
        return out

    V = prop(v0, +1)
    Z = prop(z0, -1)
    per_v = float(np.linalg.norm(V[-1] - V[0]) / np.linalg.norm(V[0]))
    per_z = float(np.linalg.norm(Z[-1] - Z[0]) / np.linalg.norm(Z[0]))
    norm = float(np.mean([Z[i] @ V[i] for i in range(n + 1)]))   # should be 1 everywhere

    Rt = np.hypot(X[:, 0], X[:, 1])
    xhat = X / Rt[:, None]
    A = float(np.mean(np.sum(V * xhat, axis=1)))
    zom2 = float(np.mean(Z[:, 1] ** 2))
    ratio = 2.0 * A * A * zom2

    # the averaging approximation, for contrast: v radial, z_om = sin Phi
    Phi = np.arctan2(X[:, 1], X[:, 0])
    naive = 2.0 * 1.0 * float(np.mean(np.sin(Phi) ** 2))

    if verbose:
        print("u=%.2f  eps=%.5f  R*=%.5f  P=%.5f  k/gamma=%.6f" % (u, eps, Rstar, P, k/GAMMA))
        print("   periodicity check  |V(P)-V(0)|/|V| = %.2e   |Z(P)-Z(0)|/|Z| = %.2e"
              % (per_v, per_z))
        print("   z.v along the orbit (should be 1) = %.8f" % norm)
        print("   A = <v . rhat>_t        = %.6f" % A)
        print("   <z_omega^2>_t           = %.6f" % zom2)
        print("   2 A^2 <z_om^2>          = %.6f   -> sigma ratio %.6f"
              % (ratio, math.sqrt(max(ratio, 0))))
        print("   (naive 2<sin^2Phi>_t    = %.6f)" % naive)
    return dict(u=u, eps=eps, Rstar=Rstar, P=P, k_over_gamma=k/GAMMA, A=A,
                zom2=zom2, var_ratio=ratio, sigma_ratio=math.sqrt(max(ratio, 0)),
                per_v=per_v, per_z=per_z, zv=norm, naive=naive)


if __name__ == "__main__":
    rows = [analyse(u) for u in (1.20, 2.00, 3.50, 5.00, 6.37)]
    print("\nMEASURED sigma_R/sqrt(T/2): 0.9792 (u=3.5), 0.9774 (u=6.37)")
    json.dump(rows, open('isostable.json', 'w'), indent=1)
