"""Validate 2 A^2 <z_om^2> before believing it, and make it predict something new.

isostable.py returned sqrt(2 A^2 <z_om^2>) = 0.9812 at u=3.5 against a measured
0.9792 +- 0.0030, and 0.9722 at u=6.37 against 0.9717 (cheap sweep).  Before that is
a result it has to (a) reproduce a limit whose answer is known and (b) predict a
number nobody has measured.

FULL RATIO.  What channels.py divides by is sigma_th = sqrt(gamma T / (2 k_avg)) with
k_avg the AVERAGED rate.  The exact variance is 2 A^2 <z_om^2> * gamma T/(2 k_true).
So the thing to compare with the measurement is

    ratio = sqrt( 2 A^2 <z_om^2> * k_avg / k_true )

and isostable.py quietly assumed k_avg = k_true.  True to 5 decimals for u >= 2,
FALSE at u = 1.2 where the Floquet exponent is 0.840 gamma.

KNOWN LIMIT.  Hold u fixed (so R* fixed) and send gamma, eps -> 0 together.  The orbit
becomes the exact harmonic circle, v -> rhat, z_om -> sin(Phi), and the averaged theory
is exact.  The ratio MUST go to 1, and it should go like O(gamma).

Iris (Opus 5), fire 265.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius
import amplitude_phase as ap
from isostable import rk4_orbit, jac

T = 5e-4


def full(u, gamma=0.05, n=20.0, theta_r=0.4, dt=2e-4):
    eps = math.pi * gamma * theta_r * u / 4.0
    esc = Grasshopper(theta_r, n)
    Rstar = limit_cycle_radius(esc, gamma, eps)
    X = rk4_orbit(esc, gamma, eps, Rstar, dt)
    ns = len(X) - 1
    P = ns * dt
    Js = np.array([jac(esc, gamma, eps, X[i, 0], X[i, 1]) for i in range(ns + 1)])
    M = np.eye(2)
    for i in range(ns):
        Ja, Jb = Js[i], Js[i + 1]; Jm = 0.5 * (Ja + Jb)
        m1 = Ja @ M; m2 = Jm @ (M + dt/2*m1); m3 = Jm @ (M + dt/2*m2); m4 = Jb @ (M + dt*m3)
        M = M + dt/6*(m1 + 2*m2 + 2*m3 + m4)
    ev, Rm = np.linalg.eig(M)
    j = int(np.argmax(np.abs(ev - 1.0)))
    k_true = -math.log(abs(float(np.real(ev[j])))) / P
    v0 = np.real(Rm[:, j]); v0 /= np.linalg.norm(v0)
    z0 = np.real(np.linalg.inv(Rm)[j, :]); z0 = z0 / float(z0 @ v0)

    def prop(vec, sign):
        out = np.empty((ns + 1, 2)); out[0] = vec; s = vec.copy()
        for i in range(ns):
            Ja, Jb = Js[i], Js[i + 1]; Jm = 0.5*(Ja + Jb)
            g = (lambda A, y: A @ y + k_true*y) if sign > 0 else (lambda A, y: -A.T @ y - k_true*y)
            k1 = g(Ja, s); k2 = g(Jm, s + dt/2*k1); k3 = g(Jm, s + dt/2*k2); k4 = g(Jb, s + dt*k3)
            s = s + dt/6*(k1 + 2*k2 + 2*k3 + k4); out[i+1] = s
        return out

    V, Z = prop(v0, +1), prop(z0, -1)
    Rt = np.hypot(X[:, 0], X[:, 1])
    A = float(np.mean(np.sum(V * (X / Rt[:, None]), axis=1)))
    zom2 = float(np.mean(Z[:, 1] ** 2))
    k_avg = ap.predict(esc, gamma, eps, T, Rstar)['k']
    shape = 2.0 * A * A * zom2
    return dict(u=u, gamma=gamma, n=n, eps=eps, Rstar=Rstar, P=P, dt=dt,
                k_true_over_gamma=k_true/gamma, k_avg_over_gamma=k_avg/gamma,
                A=A, zom2=zom2, shape=shape,
                ratio=math.sqrt(shape * k_avg / k_true))


print("1. TIMESTEP CONVERGENCE (u=3.5)")
for dt in (1e-3, 5e-4, 2e-4, 1e-4):
    r = full(3.5, dt=dt)
    print("   dt=%.0e  P=%.6f  A=%.6f  <z_om^2>=%.6f  ratio=%.6f" % (dt, r['P'], r['A'], r['zom2'], r['ratio']))

print("\n2. KNOWN LIMIT: u fixed at 3.5, gamma and eps -> 0 together.  ratio must -> 1.")
prev = None
for g in (0.05, 0.025, 0.0125, 0.00625):
    r = full(3.5, gamma=g, dt=2e-4)
    d = 1 - r['ratio']
    print("   gamma=%.5f  k_true/gamma=%.6f  A=%.6f  <z_om^2>=%.6f  ratio=%.6f  (1-ratio=%.5f%s)"
          % (g, r['k_true_over_gamma'], r['A'], r['zom2'], r['ratio'], d,
             "" if prev is None else ", halving factor %.2f" % (prev/d)))
    prev = d

print("\n3. PREDICTION for the gamma sweep already measured (u=3.5 held fixed)")
print("   measured: gamma=0.05 -> 0.9792 +- 0.0032 (4 seeds);  0.10 -> 0.9681;  0.20 -> 0.9856 (1 seed, +-0.0065)")
for g in (0.05, 0.1, 0.2):
    r = full(3.5, gamma=g, dt=2e-4)
    print("   gamma=%.2f  k_true/gamma=%.6f  ratio predicted = %.6f" % (g, r['k_true_over_gamma'], r['ratio']))

print("\n4. NEW PREDICTION, nothing measured here: the escapement smoothing n.")
print("   The whole mechanism is the sharp switch distorting the Floquet vectors,")
print("   so a smoother escapement must have a ratio closer to 1.  u=3.5, gamma=0.05.")
out = {}
for n in (5.0, 10.0, 20.0, 40.0, 80.0):
    r = full(3.5, n=n, dt=2e-4)
    out[n] = r
    print("   n=%5.1f  R*=%.6f  k_true/gamma=%.6f  k_avg/gamma=%.6f  A=%.6f  <z_om^2>=%.6f  RATIO=%.6f"
          % (n, r['Rstar'], r['k_true_over_gamma'], r['k_avg_over_gamma'], r['A'], r['zom2'], r['ratio']))
json.dump({str(k): v for k, v in out.items()}, open('isostable_n_prediction.json', 'w'), indent=1)
