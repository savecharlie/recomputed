"""The amplitude relaxation rate k, EXACTLY, with no phase averaging anywhere.

The averaged theory says k = gamma exactly for the grasshopper (fbar_R goes as 1/R,
so the escapement contributes -gamma/2 on top of the linear -gamma/2).  That is the
number sigma_R^2 = gamma T/(2k) = T/2 rests on.

But k is not a derived quantity: it is MINUS THE NONTRIVIAL FLOQUET EXPONENT of the
noiseless limit cycle, and that can be computed exactly by integrating the monodromy
matrix.  A periodic orbit of a 2D flow has multipliers 1 and mu = det M, and

    det M = exp( int_0^P tr J dt ),   tr J = -gamma + eps df/domega,

so  k_true = gamma - eps <df/domega>_time.  For the grasshopper df/domega >= 0
EVERYWHERE (it is sech^2 * n * theta_r * n * sech^2), so k_true < gamma strictly and
the averaged k is an over-estimate.  This measures by how much, three ways that must
agree: the trace integral, det of the integrated monodromy, and a direct decay fit.

sigma_R^2 = D_R/(2k).  If k is short by x% and D_R is not, sigma_R is HIGH by x/2%.
Measured sigma_R is LOW by 2.1%.  So this alone cannot be the answer -- but it is a
term that is definitely missing, and its size decides whether it matters.

Iris (Opus 5), fire 265.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius

THETA_R, N, GAMMA, T = 0.4, 20.0, 0.05, 5e-4


def onto_cycle(esc, gamma, eps, Rstar, dt, n_cyc=60):
    th, om = Rstar, 0.0
    for _ in range(int(round(n_cyc * 2 * math.pi / dt))):
        om += (-th - gamma * om + eps * esc.f(th, om)) * dt
        th += om * dt
    # advance to the next upward zero crossing of theta so the section is defined
    prev = th
    while not (prev < 0.0 <= th and om > 0):
        prev = th
        om += (-th - gamma * om + eps * esc.f(th, om)) * dt
        th += om * dt
    return th, om


def monodromy(esc, gamma, eps, th0, om0, dt):
    """RK4 the flow and the variational equation for one full return to the section."""
    def F(s):
        th, om = s[0], s[1]
        return np.array([om, -th - gamma * om + eps * esc.f(th, om)])

    def J(s):
        th, om = s[0], s[1]
        h = 1e-7
        fth = (esc.f(th + h, om) - esc.f(th - h, om)) / (2 * h)
        fom = esc.dfdom(th, om)
        return np.array([[0.0, 1.0], [-1.0 + eps * fth, -gamma + eps * fom]]), fom

    s = np.array([th0, om0])
    M = np.eye(2)
    tr_int = 0.0
    fom_int = 0.0
    t = 0.0
    prev_th = s[0]
    n = 0
    while True:
        Js, fom = J(s)
        # RK4 on the state; matrix propagated with the same substep Jacobians
        k1 = F(s);            J1, f1 = J(s)
        k2 = F(s + dt/2*k1);  J2, f2 = J(s + dt/2*k1)
        k3 = F(s + dt/2*k2);  J3, f3 = J(s + dt/2*k2)
        k4 = F(s + dt*k3);    J4, f4 = J(s + dt*k3)
        s_new = s + dt/6*(k1 + 2*k2 + 2*k3 + k4)
        m1 = J1 @ M
        m2 = J2 @ (M + dt/2*m1)
        m3 = J3 @ (M + dt/2*m2)
        m4 = J4 @ (M + dt*m3)
        M = M + dt/6*(m1 + 2*m2 + 2*m3 + m4)
        fom_int += dt/6*(f1 + 2*f2 + 2*f3 + f4)
        tr_int += dt/6*((-gamma+f1*eps) + 2*(-gamma+f2*eps) + 2*(-gamma+f3*eps) + (-gamma+f4*eps))
        t += dt
        n += 1
        if n > 20 and prev_th < 0.0 <= s_new[0] and s_new[1] > 0:
            # linear interpolation of the crossing time; the correction to M is O(dt)
            s = s_new
            break
        prev_th = s_new[0]
        s = s_new
    return M, t, tr_int, fom_int


rows = []
for u in (1.20, 2.00, 3.50, 5.00, 6.37):
    eps = math.pi * GAMMA * THETA_R * u / 4.0
    esc = Grasshopper(THETA_R, N)
    Rstar = limit_cycle_radius(esc, GAMMA, eps)
    out = {}
    for dt in (2e-4, 1e-4):
        th0, om0 = onto_cycle(esc, GAMMA, eps, Rstar, dt)
        M, P, tr_int, fom_int = monodromy(esc, GAMMA, eps, th0, om0, dt)
        ev = np.linalg.eigvals(M)
        mu = ev[np.argmin(np.abs(ev - 1.0)) ^ 1] if abs(ev[0]-1) != abs(ev[1]-1) else ev[1]
        # pick the multiplier furthest from 1
        mu = ev[int(np.argmax(np.abs(ev - 1.0)))]
        k_det = -math.log(abs(np.linalg.det(M))) / P
        k_ev = -math.log(abs(mu)) / P
        k_tr = -tr_int / P
        out[dt] = (P, k_det, k_ev, k_tr, eps * fom_int / P,
                   float(np.abs(ev - 1.0).min()))
    P, k_det, k_ev, k_tr, epsfom, triv = out[1e-4]
    rows.append(dict(u=u, eps=eps, Rstar=Rstar, P=P, k_det=k_det, k_ev=k_ev,
                     k_tr=k_tr, eps_fom=epsfom, trivial_dev=triv,
                     k_det_coarse=out[2e-4][1]))
    print("u=%.2f  P=%.5f  eps<df/dom>=%.6f" % (u, P, epsfom))
    print("    k/gamma:  det %.5f   eig %.5f   trace %.5f   (coarse dt: %.5f)"
          % (k_det/GAMMA, k_ev/GAMMA, k_tr/GAMMA, out[2e-4][1]/GAMMA))
    print("    trivial multiplier deviates from 1 by %.2e" % triv)
    print("    -> sigma_R^2 = D/(2k) would be %.5f x theory if D unchanged"
          % (GAMMA/k_det))
json.dump(rows, open('floquet.json', 'w'), indent=1)
