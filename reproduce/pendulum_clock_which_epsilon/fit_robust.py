"""Is the no-filter ruler itself breaking at gamma = 0.2?  Refit the SAME data 4 ways.

sigma_nofilter says 0.9806 +- 0.0011 at gamma=0.05 (prediction 0.9813, agree) and
0.9569 +- 0.0010 at gamma=0.20 (prediction 0.9315, 25 sigma apart).  Before that is a
statement about physics, the new instrument gets the treatment every other ruler in this
project has had.

THE TERM I KNOW IS MISSING.  A perturbation off the cycle has an isostable component AND
a phase component, and both are driven by the same noise, so they are correlated.  The
cross term in C_R(L) goes like e^{-kL} cos(j Omega L) -- an exponential TIMES an
oscillation -- and the fit basis contains e^{-kL} and cos(j Omega L) separately but not
their product.  Whatever of it is not orthogonal to e^{-kL} is absorbed into the number I
am reading out.  It scales with the wobble, which scales with gamma, so it should be
negligible at 0.05 and not at 0.2.

PREDICTION, before running: adding the mixed terms moves gamma=0.20 and leaves gamma=0.05
where it is.  If instead BOTH move, the ruler was never clean.  If NEITHER moves, the
disagreement at 0.2 is real physics and the linear-response result has a domain.

Iris (Opus 5), fire 265.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius
import amplitude_phase as ap
from sigma_nofilter import autocov, _det_omega, TR, T_DEF


def sim(u, gamma, seed, dt=2e-3, n_traj=1024, t_total=3200.0):
    esc = Grasshopper(TR, 20.0)
    eps = math.pi * gamma * TR * u / 4.0
    Rstar = limit_cycle_radius(esc, gamma, eps)
    p = ap.predict(esc, gamma, eps, T_DEF, Rstar)
    rng = np.random.default_rng(seed)
    sig = math.sqrt(2.0 * gamma * T_DEF * dt)
    th = np.full(n_traj, Rstar); om = np.zeros(n_traj)
    for _ in range(int(round((30.0 / gamma) / dt))):
        om += (-th - gamma * om + eps * esc.f(th, om)) * dt + sig * rng.standard_normal(n_traj)
        th += om * dt
    stride = max(1, int(round(0.1 / dt))); dtr = stride * dt
    nst = int(round(t_total / dt))
    R = np.empty((nst // stride, n_traj)); j = 0
    for i in range(nst):
        om += (-th - gamma * om + eps * esc.f(th, om)) * dt + sig * rng.standard_normal(n_traj)
        th += om * dt
        if (i + 1) % stride == 0 and j < R.shape[0]:
            R[j] = np.hypot(th, om); j += 1
    tau = 1.0 / p['k']
    C = autocov(R[:j], int(round(4 * tau / dtr)))
    return C, dtr, _det_omega(esc, gamma, eps, Rstar), p['sigma_R'], tau


def fit(C, dtr, Om, tau, L0_cycles=1.0, nharm=4, mixed=False):
    L = np.arange(len(C)) * dtr
    m = L >= L0_cycles * 2 * math.pi / Om
    Lm, Cm = L[m], C[m]
    def score(k):
        cols = [np.exp(-k * Lm), np.ones_like(Lm)]
        for jj in range(1, nharm + 1):
            cols += [np.cos(jj * Om * Lm), np.sin(jj * Om * Lm)]
            if mixed:
                cols += [np.exp(-k * Lm) * np.cos(jj * Om * Lm),
                         np.exp(-k * Lm) * np.sin(jj * Om * Lm)]
        A = np.column_stack(cols)
        c, *_ = np.linalg.lstsq(A, Cm, rcond=None)
        return float(np.sum((A @ c - Cm) ** 2)), float(c[0])

    # two-stage grid: same resolution as a single 801-point sweep, 6x fewer solves
    lo, hi = 0.3 / tau, 3.0 / tau
    g = np.linspace(lo, hi, 81)
    best = min(((score(k), k) for k in g), key=lambda t: t[0][0])
    d = (hi - lo) / 80
    g2 = np.linspace(best[1] - d, best[1] + d, 41)
    best = min([best] + [((score(k)), k) for k in g2], key=lambda t: t[0][0])
    return math.sqrt(max(best[0][1], 0.0))


if __name__ == "__main__":
    rows = []
    print("%-7s %-5s | %-9s %-9s %-9s %-9s" % ("gamma", "seed", "plain", "+mixed", "L0=2cyc", "harm=6"))
    for gamma, pred in ((0.05, 0.9813), (0.20, 0.9315)):
        acc = {k: [] for k in ("plain", "mixed", "L0", "h6")}
        for seed in (5, 6, 7):
            C, dtr, Om, s_th, tau = sim(3.5, gamma, seed)
            v = dict(plain=fit(C, dtr, Om, tau) / s_th,
                     mixed=fit(C, dtr, Om, tau, mixed=True) / s_th,
                     L0=fit(C, dtr, Om, tau, L0_cycles=2.0) / s_th,
                     h6=fit(C, dtr, Om, tau, nharm=6) / s_th)
            for k in acc:
                acc[k].append(v[k])
            print("%-7.2f %-5d | %-9.4f %-9.4f %-9.4f %-9.4f"
                  % (gamma, seed, v['plain'], v['mixed'], v['L0'], v['h6']), flush=True)
        print("%-7.2f %-5s | %-9.4f %-9.4f %-9.4f %-9.4f   <- mean   (isostable predicts %.4f)"
              % (gamma, "mean", *[np.mean(acc[k]) for k in ("plain", "mixed", "L0", "h6")], pred))
        rows.append(dict(gamma=gamma, predicted=pred, **{k: acc[k] for k in acc}))
    json.dump(rows, open('fit_robust.json', 'w'), indent=1)
