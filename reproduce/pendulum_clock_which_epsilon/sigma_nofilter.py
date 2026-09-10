"""sigma_R and tau with NO filter and NO correction factor anywhere.

Every sigma_R in this project has come through a one-period boxcar and then been
divided by an analytic OU-through-a-boxcar factor.  At gamma = 0.05 that factor is
0.950 and the correction is small; at gamma = 0.2 the amplitude correlation time is
tau = 5 and the window is 6.3, so the window is LONGER THAN THE THING IT IS SMOOTHING
and the factor is 0.828.  isostable.py predicts 0.9315 there and the boxcar ruler
measured 0.9856 -- eight sigma apart -- while the two smaller-gamma points agree.
Before that is a failure of the theory it has to survive a ruler with no filter in it.

THE RULER.  R_raw = x + w, x the slow amplitude, w the within-cycle wobble, which is
periodic and (after removing the time mean) zero-mean.  The lag-autocovariance is then

    C(L) = sigma^2 e^{-kL}  +  [periodic in L]  +  cross terms that vanish
                                                   because the cycle phase is uniform

so fit C(L) with an exponential PLUS a Fourier series at the known cycle frequency and
read sigma^2 off the exponential's coefficient.  Linear in everything but k, so: grid
k, solve least squares, take the best.  No smoothing, no window alignment, no factor.

VALIDATED FIRST on a synthetic signal with a sigma I type in myself and a wobble
LARGER than the slow fluctuation, at the true cycle period.

Iris (Opus 5), fire 265.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius
import amplitude_phase as ap

T_DEF, TR = 5e-4, 0.4


def autocov(x, Lmax):
    """Unbiased lag autocovariance, averaged over the trajectory axis. FFT."""
    d = x - x.mean(axis=0, keepdims=True)
    n = d.shape[0]
    nf = 1 << (2 * n - 1).bit_length()
    F = np.fft.rfft(d, n=nf, axis=0)
    ac = np.fft.irfft(F * np.conj(F), n=nf, axis=0)[:Lmax + 1]
    counts = (n - np.arange(Lmax + 1)).astype(float)
    return (ac.mean(axis=1) / counts)


def fit(C, dtr, Omega, k_lo, k_hi, L0, nharm=4):
    """Best (sigma^2, k) for C(L) = s2 e^{-kL} + Fourier(Omega) + const, L >= L0."""
    L = np.arange(len(C)) * dtr
    m = L >= L0
    Lm, Cm = L[m], C[m]
    cols = [np.ones_like(Lm)]
    for j in range(1, nharm + 1):
        cols += [np.cos(j * Omega * Lm), np.sin(j * Omega * Lm)]
    best = None
    for k in np.linspace(k_lo, k_hi, 601):
        A = np.column_stack([np.exp(-k * Lm)] + cols)
        c, res, *_ = np.linalg.lstsq(A, Cm, rcond=None)
        r = float(np.sum((A @ c - Cm) ** 2))
        if best is None or r < best[0]:
            best = (r, float(c[0]), k)
    return best[1], best[2]


# ---------------------------------------------------------------- validation
def validate():
    print("VALIDATION -- known sigma, wobble LARGER than the slow signal")
    tau, dtr, sigma = 20.0, 0.1, 0.0158113939
    Om = 2 * math.pi / 6.0468
    NT, NENS = 24000, 512
    rng = np.random.default_rng(11)
    a = math.exp(-dtr / tau); s = sigma * math.sqrt(1 - a * a)
    slow = np.empty((NT, NENS)); slow[0] = sigma * rng.standard_normal(NENS)
    for i in range(1, NT):
        slow[i] = a * slow[i - 1] + s * rng.standard_normal(NENS)
    t = (np.arange(NT) * dtr)[:, None]
    ph = rng.uniform(0, 2 * math.pi, NENS)[None, :]
    for name, amp in (("no wobble", 0.0), ("wobble 1.5x sigma", 1.5 * sigma)):
        x = slow + (amp * (np.cos(Om * t + ph) + 0.4 * np.cos(2 * Om * t + 2 * ph)) if amp else 0.0)
        C = autocov(x, int(round(4 * tau / dtr)))
        s2, k = fit(C, dtr, Om, 0.02, 0.10, L0=2 * math.pi / Om)
        print("   %-20s  sigma/true = %.4f   tau = %.2f (true 20.0)"
              % (name, math.sqrt(max(s2, 0)) / sigma, 1 / k))


# ---------------------------------------------------------------- the system
def run(u, gamma, dt=2e-3, n_traj=1024, t_total=4800.0, t_burn=None, seed=5, T=T_DEF, n=20.0):
    esc = Grasshopper(TR, n)
    eps = math.pi * gamma * TR * u / 4.0
    Rstar = limit_cycle_radius(esc, gamma, eps)
    p = ap.predict(esc, gamma, eps, T, Rstar)
    if t_burn is None:
        t_burn = 30.0 / gamma
    rng = np.random.default_rng(seed)
    sig = math.sqrt(2.0 * gamma * T * dt)
    th = np.full(n_traj, Rstar); om = np.zeros(n_traj)
    for _ in range(int(round(t_burn / dt))):
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
    R = R[:j]
    tau = 1.0 / p['k']
    C = autocov(R, int(round(4 * tau / dtr)))
    # cycle frequency from the noiseless orbit, via the deterministic period
    Om = _det_omega(esc, gamma, eps, Rstar)
    s2, k = fit(C, dtr, Om, 0.3 / tau, 3.0 / tau, L0=2 * math.pi / Om)
    return dict(u=u, gamma=gamma, seed=seed, sigma_th=p['sigma_R'], tau_th=tau,
                sigma_meas=math.sqrt(max(s2, 0)), tau_meas=1.0 / k,
                ratio=math.sqrt(max(s2, 0)) / p['sigma_R'], tau_ratio=(1.0 / k) / tau)


def _det_omega(esc, gamma, eps, Rstar, dt=2e-4):
    th, om = Rstar, 0.0
    for _ in range(int(round(60 * 2 * math.pi / dt))):
        om += (-th - gamma * om + eps * esc.f(th, om)) * dt; th += om * dt
    prev = th; t = 0.0; cr = 0
    while True:
        om += (-th - gamma * om + eps * esc.f(th, om)) * dt; th += om * dt
        if prev < 0 <= th and om > 0:
            cr += 1
            if cr == 1:
                t = 0.0
            else:
                return 2 * math.pi / t
        prev = th; t += dt


if __name__ == "__main__":
    validate()
    print("\nTHE SYSTEM.  isostable.py predicts:  gamma=0.05 -> 0.9813,  gamma=0.20 -> 0.9315")
    print("boxcar ruler measured:               gamma=0.05 -> 0.9792,  gamma=0.20 -> 0.9856")
    rows = []
    for gamma, seeds in ((0.05, (5, 6, 7)), (0.20, (5, 6, 7))):
        rs = []
        for sd in seeds:
            r = run(3.5, gamma, seed=sd)
            rs.append(r['ratio']); rows.append(r)
            print("   gamma=%.2f seed=%d   sigma/th = %.4f   tau/th = %.4f"
                  % (gamma, sd, r['ratio'], r['tau_ratio']), flush=True)
        a = np.array(rs)
        print("   gamma=%.2f  MEAN %.4f  s.d. %.4f  s.e.m %.4f\n"
              % (gamma, a.mean(), a.std(ddof=1), a.std(ddof=1) / math.sqrt(len(a))))
    json.dump(rows, open('sigma_nofilter.json', 'w'), indent=1)
