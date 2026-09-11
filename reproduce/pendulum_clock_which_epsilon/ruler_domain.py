"""Where does the no-filter estimator actually work?  Synthetic only, sigma typed in.

fit_robust.py refit the SAME gamma=0.20 data four ways and got 0.9317 / 0.9565 / 1.0563 /
1.2056 -- a 30% spread across fit bases -- while the same four bases agree to four decimals
at gamma=0.05.  Fire 265's own written-down criterion was: "if the mixed terms move gamma=0.2
and not gamma=0.05, the ruler was the problem."  They did.  So the ruler is the problem, and
the question is no longer about the clock at all, it is: what is this estimator's DOMAIN?

THE THING I MISSED, and it is the reason this file exists.  sigma_nofilter.validate() ran at
tau = 20 with a cycle period 6.0468 -- 3.3 cycles per correlation time, which is the gamma=0.05
regime.  It was then pointed at gamma = 0.20, where tau = 5 and there are 0.83 cycles per
correlation time, and it had never been checked there.  I validated the instrument in one
regime and deployed it in another.  That is the failure this project keeps catching, arriving
in the one place I thought was clean.

The controlling dimensionless number is N_c = tau / P, cycles per correlation time.  The fit
separates exp(-L/tau) from a Fourier series in Omega*L over L in [L0, 4 tau].  Large N_c: the
exponential decays across many cycles and is easy to tell from them.  N_c ~ 1: the exponential
is gone before one oscillation completes, and the two are nearly collinear.

PREDICTION, written before the first run:
 (a) every basis is unbiased for N_c >~ 2 and the spread between bases is small there;
 (b) the plain basis biases HIGH as N_c falls below ~1, because the exponential's short tail
     becomes hard to tell from the constant + the slow part of the Fourier series, so wobble
     power is absorbed into the coefficient I read out;
 (c) h6 is NOT systematically better on synthetic data.  There is no reason more harmonics
     should fix a collinearity between an exponential and a constant.  If h6 also fails here,
     then its landing on the isostable prediction at gamma=0.2 (0.9317 vs 0.9315) is a
     COINCIDENCE and must never be cited as confirmation.

(c) is the one that matters.  Picking the basis that agrees with my theory is exactly how a
person fools themselves, and the only defence is to choose the basis on data whose answer I
already know.

Iris (Opus 5), fire 266.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

P = 6.0468                     # cycle period of the real orbit at u = 3.5
OM = 2 * math.pi / P
DTR = 0.1                      # the real recording stride
NENS = 1024                    # the real ensemble size
BLOCK = 128                    # ensemble chunk, to keep peak memory ~150 MB


def autocov_chunked(make_block, n_blocks, NT, Lmax):
    """Unbiased lag autocovariance averaged over the ensemble, built a block at a time."""
    nf = 1 << (2 * NT - 1).bit_length()
    acc = np.zeros(Lmax + 1)
    ntraj = 0
    for b in range(n_blocks):
        x = make_block(b)
        d = x - x.mean(axis=0, keepdims=True)
        F = np.fft.rfft(d, n=nf, axis=0)
        ac = np.fft.irfft(F * np.conj(F), n=nf, axis=0)[:Lmax + 1]
        acc += ac.sum(axis=1)
        ntraj += x.shape[1]
        del x, d, F, ac
    counts = (NT - np.arange(Lmax + 1)).astype(float)
    return acc / ntraj / counts


def synth_block(tau, sigma, wobble, NT, seed, b, ntraj=BLOCK):
    """OU(sigma, tau) sampled at DTR + a periodic wobble at period P, random phase per traj."""
    rng = np.random.default_rng((seed, b))
    a = math.exp(-DTR / tau)
    s = sigma * math.sqrt(1 - a * a)
    x = np.empty((NT, ntraj))
    x[0] = sigma * rng.standard_normal(ntraj)
    for i in range(1, NT):
        x[i] = a * x[i - 1] + s * rng.standard_normal(ntraj)
    if wobble > 0:
        t = (np.arange(NT) * DTR)[:, None]
        ph = rng.uniform(0, 2 * math.pi, ntraj)[None, :]
        x += wobble * sigma * (np.cos(OM * t + ph) + 0.4 * np.cos(2 * OM * t + 2 * ph))
    return x


def fit(C, tau_guess, L0_cycles=1.0, nharm=4, mixed=False):
    """Best sigma for C(L) = s2 e^{-kL} + Fourier(Omega) + const, over L >= L0.  As fit_robust."""
    L = np.arange(len(C)) * DTR
    m = L >= L0_cycles * P
    Lm, Cm = L[m], C[m]

    def score(k):
        cols = [np.exp(-k * Lm), np.ones_like(Lm)]
        for jj in range(1, nharm + 1):
            cols += [np.cos(jj * OM * Lm), np.sin(jj * OM * Lm)]
            if mixed:
                cols += [np.exp(-k * Lm) * np.cos(jj * OM * Lm),
                         np.exp(-k * Lm) * np.sin(jj * OM * Lm)]
        A = np.column_stack(cols)
        c, *_ = np.linalg.lstsq(A, Cm, rcond=None)
        return float(np.sum((A @ c - Cm) ** 2)), float(c[0])

    lo, hi = 0.3 / tau_guess, 3.0 / tau_guess
    g = np.linspace(lo, hi, 81)
    best = min(((score(k), k) for k in g), key=lambda t: t[0][0])
    d = (hi - lo) / 80
    g2 = np.linspace(best[1] - d, best[1] + d, 41)
    best = min([best] + [(score(k), k) for k in g2], key=lambda t: t[0][0])
    return math.sqrt(max(best[0][1], 0.0)), best[1]


def measure(tau, wobble, seed, sigma=0.0158113939, t_total=3200.0):
    NT = int(round(t_total / DTR))
    Lmax = int(round(4 * tau / DTR))
    C = autocov_chunked(lambda b: synth_block(tau, sigma, wobble, NT, seed, b),
                        NENS // BLOCK, NT, Lmax)
    out = {}
    for name, kw in (("plain", {}), ("mixed", dict(mixed=True)),
                     ("L0", dict(L0_cycles=2.0)), ("h6", dict(nharm=6))):
        s, k = fit(C, tau, **kw)
        out[name] = s / sigma
        out["tau_" + name] = (1.0 / k) / tau
    return out


if __name__ == "__main__":
    NAMES = ("plain", "mixed", "L0", "h6")
    rows = []
    print("SCAN A -- cycles per correlation time, wobble fixed at 1.5 sigma")
    print("  the real runs sit at N_c = 3.31 (gamma=0.05), 1.65 (0.10), 0.83 (0.20)\n")
    print("%-7s %-6s %-5s | %-8s %-8s %-8s %-8s" % ("N_c", "tau", "seed", *NAMES))
    for Nc in (4.0, 3.31, 2.0, 1.65, 1.2, 0.83, 0.6):
        tau = Nc * P
        acc = {n: [] for n in NAMES}
        for seed in (101, 102, 103):
            v = measure(tau, 1.5, seed)
            for n in NAMES:
                acc[n].append(v[n])
            print("%-7.2f %-6.1f %-5d | %-8.4f %-8.4f %-8.4f %-8.4f"
                  % (Nc, tau, seed, *[v[n] for n in NAMES]), flush=True)
        print("%-7.2f %-6.1f %-5s | %-8.4f %-8.4f %-8.4f %-8.4f   <- mean\n"
              % (Nc, tau, "mean", *[float(np.mean(acc[n])) for n in NAMES]))
        rows.append(dict(scan="A", N_c=Nc, tau=tau, wobble=1.5,
                         **{n: acc[n] for n in NAMES}))
        json.dump(rows, open('ruler_domain.json', 'w'), indent=1)

    print("\nSCAN B -- wobble amplitude, at the two N_c the real runs actually sit at")
    print("%-7s %-8s %-5s | %-8s %-8s %-8s %-8s" % ("N_c", "wobble", "seed", *NAMES))
    for Nc in (3.31, 0.83):
        for wob in (0.0, 0.5, 1.5, 3.0, 6.0):
            acc = {n: [] for n in NAMES}
            for seed in (201, 202):
                v = measure(Nc * P, wob, seed)
                for n in NAMES:
                    acc[n].append(v[n])
            print("%-7.2f %-8.1f %-5s | %-8.4f %-8.4f %-8.4f %-8.4f"
                  % (Nc, wob, "mean", *[float(np.mean(acc[n])) for n in NAMES]), flush=True)
            rows.append(dict(scan="B", N_c=Nc, tau=Nc * P, wobble=wob,
                             **{n: acc[n] for n in NAMES}))
            json.dump(rows, open('ruler_domain.json', 'w'), indent=1)
        print()
