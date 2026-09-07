"""alpha is NOT fitted in the paper. It is read off an eye-guide. Here is what
the same data can actually say about it.

The structural point that makes this well-posed, and which the paper does not
state: in the IMM the exploration probability Pnew = q n^-beta depends on n
ALONE, so the discovery curve n(t) carries no information about alpha, and alpha
carries none about n(t). (beta, q) are therefore fully determined by Fig. 4a, and
alpha is then determined by how the returns are distributed -- which is exactly
what Fig. 5 measures. The paper fits the first half and eyeballs the second.

Inputs, all printed in the paper except n(200 d), which is read off Fig. 4a:
    FT: beta 0.72, 3.74 visits/day, p1 = 0.28, n(200 d) ~ 57
    BN: beta 0.37, 5.49 visits/day, p1 = 0.11, n(200 d) ~ 190
q is set per species so that the simulated n at t = 200*lambda steps matches.
"""
import numpy as np
from multiprocessing import Pool

SPECIES = {
    "FT": dict(beta=0.72, lam=3.74, p1=0.28, n200=57.0),
    "BN": dict(beta=0.37, lam=5.49, p1=0.11, n200=190.0),
}
DAYS = 200


def run(T, alpha, beta, q, rng, ranks=10, checks=()):
    m = np.zeros(T + 2)
    m[0] = 1.0
    n = 1
    out = {}
    cs = set(checks)
    for t in range(1, T + 1):
        if rng.random() < q * n ** (-beta):
            m[n] = 1.0
            n += 1
        else:
            w = m[:n] ** alpha
            j = rng.choice(n, p=w / w.sum())
            m[j] += 1.0
        if t in cs:
            k = min(ranks, n)
            top = np.sort(np.partition(m[:n], n - k)[n - k:])[::-1]
            out[t] = np.pad(top, (0, ranks - k)) / t
    k = min(ranks, n)
    top = np.sort(np.partition(m[:n], n - k)[n - k:])[::-1]
    return n, np.pad(top, (0, ranks - k)) / T, out


def fit_q(beta, T, n_target, seed=3, R=40):
    """q from the discovery curve alone (alpha is irrelevant here -- verified)."""
    lo, hi = 0.05, 5.0
    for _ in range(22):
        q = 0.5 * (lo + hi)
        rng = np.random.default_rng(seed)
        ns = [run(T, 1.0, beta, q, rng, ranks=1)[0] for _ in range(R)]
        (lo, hi) = (q, hi) if np.mean(ns) < n_target else (lo, q)
    return 0.5 * (lo + hi)


def job(a):
    sp, alpha, q, T, R, seed = a
    b = SPECIES[sp]["beta"]
    rng = np.random.default_rng(seed)
    P = np.array([run(T, alpha, b, q, rng)[1] for _ in range(R)])
    return sp, alpha, P.mean(0), P[:, 0].std()


if __name__ == "__main__":
    ALPHAS = [round(a, 2) for a in np.arange(0.50, 1.51, 0.05)]
    tasks, meta = [], {}
    for sp, d in SPECIES.items():
        T = int(round(DAYS * d["lam"]))
        q = fit_q(d["beta"], T, d["n200"])
        meta[sp] = (T, q)
        print(f"{sp}: T = {T} steps, beta = {d['beta']}, fitted q = {q:.3f}", flush=True)
        for a in ALPHAS:
            tasks.append((sp, a, q, T, 40, 900 + int(a * 100)))
    with Pool(5) as pool:
        res = pool.map(job, tasks)

    print("\nmean top-site share p1 at day 200, vs alpha  (paper: FT 0.28, BN 0.11)")
    print("  alpha " + "".join(f"{a:7.2f}" for a in ALPHAS))
    store = {}
    for sp in SPECIES:
        row = {a: (p, s) for s2, a, p, s in res if s2 == sp}
        store[sp] = row
        print(f"  {sp}    " + "".join(f"{row[a][0][0]:7.3f}" for a in ALPHAS))
        print(f"   sd    " + "".join(f"{row[a][1]:7.3f}" for a in ALPHAS))
    for sp, d in SPECIES.items():
        xs = np.array(ALPHAS); ys = np.array([store[sp][a][0][0] for a in xs])
        i = int(np.argmin(np.abs(ys - d["p1"])))
        print(f"\n{sp}: observed p1 = {d['p1']} -> closest simulated alpha = {xs[i]} "
              f"(p1_sim {ys[i]:.3f})")
        sd = store[sp][xs[i]][1]
        ok = xs[(ys > d['p1'] - sd) & (ys < d['p1'] + sd)]
        if len(ok):
            print(f"   alphas within 1 individual-sd ({sd:.3f}) of the observed p1: "
                  f"{ok.min():.2f} .. {ok.max():.2f}")
    np.save("alpha_scan.npy", np.array([[store[sp][a][0] for a in ALPHAS] for sp in SPECIES]))
