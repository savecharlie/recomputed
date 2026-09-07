"""Before saying anything about the slope diagnostic, settle the ruler.

"m_i(t) = visits to the rank-i site up to step t" admits two readings:
  (A) re-rank at every t   -- m_1(t) is a running maximum over a growing set,
      which is biased UP at early t and so biases the fitted slope DOWN;
  (B) rank once at the end and follow those sites' histories -- the natural
      reading of "cumulative visit counts to the rank-i most visited sites",
      and biased the other way.
The paper does not say which. The alpha inference from p1 is unaffected (the two
agree at t = T, which is where p1 is read), but the SLOPE is not, and the slope
is the paper's stated alpha diagnostic. So measure both.
"""
import numpy as np
from alpha_scan import fit_q, SPECIES, DAYS

R = 60


def both(T, alpha, beta, q, rng, cps, ranks=3):
    m = np.zeros(T + 2); m[0] = 1.0; n = 1
    A = np.zeros((len(cps), ranks))              # re-ranked
    hist = np.zeros((len(cps), T + 2))           # full snapshot for fixed ranking
    ci = 0
    for t in range(1, T + 1):
        if rng.random() < q * n ** (-beta):
            m[n] = 1.0; n += 1
        else:
            w = m[:n] ** alpha
            m[rng.choice(n, p=w / w.sum())] += 1.0
        if ci < len(cps) and t == cps[ci]:
            k = min(ranks, n)
            top = np.sort(np.partition(m[:n], n - k)[n - k:])[::-1]
            A[ci, :k] = top
            hist[ci, :n] = m[:n]
            ci += 1
    order = np.argsort(m[:n])[::-1][:ranks]      # ranks fixed at the END
    B = hist[:, order]
    return A, B


if __name__ == "__main__":
    for sp, d in SPECIES.items():
        T = int(round(DAYS * d["lam"])); beta = d["beta"]
        q = fit_q(beta, T, d["n200"])
        lo = int(round(10 * d["lam"]))
        cps = sorted(set(np.unique(np.round(np.logspace(np.log10(lo), np.log10(T), 24))).astype(int)))
        print(f"\n{sp}  T={T}  beta={beta}  q={q:.3f}   (slope of a true power law t^1 is 1.000)")
        print("   alpha |  re-ranked each t        |  ranks fixed at t=T")
        print("         |  i=1     i=2     i=3     |  i=1     i=2     i=3")
        for alpha in (0.8, 0.9, 1.0, 1.1, 1.2):
            rng = np.random.default_rng(7)
            SA = np.zeros((len(cps), 3)); SB = np.zeros((len(cps), 3))
            for _ in range(R):
                a, b = both(T, alpha, beta, q, rng, cps)
                SA += a; SB += b
            SA /= R; SB /= R
            f = lambda Y: [np.polyfit(np.log(cps), np.log(np.maximum(Y[:, i], 1e-9)), 1)[0] for i in range(3)]
            sa, sb = f(SA), f(SB)
            print(f"    {alpha:4.2f} |  " + " ".join(f"{x:.3f}" for x in sa)
                  + "  |  " + " ".join(f"{x:.3f}" for x in sb), flush=True)
