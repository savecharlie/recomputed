"""How noisy is the D_cor estimator itself?  Twenty synthetic records, exact known answer.

cov_shape2 showed the same estimator returning 0.945 at one seed and 1.121 at another, so
the 5% 'deficit' in CHANNELS_FINDINGS may be scatter and not bias.  This measures the
estimator's own distribution on a Gaussian OU with the SAME tau, sigma, dtr, record length
and ensemble size as the pendulum runs, where D_cor = sigma^2 tau exactly.

The lag products are done by FFT, which is algebraically identical to
(gd[:-L]*gd[L:]).mean() summed over trajectories -- verified against the direct loop below.

Iris (Opus 5), fire 264.
"""
import numpy as np, math, json, sys

tau, dtr, sigma = 20.0, 0.1, 0.015811393924397392
NT, NENS = 12000, 512
LC = int(round(8.0*tau/dtr))


def autocov(gd, Lc):
    """cov[L] = sum_traj sum_t gd[t]gd[t+L] / ((NT-L)*NENS), same as the direct loop."""
    n = gd.shape[0]
    nf = 1 << (2*n - 1).bit_length()
    F = np.fft.rfft(gd, n=nf, axis=0)
    ac = np.fft.irfft(F*np.conj(F), n=nf, axis=0)[:Lc+1].sum(axis=1)
    return ac / (np.arange(n, n-Lc-1, -1) * gd.shape[1])


def one(seed):
    rng = np.random.default_rng(seed)
    a = math.exp(-dtr/tau); s = sigma*math.sqrt(1-a*a)
    R = np.empty((NT, NENS)); R[0] = sigma*rng.standard_normal(NENS)
    for i in range(1, NT):
        R[i] = a*R[i-1] + s*rng.standard_normal(NENS)
    gd = R - R.mean()
    cov = autocov(gd, LC)
    zc = np.nonzero(cov <= 0)[0]
    cut = int(zc[0]) if len(zc) else len(cov)
    D = float(np.trapz(cov[:cut], dx=dtr))
    return D/(sigma**2*tau), cut*dtr, float(cov[0])/sigma**2


if __name__ == '__main__':
    # instrument validation: FFT path must equal the direct loop it replaces
    rng = np.random.default_rng(0); x = rng.standard_normal((400, 5))
    xd = x - x.mean(); ref = np.array([float((xd[:-L]*xd[L:]).mean()) for L in range(1, 21)])
    got = autocov(xd, 20)[1:21]
    assert np.allclose(ref, got, rtol=1e-10), np.abs(ref-got).max()
    print("FFT autocov == direct loop  (max dev %.2e)" % np.abs(ref-got).max())

    rows = [one(100+i) for i in range(20)]
    r = np.array([x[0] for x in rows])
    print("\nD_cor estimator on 20 records of a process whose D_cor is exact:")
    print("  ratios: " + " ".join("%.3f" % x for x in r))
    print("  mean %.4f   s.d. %.4f (%.1f%%)   s.e.m %.4f   min %.3f  max %.3f"
          % (r.mean(), r.std(ddof=1), 100*r.std(ddof=1)/r.mean(),
             r.std(ddof=1)/math.sqrt(len(r)), r.min(), r.max()))
    c = np.array([x[1] for x in rows]); v = np.array([x[2] for x in rows])
    print("  cut lag: mean %.1f (%.2f tau)  range %.1f-%.1f" % (c.mean(), c.mean()/tau, c.min(), c.max()))
    print("  var(0)/th: mean %.4f  s.d. %.4f" % (v.mean(), v.std(ddof=1)))
    json.dump(dict(ratios=r.tolist(), cut=c.tolist(), var0=v.tolist(),
                   mean=float(r.mean()), sd=float(r.std(ddof=1))),
              open('estimator_scatter.json', 'w'), indent=1)
