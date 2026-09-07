"""FT's m_1(t) is CONCAVE on log-log: the local slope falls monotonically from
1.09 to 0.82 across the plotted window (digitize2/curvature read of Fig. 4b).
BN's is straight at about 1.0. Does the IMM at alpha ~ 1 do that, or is the
curvature a deviation from the model?

Same eight sample times as the figure's marker columns, which are equally spaced
in log t (ratio 1.43 per column, starting at day 10).
"""
import numpy as np
from alpha_scan import fit_q, SPECIES, DAYS

R = 400
DAYS_AT = [13.7 * 1.394 ** k for k in range(8)]   # measured off the axis, see axis_calib.txt


def m1_curve(T, alpha, beta, q, cps, rng, R=R):
    acc = np.zeros(len(cps))
    for _ in range(R):
        m = np.zeros(T + 2); m[0] = 1.0; n = 1; ci = 0
        for t in range(1, T + 1):
            if rng.random() < q * n ** (-beta):
                m[n] = 1.0; n += 1
            else:
                w = m[:n] ** alpha
                m[rng.choice(n, p=w / w.sum())] += 1.0
            if ci < len(cps) and t == cps[ci]:
                acc[ci] += m[:n].max(); ci += 1
    return acc / R


for sp, d in SPECIES.items():
    T = int(round(DAYS * d["lam"]))
    q = fit_q(d["beta"], T, d["n200"])
    cps = sorted(set(int(round(dd * d["lam"])) for dd in DAYS_AT if dd * d["lam"] <= T))
    print(f"\n{sp}: beta={d['beta']} q={q:.3f} T={T}; sample steps {cps}")
    for alpha in (0.9, 1.0, 1.1):
        rng = np.random.default_rng(21)
        y = m1_curve(T, alpha, d["beta"], q, cps, rng)
        ls = [np.log(y[i] / y[i - 1]) / np.log(cps[i] / cps[i - 1]) for i in range(1, len(cps))]
        g = np.polyfit(np.log(cps), np.log(y), 1)[0]
        print(f"  alpha {alpha:.2f}  global {g:.3f}   local " +
              " ".join(f"{s:5.3f}" for s in ls), flush=True)
print("\nmeasured off Fig. 4 (top edge of the i=1 marker):")
print("  FT   global 0.958   local 1.091 1.107 0.972 0.911 0.962 0.852 0.818")
print("  BN   global 1.062   local 0.938 0.906 1.178 0.911 1.102 1.330 1.012")
