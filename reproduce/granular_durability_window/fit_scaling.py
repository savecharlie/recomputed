"""Two independent routes to the exponent nu in  ln N_f = A + C * alpha^nu / G.

The paper asserts nu = 2 from a mean-field pair argument (a contact is cohesive
only if BOTH grains are coated), and fits p0 = 0.0037, C1 = 26.8 to Fig. 4.
Route 1 is a global fit of (A, C, nu) to every digitized point.
Route 2 never assumes a shared intercept: for each alpha separately, slope of
ln N_f against 1/G is C*alpha^nu, so a log-log fit of those five slopes against
alpha returns nu on its own.  Two routes can disagree; one can only be confirmed.
"""
import json, numpy as np
from scipy.optimize import curve_fit

P = json.load(open("fig3_points.json"))
G = np.array([p["G"] for p in P]); A = np.array([p["alpha"] for p in P])
N = np.array([p["Nf"] for p in P]); L = np.log(N)

# snap G to the discrete drive settings
canon = np.array([0.41, 0.89, 1.40, 2.40, 3.07])
Gs = canon[np.argmin(np.abs(G[:, None] - canon[None, :]), 1)]
print("max |G - snapped| =", np.abs(G - Gs).max().round(3))
for al in sorted(set(A)):
    m = A == al
    print(f"  alpha={al}  n={m.sum():2d}  " +
          "  ".join(f"G={g}:{(Gs[m]==g).sum()}" for g in canon if (Gs[m]==g).any()))

# ---------------- Route 1: global three-parameter fit ------------------------
f = lambda X, A0, C, nu: A0 + C * X[0] ** nu / X[1]
X = np.vstack([A, Gs])
p, cov = curve_fit(f, X, L, p0=[5.6, 26.8, 2.0], maxfev=40000)
e = np.sqrt(np.diag(cov))
res = L - f(X, *p)
print(f"\nRoute 1  nu = {p[2]:.3f} +- {e[2]:.3f}   C = {p[1]:.2f} +- {e[1]:.2f}"
      f"   1/p0 = {np.exp(p[0]):.1f}  (p0 = {np.exp(-p[0]):.5f})")
print(f"         rms residual in ln N_f = {res.std(ddof=3):.3f}"
      f"  ({100*(np.exp(res.std(ddof=3))-1):.0f}% in N_f)")

# nu fixed at 2, for comparison
f2 = lambda X, A0, C: A0 + C * X[0] ** 2 / X[1]
p2, _ = curve_fit(f2, X, L, p0=[5.6, 26.8], maxfev=40000)
r2 = L - f2(X, *p2)
print(f"nu=2 fixed:  C = {p2[1]:.2f}   1/p0 = {np.exp(p2[0]):.1f}"
      f"  (p0 = {np.exp(-p2[0]):.5f})   rms = {r2.std(ddof=2):.3f}")

# ---------------- Route 2: per-alpha slopes, no shared intercept -------------
print("\nRoute 2 — slope of ln N_f vs 1/G, fitted inside each alpha alone")
sl, se, als = [], [], []
for al in sorted(set(A)):
    m = A == al
    x = 1.0 / Gs[m]; y = L[m]
    if len(set(x)) < 2:
        print(f"  alpha={al}: only one G, skipped"); continue
    k, b = np.polyfit(x, y, 1)
    n = len(x); yh = k * x + b
    s = np.sqrt(((y - yh) ** 2).sum() / (n - 2) / ((x - x.mean()) ** 2).sum())
    print(f"  alpha={al}: slope C*alpha^nu = {k:7.2f} +- {s:5.2f}   "
          f"intercept ln(1/p0) = {b:5.2f}   n={n}")
    sl.append(k); se.append(s); als.append(al)
sl = np.array(sl); se = np.array(se); als = np.array(als)
w = 1 / (se / sl) ** 2                      # weights in log space
nu, lnC = np.polyfit(np.log(als), np.log(sl), 1, w=np.sqrt(w))
# bootstrap the interval
rng = np.random.default_rng(587); boot = []
for _ in range(20000):
    s = sl + rng.normal(0, se)
    if (s > 0).all():
        boot.append(np.polyfit(np.log(als), np.log(s), 1, w=np.sqrt(w))[0])
boot = np.array(boot)
print(f"\nRoute 2  nu = {nu:.3f}   68% [{np.percentile(boot,16):.2f}, "
      f"{np.percentile(boot,84):.2f}]   95% [{np.percentile(boot,2.5):.2f}, "
      f"{np.percentile(boot,97.5):.2f}]   C = {np.exp(lnC):.2f}")
json.dump(dict(route1=dict(nu=p[2], nu_err=e[2], C=p[1], p0=float(np.exp(-p[0]))),
               route1_nu2=dict(C=p2[1], p0=float(np.exp(-p2[0]))),
               route2=dict(nu=float(nu), lo68=float(np.percentile(boot,16)),
                           hi68=float(np.percentile(boot,84)),
                           lo95=float(np.percentile(boot,2.5)),
                           hi95=float(np.percentile(boot,97.5)),
                           C=float(np.exp(lnC)),
                           slopes={float(a): float(s) for a, s in zip(als, sl)})),
          open("fit_results.json", "w"), indent=1)
