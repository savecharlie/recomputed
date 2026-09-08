"""Is ln N_f really a function of the single variable alpha^2/G?

Eq. 2 of the paper is a PRODUCT form: one coefficient multiplying alpha^2/G.
That is a strong claim -- it says alpha and 1/G enter only through their product.
Fit three nested models to the Fig. 3 data (digitized) and compare:
   M1  A + P*(a^2/G)                 the paper's form
   M2  A + B*a^2 + C*(1/G)           additive, no interaction at all
   M3  A + B*a^2 + C*(1/G) + D*(a^2/G)
Also fit M1 on the Fig-4 window (G >= 2.4) alone, to see whether the paper's
C1 = 26.8 is recovered there and nowhere else.
"""
import json, numpy as np
P = json.load(open("fig3_points.json"))
canon = np.array([0.41, 0.89, 1.40, 2.40, 3.07])
al = np.array([p["alpha"] for p in P])
G  = canon[np.argmin(np.abs(np.array([p["G"] for p in P])[:, None] - canon), 1)]
L  = np.log(np.array([p["Nf"] for p in P]))
a2 = al ** 2; ig = 1.0 / G; pr = a2 * ig
n  = len(L)

def ols(cols, names, mask=None):
    m = np.ones(n, bool) if mask is None else mask
    X = np.column_stack([np.ones(m.sum())] + [c[m] for c in cols])
    y = L[m]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    k = X.shape[1]; N = m.sum()
    s2 = r @ r / (N - k)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    aic = N * np.log(r @ r / N) + 2 * k
    ss_tot = ((y - y.mean()) ** 2).sum()
    return beta, se, np.sqrt(r @ r / (N - k)), aic, 1 - (r @ r) / ss_tot, names

def show(tag, res):
    beta, se, rmse, aic, r2, names = res
    print(f"{tag}")
    for nm, b, s in zip(["const"] + names, beta, se):
        print(f"    {nm:>12} = {b:9.3f} +- {s:6.3f}")
    print(f"    rmse(lnNf) {rmse:.3f}   R2 {r2:.3f}   AIC {aic:.1f}\n")

print(f"n = {n} digitized points, G in [{G.min()}, {G.max()}]\n")
show("M1  paper's product form, ALL of Fig 3", ols([pr], ["C1 (a^2/G)"]))
show("M2  additive, no interaction",           ols([a2, ig], ["B (a^2)", "C (1/G)"]))
show("M3  full",                               ols([a2, ig, pr], ["B (a^2)", "C (1/G)", "D (a^2/G)"]))

hi = G >= 2.4
print(f"--- restricted to the Fig-4 window, G >= 2.4  (n = {hi.sum()}) ---")
show("M1 on G >= 2.4 only", ols([pr], ["C1 (a^2/G)"], hi))
lo = ~hi
print(f"--- the data Fig 4 leaves out, G <= 1.4  (n = {lo.sum()}) ---")
show("M1 on G <= 1.4 only", ols([pr], ["C1 (a^2/G)"], lo))
