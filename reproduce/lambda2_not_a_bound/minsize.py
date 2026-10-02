"""How small can the counterexample be?

The paper's violation is asymptotic in L and no numbers are printed at any
finite size.  At a fixed L the construction has four free parameters
(eps, d, kappa, K) -- q appears nowhere except through eps = L^-q -- so the
question "what is the smallest Markov chain that breaks the OSB conjecture"
is a bounded four-parameter maximisation of

    ratio(L; eps, d, kappa, K) = (lambda_I^2 / lambda_R) / sigma_ss

subject to  all rates >= 0  and the conjecture's own precondition lambda_I >= lambda_R.

Constraints, exactly:
  r^-_A > 0           <=>  eps > 2 tan(pi/L)                  (hard floor on eps)
  Q off-diagonals >=0 <=>  kappa >= kappa_min(eps,d)  and  h >= v
  d > 2               so that lambda_+- of Q is complex and sits below -1
  eps < 1             so that a = 1-eps > 0
"""
import sys, json, numpy as np
sys.path.insert(0, ".")
import model as M
from scipy.optimize import minimize

RNG = np.random.default_rng(20261002)      # fixed seed: byte-determinism


def ratio_at(L, eps, d, kmul, K):
    """-ratio, for a minimiser.  Returns +inf outside the feasible set."""
    if not (0.0 < eps < 1.0) or d <= 2.0 or kmul < 1.0 or K <= 0.0:
        return np.inf, None
    if eps <= 2.0 * np.tan(np.pi / L):
        return np.inf, None
    P = M.rates_eps(L, eps, d=d, kappa=None, K=K)
    P = M.rates_eps(L, eps, d=d, kappa=kmul * P["kappa_min"], K=K)
    if P["p"] >= 1.0:
        return np.inf, None
    ok, _ = M.valid(P)
    if not ok:
        return np.inf, None
    r = M.osb(P)
    if not r["precondition"] or not np.isfinite(r["ratio"]):
        return np.inf, None
    return -r["ratio"], (P, r)


def best_at(L, n_starts=40):
    """Maximise the violation ratio at this L.  Coarse random starts + Nelder-Mead."""
    lo = 2.0 * np.tan(np.pi / L)
    if lo >= 1.0:
        return None                       # no admissible eps at all
    best = (np.inf, None, None)
    starts = []
    for _ in range(n_starts):
        e = lo + (1.0 - lo) * RNG.random() ** 2        # bias toward the floor
        starts.append([e, 2.0 + 8.0 * RNG.random(), 1.0 + 2.0 * RNG.random(),
                       10.0 ** (RNG.random() * 3)])
    starts.append([min(0.999, lo * 1.001), 3.0, 1.05, 10.0])
    starts.append([min(0.999, (lo + 1.0) / 2), 3.0, 1.05, 10.0])
    for x0 in starts:
        f = lambda x: ratio_at(L, *x)[0]
        res = minimize(f, x0, method="Nelder-Mead",
                       options=dict(maxiter=600, xatol=1e-10, fatol=1e-12))
        val, pack = ratio_at(L, *res.x)
        if val < best[0]:
            best = (val, res.x, pack)
    if best[2] is None:
        return None
    return dict(L=L, ratio=-best[0], x=list(best[1]), P=best[2][0], r=best[2][1])


if __name__ == "__main__":
    print(__doc__)
    print(f"{'L':>5} {'states':>7} {'eps_floor':>10} {'best eps':>10} {'d':>7}"
          f" {'kappa/min':>9} {'K':>9} {'sigma_ss':>11} {'rhs':>8} {'max ratio':>10}")
    first = None
    out = []
    for L in list(range(7, 60)) + list(range(60, 121, 2)):
        b = best_at(L)
        if b is None:
            print(f"{L:5d} {3*L:7d}  -- no admissible point (eps floor "
                  f"{2*np.tan(np.pi/L):.4f} >= 1)")
            continue
        P, r = b["P"], b["r"]
        out.append(dict(L=L, ratio=b["ratio"], eps=P["eps"], d=P["d"],
                        kappa=P["kappa"], K=P["K"], sigma=r["sigma"],
                        rhs=r["rhs"], q=P["q"]))
        flag = "  <-- VIOLATED" if b["ratio"] > 1 else ""
        print(f"{L:5d} {3*L:7d} {2*np.tan(np.pi/L):10.4f} {P['eps']:10.5f}"
              f" {P['d']:7.3f} {b['x'][2]:9.4f} {P['K']:9.3f} {r['sigma']:11.4e}"
              f" {r['rhs']:8.4f} {b['ratio']:10.4f}{flag}")
        if b["ratio"] > 1 and first is None:
            first = b
    print()
    if first:
        P, r = first["P"], first["r"]
        print("=" * 76)
        print(f"SMALLEST VIOLATING SYSTEM FOUND: L = {P['L']}  ->  {3*P['L']} states")
        print("=" * 76)
        for k in ("eps", "d", "kappa", "K", "p", "a", "u", "h", "v", "b",
                  "rA_p", "rA_m", "rBC", "q"):
            print(f"  {k:6s} = {P[k]!r}")
        print(f"  lambda_2       = {r['lam2']}")
        print(f"  lambda_R       = {r['lR']:.10f}")
        print(f"  lambda_I       = {r['lI']:.10f}   (precondition lambda_I>=lambda_R: "
              f"{r['precondition']})")
        print(f"  sigma_ss       = {r['sigma']:.10f}")
        print(f"  lambda_I^2/lR  = {r['rhs']:.10f}")
        print(f"  VIOLATION      = {r['ratio']:.6f}x")
    json.dump(out, open("minsize.json", "w"), indent=1)
    print("\nwrote minsize.json")
