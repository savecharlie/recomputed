"""The amplitude the paper asserts but never computes.

arXiv:2609.34352v1 explains its counterexample like this:

  "the slowest relaxation of the time correlation function is described as
   A e^{-lambda_R t} sin(lambda_I t) + ...  However, this description does not
   guarantee that its amplitude A is nonnegligible, and our counterexample
   indeed has vanishing amplitude A on this mode."

A is never defined and never evaluated.  Definition used here, chosen to be
exactly the A in that sentence, normalised by the t=0 value of the same
correlation function so it is dimensionless and observable-comparable:

    C_f(t) = <f(t) f(0)> - <f>^2 = sum_{n>=1} e^{lambda_n t} W_n(f)
    W_n(f) = (f . r_n)(l_n . (f * p^ss)),   R r_n = lam_n r_n,  l_n R = lam_n l_n,
                                            l_m . r_n = delta_mn
    C_f(0) = Var_pss(f) = sum_{n>=1} W_n(f)                      [sum rule]

the lambda_2 conjugate pair contributes 2|W_2| e^{-lam_R t} cos(lam_I t + arg W_2),
hence

    A(f) := 2 |W_2(f)| / Var_pss(f)   in [0, 1]-ish, = 1 if the mode is the
                                      whole correlation function.

THE k=1 SECTOR.  R is shift invariant, so only the k=+-1 Fourier component of f
couples to the lambda_2 mode (which lives at k=1).  Writing

    f_{x,a} = Re[ g_a e^{2 pi i x / L} ],   g in C^3,

and with phi, psi the right/left 3-vector eigenvectors of the block R-hat^1 for
lambda_2 (psi . phi = 1) and pi = (p, (1-p)/2, (1-p)/2) the per-layer stationary
weights, the algebra collapses to

    W_2 = (1/4) (sum_a phi_a conj(g_a)) (sum_a psi_a pi_a g_a)
    Var = (1/2) sum_a pi_a |g_a|^2
    A   = |sum_a phi_a conj(g_a)| |sum_a psi_a pi_a g_a| / sum_a pi_a |g_a|^2

-- three complex dimensions at ANY L.  Checked against the full 3L-dimensional
computation below before it is used for anything.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import sys, numpy as np
sys.path.insert(0, ".")
import model as M
from scipy.linalg import expm
from scipy.optimize import minimize

FAIL = []
def chk(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + ("   " + detail if detail else ""))
    if not cond: FAIL.append(name)


# ---------- full space (validation only; O((3L)^3)) --------------------------

def modes_full(P):
    R = M.R_full(P)
    w, V = np.linalg.eig(R)
    Li = np.linalg.inv(V)
    o = np.argsort(-w.real)
    return R, w[o], V[:, o], Li[o, :]

def A_full(P, f, cache={}):
    key = (P["L"], P["eps"], P["d"], P["kappa"], P["K"])
    if key not in cache:
        cache.clear(); cache[key] = modes_full(P)
    R, w, V, Li = cache[key]
    p = M.pss(P)
    W = (V.T @ f) * (Li @ (f * p))
    var = float(p @ (f * f) - (p @ f) ** 2)
    sumrule = abs(W[1:].sum().real - var) / max(var, 1e-300)
    return 2.0 * abs(W[1]) / var, var, sumrule, w[1]


# ---------- k=1 sector (the instrument; 3x3 at any L) -----------------------

def sector(P):
    """(phi, psi, pi, lam2) for the lambda_2 mode of the k=1 block."""
    B = M.R_block(P, 1)
    w, V = np.linalg.eig(B)
    j = int(np.argmax(w.real))                  # largest real part within k=1
    phi = V[:, j]
    psi = np.linalg.inv(V)[j, :]                # psi . phi = 1 by construction
    p = P["p"]
    pi = np.array([p, (1 - p) / 2, (1 - p) / 2])
    return phi, psi, pi, w[j]

def A_sector(P, g, s=None):
    phi, psi, pi, lam = s if s else sector(P)
    # A = 2|W_2|/Var with W_2 = (1/4) N and Var = (1/2) D  ==>  A = N/D exactly.
    # (The first version of this line carried a spurious 1/2 and validation 2
    #  caught it as a clean factor of two at 18/18 random g.  Kept in the note
    #  because the factor is the whole content of the formula.)
    num = abs(phi @ np.conj(g)) * abs((psi * pi) @ g)
    den = float(np.real(pi @ (np.abs(g) ** 2)))
    return num / den

def best_A(P, n=60, seed=11):
    """max_g A(g) over C^3.  Many starts; the spread is reported, never hidden."""
    s = sector(P)
    rng = np.random.default_rng(seed)
    f = lambda z: -A_sector(P, z[:3] + 1j * z[3:], s)
    vals = []
    for _ in range(n):
        z0 = rng.standard_normal(6)
        r = minimize(f, z0, method="Nelder-Mead",
                     options=dict(maxiter=4000, xatol=1e-13, fatol=1e-15))
        vals.append(-r.fun)
    return float(np.max(vals)), np.array(vals)


if __name__ == "__main__":
    print(__doc__)
    print("VALIDATION 1  sum rule and the propagator, in the full space")
    P = M.rates(100, q=0.55)
    L = P["L"]; rng = np.random.default_rng(1)
    rs = [A_full(P, rng.standard_normal(3 * L))[2] for _ in range(5)]
    chk("sum_n W_n(f) == Var(f) on 5 random f", max(rs) < 1e-9, f"max rel {max(rs):.2e}")
    R, w, V, Li = modes_full(P); p = M.pss(P)
    f = rng.standard_normal(3 * L)
    W = (V.T @ f) * (Li @ (f * p))
    err = 0.0
    for t in (0.0, 0.1, 0.5, 1.0, 3.0):
        direct = float(f @ (expm(R * t) @ (f * p))) - float(p @ f) ** 2
        spect = float(np.real((np.exp(w[1:] * t) * W[1:]).sum()))
        err = max(err, abs(direct - spect) / max(abs(direct), 1e-12))
    chk("C_f(t) from expm(Rt) == C_f(t) from the modes, 5 times t", err < 1e-8,
        f"max rel {err:.2e}")

    print("\nVALIDATION 2  the 3x3 sector formula == the full 3L computation")
    worst = 0.0
    for Lv, q in ((60, 0.55), (100, 0.55), (160, 0.6)):
        P = M.rates(Lv, q=q)
        s = sector(P)
        ph = np.exp(2j * np.pi * np.arange(Lv) / Lv)
        for _ in range(6):
            g = rng.standard_normal(3) + 1j * rng.standard_normal(3)
            f = np.zeros(3 * Lv)
            for a in range(3):
                f[a::3] = np.real(g[a] * ph)
            af = A_full(P, f)[0]
            asec = A_sector(P, g, s)
            worst = max(worst, abs(af - asec) / max(af, 1e-300))
    chk("18 random g at 3 sizes agree", worst < 1e-7, f"max rel {worst:.2e}")

    print("\nVALIDATION 3  lambda_2 of the k=1 block IS lambda_2 of R")
    bad = []
    for Lv, q in ((60, 0.55), (100, 0.55), (1559, 0.75), (10000, 0.6)):
        P = M.rates(Lv, q=q)
        lam_sec = sector(P)[3]
        lam_all = M.osb(P)["lam2"]
        if abs(lam_sec - lam_all) > 1e-8 and abs(np.conj(lam_sec) - lam_all) > 1e-8:
            bad.append((Lv, lam_sec, lam_all))
    chk("the second-largest eigenvalue of R lives in the k=1 sector", not bad, str(bad))

    if FAIL:
        print(f"\nVALIDATION FAILED: {FAIL} -- reporting nothing."); sys.exit(1)
    print("\nvalidation complete.\n")

    print("=" * 78)
    print("A, THE AMPLITUDE, for three observables and for the best possible one")
    print("=" * 78)
    print("  ring   : f = cos(2 pi x / L), the same in all three layers")
    print("  ringA  : f = cos(2 pi x / L) in layer A only, 0 in B and C")
    print("  best   : max over all observables\n")
    print(f"{'L':>7} {'eps':>10} {'p':>11} {'sigma_ss':>10} | {'A(ring)':>11}"
          f" {'A(ringA)':>10} {'A(best)':>9} {'spread':>8}")
    rows = []
    for Lv in (100, 200, 400, 800, 1600, 3200):
        P = M.rates(Lv, q=0.55)
        s = sector(P)
        g_ring = np.array([1.0, 1.0, 1.0]); g_A = np.array([1.0, 0.0, 0.0])
        a1 = A_sector(P, g_ring, s); a2 = A_sector(P, g_A, s)
        ab, vals = best_A(P)
        r = M.osb(P)
        rows.append(dict(L=Lv, eps=P["eps"], p=P["p"], sigma=r["sigma"],
                         ring=a1, ringA=a2, best=ab, rhs=r["rhs"]))
        print(f"{Lv:7d} {P['eps']:10.3e} {P['p']:11.3e} {r['sigma']:10.3e} |"
              f" {a1:11.4e} {a2:10.6f} {ab:9.6f} {vals.std():8.1e}")

    print("\nSCALING in eps (power-law fit over the six sizes):")
    e = np.array([r["eps"] for r in rows])
    for k in ("ring", "ringA", "best", "sigma", "p"):
        y = np.array([r[k] if k != "p" else r["p"] for r in rows])
        sl = np.polyfit(np.log(e), np.log(y), 1)[0]
        print(f"  {k:6s} ~ eps^{sl:+.4f}")
    import json; json.dump(rows, open("amplitude.json", "w"), indent=1)
    print("\nwrote amplitude.json")
