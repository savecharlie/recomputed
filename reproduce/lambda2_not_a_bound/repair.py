"""Two different things are called "the amplitude", and the paper's sentence is
ambiguous between them.  Measuring both, and testing the obvious repair.

  RELATIVE amplitude   A_rel = max_f  2|W_2(f)| / Var_pss(f)
      "what fraction of this observable's fluctuation is the coherent mode?"
      -> measured: 1 + O(eps^2).  The mode is the WHOLE correlation function
         of the best observable.  NOT vanishing.

  ABSOLUTE amplitude   A_abs = max{ 2|W_2(f)| : f real, RMS(f) = 1 under the
                               UNIFORM (counting) measure on the 3L states }
      "how large is the oscillating signal, for an observable of fixed size that
       does not get to put all its weight on a rare layer?"
      -> this is the one that vanishes, and it vanishes like the stationary
         weight p of the layer the mode lives on.

The difference between the two is entirely the measure the observable is
normalised in: p^ss (which is itself concentrated away from layer A) versus
uniform.  A_rel divides out the very suppression that is the whole mechanism.

The repair the paper's own discussion points at -- weight the bound by how much
of the stationary state the second eigenmode actually occupies -- is therefore
a statement about A_abs, not A_rel.  Tested below:

      sigma_ss  >=  A_abs * lambda_I^2 / lambda_R   ?

WHY THE L2 NORM AND NOT THE SUP NORM.  The first version of this file used the
sup-norm ball |f| <= 1 and restricted the search to the k=1 Fourier sector, on
the grounds that W_2's two linear functionals are supported there.  That is
WRONG, and validation caught it: a full-space vertex search beat the sector
answer by 1.6223x at L=60.  The reason is exactly (4/pi)^2 = 1.6211 -- a square
wave of unit height has a larger fundamental Fourier coefficient (4/pi) than a
cosine of unit height, W_2 is a PRODUCT of two such functionals, and the sup-norm
ball is not Fourier-diagonal so the maximiser leaves the sector.  Under the
UNIFORM L2 norm Parseval makes the sector restriction exact, so the 3x3 problem
is the whole problem.  The sup-norm ratio is reported below as a separate
measurement, because its value is a fact about square waves and not about this
model.
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "2")
import sys, json, numpy as np
sys.path.insert(0, ".")
import model as M
import amplitude as A
from scipy.optimize import minimize

FAIL = []
def chk(n, c, d=""):
    print(("  PASS " if c else "  FAIL ") + n + ("   " + d if d else ""))
    if not c: FAIL.append(n)


def W2_sector(P, g, s=None):
    phi, psi, pi, lam = s if s else A.sector(P)
    return 0.25 * (phi @ np.conj(g)) * ((psi * pi) @ g)

def A_abs(P, n=40, seed=5):
    """max 2|W_2(f)| over real f with RMS(f)=1 in the uniform measure.

    f_{x,a} = Re[g_a e^{2 pi i x/L}] has sum_{x,a} f^2 = (L/2) sum_a |g_a|^2, so
    RMS 1 over 3L states means sum_a |g_a|^2 = 6.  Parseval makes the k=1 sector
    restriction exact for this norm."""
    s = A.sector(P)
    rng = np.random.default_rng(seed)
    def neg(z):
        g = z[:3] + 1j * z[3:]
        nrm = np.linalg.norm(g)
        if nrm < 1e-14:
            return 0.0
        g = g * (np.sqrt(6.0) / nrm)
        return -2.0 * abs(W2_sector(P, g, s))
    vals = []
    for _ in range(n):
        r = minimize(neg, rng.standard_normal(6), method="Nelder-Mead",
                     options=dict(maxiter=4000, xatol=1e-13, fatol=1e-15))
        vals.append(-r.fun)
    return float(np.max(vals)), np.array(vals)

def A_abs_full_L2(P, n=30, seed=9):
    """Same quantity under the uniform L2 norm, in the FULL 3L space, by
    projected gradient-free ascent.  Must not beat the sector answer."""
    R, w, V, Li = A.modes_full(P)
    p = M.pss(P); L = P["L"]
    a_vec, b_vec = V[:, 1], Li[1, :] * p
    rng = np.random.default_rng(seed)
    best = 0.0
    for _ in range(n):
        def neg(f):
            nrm = np.linalg.norm(f)
            if nrm < 1e-14: return 0.0
            f = f * (np.sqrt(3.0 * L) / nrm)
            return -2.0 * abs((a_vec @ f) * (b_vec @ f))
        r = minimize(neg, rng.standard_normal(3 * L), method="Powell",
                     options=dict(maxiter=20000, xtol=1e-10, ftol=1e-12))
        best = max(best, -r.fun)
    return best


def A_abs_sup_full(P, n=400, seed=9):
    """The SUP-NORM version, full space, vertex search.  Reported only as the
    measurement that killed the sector restriction for that norm."""
    R, w, V, Li = A.modes_full(P)
    p = M.pss(P); L = P["L"]
    a_vec = V[:, 1]
    b_vec = Li[1, :] * p
    rng = np.random.default_rng(seed)
    best = 0.0
    for _ in range(n):
        f = rng.choice([-1.0, 1.0], size=3 * L)
        cur = 2 * abs((a_vec @ f) * (b_vec @ f))
        improved = True
        while improved:
            improved = False
            for i in rng.permutation(3 * L):
                f[i] = -f[i]
                nv = 2 * abs((a_vec @ f) * (b_vec @ f))
                if nv > cur + 1e-15:
                    cur = nv; improved = True
                else:
                    f[i] = -f[i]
        best = max(best, cur)
    return best


def sigma_parts(P):
    """sigma_ss split the paper's way: (i) the driven ring in layer A, and
    (ii) the inter-layer circulation A->B->C->A, summed over sites."""
    L, p = P["L"], P["p"]
    cell = np.array([p, (1 - p) / 2, (1 - p) / 2]) / L
    rp, rm = P["rA_p"], P["rA_m"]
    sig_A = L * (rp - rm) * cell[0] * np.log(rp / rm)
    Q = M.Q_matrix(P); sig_in = 0.0
    for a in range(3):
        for b in range(3):
            if a != b:
                f, g = Q[a, b] * cell[b], Q[b, a] * cell[a]
                if f > 0: sig_in += L * f * np.log(f / g)
    return sig_A, sig_in


if __name__ == "__main__":
    print(__doc__)
    print("VALIDATION  under the UNIFORM L2 norm the k=1 restriction is exact")
    for Lv in (60, 100):
        P = M.rates(Lv, q=0.55)
        if not M.valid(P)[0]:
            print(f"  (L={Lv} outside the model, skipped)"); continue
        s_sec = A_abs(P)[0]; s_full = A_abs_full_L2(P)
        chk(f"L={Lv} full-space search does not beat the sector",
            s_full <= s_sec * (1 + 1e-5),
            f"sector {s_sec:.8e}  full {s_full:.8e}  ratio {s_full/s_sec:.8f}")

    print("\nMEASUREMENT  for the SUP norm the sector is NOT exact; the factor is")
    print("             a fact about square waves, (4/pi)^2 = 1.621139, not about")
    print("             this model.  This is why the L2 norm is used above.")
    for Lv in (40, 60, 80):
        P = M.rates(Lv, q=0.55)
        if not M.valid(P)[0]: continue
        # sector value under the SUP ball: |g_a| <= 1
        s = A.sector(P); rng2 = np.random.default_rng(3)
        def negsup(z):
            g = z[:3] + 1j * z[3:]
            g = g / np.maximum(np.abs(g), 1.0)
            return -2.0 * abs(W2_sector(P, g, s))
        sec = max(-minimize(negsup, rng2.standard_normal(6), method="Nelder-Mead",
                  options=dict(maxiter=4000)).fun for _ in range(30))
        full = A_abs_sup_full(P, n=50)
        print(f"  L={Lv:3d}  sup-norm sector {sec:.6e}  full {full:.6e}"
              f"  ratio {full/sec:.6f}   ((4/pi)^2 = 1.621139)")
    print("\nVALIDATION  sigma_A + sigma_in == sigma_ss (the paper's own split)")
    worst = 0
    for Lv, q in ((100, 0.55), (1000, 0.55), (1559, 0.75)):
        P = M.rates(Lv, q=q)
        sa, si = sigma_parts(P); tot = M.sigma_ss_fast(P)
        worst = max(worst, abs(sa + si - tot) / tot)
    chk("split sums to the total at 3 sizes", worst < 1e-12, f"max rel {worst:.1e}")
    if FAIL:
        print(f"\nVALIDATION FAILED {FAIL}"); sys.exit(1)

    print("\n" + "=" * 92)
    print("A_rel vs A_abs, and the repaired bound, at q = 0.55")
    print("=" * 92)
    print(f"{'L':>7} {'eps':>10} {'p':>10} {'A_rel':>10} {'A_abs':>11} {'A_abs/p':>8}"
          f" | {'sig_A':>10} {'sig_in':>10} {'sigma':>10} {'rhs':>7}"
          f" {'A_abs*rhs':>10} {'margin':>9}")
    rows = []
    for Lv in (100, 200, 400, 800, 1600, 3200, 6400):
        P = M.rates(Lv, q=0.55)
        ar = A.best_A(P, n=25)[0]
        aa = A_abs(P)[0]
        sa, si = sigma_parts(P)
        r = M.osb(P)
        rep = aa * r["rhs"]
        rows.append(dict(L=Lv, eps=P["eps"], p=P["p"], A_rel=ar, A_abs=aa,
                         sig_A=sa, sig_in=si, sigma=r["sigma"], rhs=r["rhs"],
                         repaired=rep, margin=r["sigma"] / rep))
        print(f"{Lv:7d} {P['eps']:10.3e} {P['p']:10.3e} {ar:10.6f} {aa:11.4e}"
              f" {aa/P['p']:8.4f} | {sa:10.3e} {si:10.3e} {r['sigma']:10.3e}"
              f" {r['rhs']:7.4f} {rep:10.3e} {r['sigma']/rep:9.1f}x")
    print()
    e = np.array([r["eps"] for r in rows])
    print("SCALING (power-law fit in eps over the seven sizes):")
    for k in ("A_rel", "A_abs", "p", "sig_A", "sig_in", "sigma", "repaired", "margin"):
        y = np.array([r[k] for r in rows])
        sl = np.polyfit(np.log(e), np.log(y), 1)[0]
        print(f"  {k:9s} ~ eps^{sl:+.4f}")
    arel = np.array([r["A_rel"] for r in rows])
    sl = np.polyfit(np.log(e), np.log(arel - 1.0), 1)[0]
    print(f"  A_rel - 1 ~ eps^{sl:+.4f}   (A_rel approaches 1 FROM ABOVE: the mode"
          f" weights are not all positive in a non-reversible chain)")
    print()
    ok = all(r["margin"] > 1 for r in rows)
    print(f"REPAIRED BOUND  sigma_ss >= A_abs * lambda_I^2/lambda_R :"
          f"  {'HOLDS' if ok else 'VIOLATED'} at all seven sizes,"
          f" margin {min(r['margin'] for r in rows):.0f}x to"
          f" {max(r['margin'] for r in rows):.0f}x and GROWING like"
          f" eps^{np.polyfit(np.log(e), np.log([r['margin'] for r in rows]), 1)[0]:+.2f}")
    json.dump(rows, open("repair.json", "w"), indent=1)
    print("wrote repair.json")
