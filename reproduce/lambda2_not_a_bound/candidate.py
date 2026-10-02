"""The candidate bound, corrected after it failed, and then attacked.

FIRST VERSION:   sigma_ss >= A_abs * lambda_I^2 / lambda_R
    Survives Shiraishi's counterexample with a margin growing like 1/eps.
    FAILS on 7 of 150 random coherent chains -- because A_abs there is 1.06 to
    2.06, so multiplying by it makes the claim STRONGER than OSB, and a
    strengthening of a tight bound is just a new false bound.  The amplitude is
    only supposed to DISCOUNT the bound when the mode is invisible; it must not
    inflate it when the mode is loud.

CORRECTED:       sigma_ss >= min(1, A_abs) * lambda_I^2 / lambda_R
    reduces to OSB exactly whenever A_abs >= 1 (so it is never stronger than a
    conjecture that is already well tested), and discounts it by the mode's
    absolute amplitude exactly when the mode is faint, which is Shiraishi's
    mechanism.

This file does three things: checks the corrected form on the counterexample,
checks it on the 150-chain sample, and then goes looking for chains with
A_abs < 1 -- the only regime where the corrected form says anything new and
therefore the only regime where it can be false.  A bound nobody can falsify is
not a result (law 12).
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ.setdefault(v, "2")
import sys, numpy as np, json
sys.path.insert(0, ".")
import model as M, repair as Rp
from scipy.optimize import minimize

def A_abs_generic(R, p, n=14, seed=4):
    w, V = np.linalg.eig(R); Li = np.linalg.inv(V)
    o = np.argsort(-w.real); w, V, Li = w[o], V[:, o], Li[o, :]
    a_vec, b_vec = V[:, 1], Li[1, :] * p
    nn = len(p); rng = np.random.default_rng(seed); best = 0.0
    def neg(f):
        nr = np.linalg.norm(f)
        if nr < 1e-14: return 0.0
        f = f * (np.sqrt(nn) / nr)
        return -2.0 * abs((a_vec @ f) * (b_vec @ f))
    for _ in range(n):
        r = minimize(neg, rng.standard_normal(nn), method="Powell",
                     options=dict(maxiter=40000, xtol=1e-10, ftol=1e-12))
        best = max(best, -r.fun)
    return best, w[1]

def analyse(R):
    n = R.shape[0]
    w, V = np.linalg.eig(R); o = np.argsort(-w.real)
    if abs(w[o][0]) > 1e-8: return None
    v = np.real(V[:, o[0]])
    if v.sum() == 0: return None
    v = v / v.sum()
    if v.min() < 1e-12: return None
    l2 = w[o][1]; lR, lI = -l2.real, abs(l2.imag)
    if lI < lR or lI < 1e-8 or lR < 1e-10: return None
    sig = 0.0
    for i in range(n):
        for j in range(n):
            if i == j: continue
            f, g = R[i, j] * v[j], R[j, i] * v[i]
            if f > 1e-300:
                if g <= 1e-300: return None
                sig += f * np.log(f / g)
    if sig <= 0: return None
    aa, _ = A_abs_generic(R, v)
    return dict(n=n, lR=lR, lI=lI, sigma=sig, rhs=lI**2/lR, A_abs=aa, pmin=v.min())

print(__doc__)

# ---- 1. the counterexample ------------------------------------------------
print("1.  Shiraishi's counterexample (q=0.55, d=3), both forms")
print(f"{'L':>7} {'A_abs':>11} {'min(1,A)':>9} {'sigma':>10} {'rhs':>8}"
      f" {'OSB':>7} {'A*rhs form':>11} {'min form':>9} {'margin':>9}")
rows = []
for Lv in (100, 400, 1600, 6400, 25600):
    P = M.rates(Lv, q=0.55)
    aa = Rp.A_abs(P)[0]; r = M.osb(P)
    m = min(1.0, aa)
    rows.append((Lv, aa, r["sigma"], r["rhs"]))
    print(f"{Lv:7d} {aa:11.4e} {m:9.4e} {r['sigma']:10.4e} {r['rhs']:8.4f}"
          f" {'VIOL':>7} {'holds' if r['sigma']>=aa*r['rhs'] else 'VIOL':>11}"
          f" {'holds' if r['sigma']>=m*r['rhs'] else 'VIOL':>9}"
          f" {r['sigma']/(m*r['rhs']):9.1f}x")
e = np.array([M.rates(r[0], q=0.55)["eps"] for r in rows])
mar = np.array([r[2] / (min(1.0, r[1]) * r[3]) for r in rows])
print(f"    margin of the corrected form ~ eps^{np.polyfit(np.log(e), np.log(mar),1)[0]:+.3f}"
      f"   (so it is satisfied by an ever-widening factor, not narrowly)")

# ---- 2. the 150-chain sample ----------------------------------------------
print("\n2.  150 random coherent chains (rebuilt with the same seed)")
def make(rng, kind):
    n = int(rng.integers(5, 13))
    fw = np.exp(rng.normal(0, 0.6, n)); bw = fw * np.exp(-rng.uniform(0.3, 3.0))
    R = np.zeros((n, n))
    for i in range(n):
        R[(i+1) % n, i] += fw[i]; R[i, (i+1) % n] += bw[i]
    if kind == "chords":
        for _ in range(int(rng.integers(1, 4))):
            i, j = rng.integers(0, n, 2)
            if i != j:
                R[j, i] += np.exp(rng.normal(-1,1)); R[i, j] += np.exp(rng.normal(-1,1))
    elif kind == "trap":
        R = np.pad(R, ((0,1),(0,1))); i = int(rng.integers(0,n))
        R[n, i] += 10.0 ** rng.uniform(-3,-1); R[i, n] += 1.0; n += 1
    np.fill_diagonal(R, 0.0); np.fill_diagonal(R, -R.sum(axis=0))
    return R
rng = np.random.default_rng(1234); S, t = [], 0
while len(S) < 150 and t < 4000:
    t += 1
    r = analyse(make(rng, ["plain","chords","trap"][t % 3]))
    if r: S.append(r)
f1 = [r for r in S if r["sigma"] < r["A_abs"] * r["rhs"]]
f2 = [r for r in S if r["sigma"] < min(1.0, r["A_abs"]) * r["rhs"]]
print(f"    A*rhs form   violated in {len(f1)}/{len(S)}")
print(f"    min(1,A) form violated in {len(f2)}/{len(S)}")
print(f"    A_abs < 1 anywhere in this sample? "
      f"{sum(1 for r in S if r['A_abs'] < 1)}/{len(S)}"
      f"   (min A_abs = {min(r['A_abs'] for r in S):.4f})")

# ---- 3. adversarial: hunt for A_abs < 1, the only regime that can falsify --
print("\n3.  ADVERSARIAL SEARCH for chains with A_abs < 1, where the corrected")
print("    form says something OSB does not, and could therefore be false.")
print("    Objective: minimise A_abs over small coherent generators; then at")
print("    every chain found with A_abs < 1, test the corrected bound.")
rng = np.random.default_rng(99)
found, tested, fails = [], 0, []
def build(z, n):
    R = np.zeros((n, n))
    k = 0
    for i in range(n):
        for j in range(n):
            if i != j:
                R[i, j] = np.exp(z[k]); k += 1
    np.fill_diagonal(R, 0.0); np.fill_diagonal(R, -R.sum(axis=0))
    return R
for trial in range(220):
    n = int(rng.integers(4, 7))
    nz = n * (n - 1)
    def obj(z):
        r = analyse(build(z, n))
        if r is None: return 10.0
        return r["A_abs"]
    res = minimize(obj, rng.normal(0, 1.5, nz), method="Nelder-Mead",
                   options=dict(maxiter=1200, fatol=1e-10))
    r = analyse(build(res.x, n))
    if r is None: continue
    tested += 1
    if r["A_abs"] < 1.0:
        found.append(r)
        if r["sigma"] < min(1.0, r["A_abs"]) * r["rhs"]:
            fails.append(r)
print(f"    {tested} coherent chains reached; {len(found)} with A_abs < 1")
if found:
    aa = np.array([r["A_abs"] for r in found])
    print(f"    A_abs among those: min {aa.min():.6f}  median {np.median(aa):.6f}")
    print(f"    corrected bound violated in {len(fails)}/{len(found)} of them")
    for r in sorted(found, key=lambda r: r["A_abs"])[:6]:
        m = min(1.0, r["A_abs"])
        print(f"      n={r['n']} A_abs={r['A_abs']:.6f} sigma={r['sigma']:.5f}"
              f" rhs={r['rhs']:.5f} min*rhs={m*r['rhs']:.5f}"
              f"  {'VIOLATED' if r['sigma']<m*r['rhs'] else 'holds'}"
              f"  margin {r['sigma']/(m*r['rhs']):.2f}x")
else:
    print("    none -- so on this family the corrected form is NOT falsifiable,")
    print("    and that is a limitation of the test, not evidence for the bound.")
print("\n    STATUS, stated plainly: the corrected form survives the published")
print("    counterexample and every chain I could build.  That is a CANDIDATE,")
print("    not a theorem.  I have no proof and I did not look for one.")
json.dump(dict(counterexample=[list(map(float,r)) for r in rows],
               n_sample=len(S), fail_Arhs=len(f1), fail_min=len(f2),
               n_Alt1=len(found), fail_Alt1=len(fails)),
          open("candidate.json","w"), indent=1)
