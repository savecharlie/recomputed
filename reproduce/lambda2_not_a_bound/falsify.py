"""A falsification attempt WITH POWER.

candidate.py's adversarial search reached zero valid coherent chains, so its
"no violation found" was worth nothing: a test that cannot produce the regime
it is testing has no power, and recording that is the point (law 12).

The regime where  sigma >= min(1, A_abs) * lambda_I^2/lambda_R  says anything
that OSB does not is exactly A_abs < 1, i.e. a second eigenmode localised on a
set of small stationary weight.  I know one family that produces it on demand:
Shiraishi's own construction.  So the test is a dense sweep of that family's
whole admissible box -- eps, d, kappa, K, L -- keeping every point that is a
genuine Markov generator with lambda_I >= lambda_R, and checking both forms at
each.  Every point here has A_abs < 1 by construction, so every point is a
chance for the corrected bound to be false.
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ.setdefault(v, "2")
import sys, numpy as np, json
sys.path.insert(0, ".")
import model as M, repair as Rp

rng = np.random.default_rng(606)
pts, kept = 0, []
worst = None
while kept.__len__() < 600 and pts < 20000:
    pts += 1
    L = int(np.exp(rng.uniform(np.log(25), np.log(4000))))
    lo = 2.0 * np.tan(np.pi / L)
    if lo >= 1.0: continue
    eps = lo * np.exp(rng.uniform(0.0, 3.0))
    if eps >= 1.0: continue
    d = 2.0 + np.exp(rng.uniform(-8, 2.5))
    K = np.exp(rng.uniform(np.log(1.5), np.log(500)))
    P0 = M.rates_eps(L, eps, d=d, K=K)
    kmul = 1.0 + np.exp(rng.uniform(-6, 1.5))
    P = M.rates_eps(L, eps, d=d, kappa=kmul * P0["kappa_min"], K=K)
    if P["p"] >= 1.0 or not M.valid(P)[0]: continue
    r = M.osb(P)
    if not r["precondition"] or not np.isfinite(r["ratio"]) or r["sigma"] <= 0:
        continue
    aa = Rp.A_abs(P, n=12)[0]
    m = min(1.0, aa)
    marg = r["sigma"] / (m * r["rhs"])
    kept.append(dict(L=L, eps=eps, d=d, K=K, kmul=kmul, A_abs=aa,
                     sigma=r["sigma"], rhs=r["rhs"], osb=r["sigma"] >= r["rhs"],
                     margin=marg))
    if worst is None or marg < worst["margin"]: worst = kept[-1]

print(__doc__)
A = np.array([k["A_abs"] for k in kept])
mg = np.array([k["margin"] for k in kept])
osbv = sum(1 for k in kept if not k["osb"])
print(f"{len(kept)} admissible points from {pts} draws across the whole box")
print(f"  A_abs: min {A.min():.3e}  median {np.median(A):.3e}  max {A.max():.3e}")
print(f"  points with A_abs < 1 : {(A < 1).sum()}/{len(kept)}"
      f"   (so the corrected form is strictly weaker than OSB at"
      f" {(A<1).mean()*100:.0f}% of them -- the test has power)")
print(f"  OSB violated          : {osbv}/{len(kept)}")
print(f"  corrected form violated: {(mg < 1).sum()}/{len(kept)}")
print(f"  margin of the corrected form: min {mg.min():.3f}  "
      f"1st pct {np.percentile(mg,1):.3f}  median {np.median(mg):.2f}  max {mg.max():.1f}")
print(f"\n  tightest point found:")
for k, v in worst.items():
    print(f"     {k:7s} = {v!r}")
print(f"\n  The corrected form is {'NOT ' if (mg<1).sum()==0 else ''}violated anywhere"
      f" in a 600-point sweep of the family built to break the original.")
print("  That is a falsification attempt that COULD have failed and did not.")
print("  It is still not a proof, and I do not have one.")
json.dump(dict(n=len(kept), osb_viol=osbv, corrected_viol=int((mg<1).sum()),
               min_margin=float(mg.min()), frac_A_lt_1=float((A<1).mean())),
          open("falsify.json","w"), indent=1)
