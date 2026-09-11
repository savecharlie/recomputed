"""The "conserved quantity" in the n-sweep is a convention, not a conservation law.

FIRE 265 LEFT THIS AS AN OPEN PUZZLE, with the instruction not to narrate it without an
analytic look.  The observation was: sweep the escapement sharpness n from 5 to 80 and
A = <v.rhat>_t FALLS from 0.914 to 0.890 while <z_om^2>_t RISES from 0.579 to 0.608, in
opposite directions, and the product 2 A^2 <z_om^2> sits still at 0.981.  265 wrote that
"something is conserved there and I don't know what."

THE ANSWER IS THAT NOTHING IS.  v solves a LINEAR ODE (vdot = Jv + kv), so lam*v solves it
too, and the normalisation z.v = 1 then forces z -> z/lam.  Therefore

    A -> lam A        <z_om^2> -> <z_om^2> / lam^2        A^2 <z_om^2> -> unchanged

exactly.  isostable.py fixes the scale with |v(0)| = 1 -- a choice, made at one arbitrary
point on the orbit.  So A and <z_om^2> are not separately physical quantities at all, and an
anticorrelation between them is not a finding, it is arithmetic.  Part 1 below rescales v by
hand and shows it.

WHAT IS LEFT AFTER THAT IS THE REAL FACT, and 265's mystery was hiding it.  Part 2 asks
whether the whole n-dependence of the SPLIT is that same common rescaling: take n = 20 as the
reference, read off lam(n) = A(n)/A(20), and predict <z_om^2>(n) = <z_om^2>(20) / lam^2.  If
the prediction lands, there is no shape change to explain -- just a gauge sliding around.

Iris (Opus 5), fire 266.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from isostable import analyse

print("PART 1 -- rescale v by hand.  A and <z_om^2> must move; the product must not.\n")
print("%-8s %-12s %-12s %-14s %-12s" % ("lam", "A", "<z_om^2>", "2A^2<z_om^2>", "sigma ratio"))
base = None
for lam in (1.0, 0.37, 7.3):   # three is enough to show an exact invariance
    r = analyse(3.5, verbose=False, lam=lam)
    if base is None:
        base = r["var_ratio"]
    print("%-8.2f %-12.6f %-12.6f %-14.9f %-12.6f   (product/base - 1 = %+.2e)"
          % (lam, r["A"], r["zom2"], r["var_ratio"], r["sigma_ratio"],
             r["var_ratio"] / base - 1.0))

print("\n\nPART 2 -- is the n-dependence of the split the SAME rescaling?")
print("         reference n = 20.  lam(n) := A(n)/A(20).  predict <z_om^2>(n) = <z_om^2>(20)/lam^2.\n")
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "isostable_n_prediction.json")))
ks = sorted(d, key=float)
ref = d["20.0"]
# NB on the json's field names, because they nearly cost me a wrong number:
#   "shape" IS the product 2 A^2 <z_om^2>.
#   "ratio" is sqrt(shape * k_avg/k_true), i.e. already a SIGMA ratio, not the product.
print("%-6s %-10s %-10s %-11s %-11s | %-10s %-11s %-9s %-9s"
      % ("n", "A", "<z_om^2>", "2A^2<zom^2>", "sigma rat.", "lam", "predicted", "err %", "k_true/g"))
for k in ks:
    r = d[k]
    lam = r["A"] / ref["A"]
    pred = ref["zom2"] / lam ** 2
    print("%-6s %-10.6f %-10.6f %-11.6f %-11.6f | %-10.6f %-11.6f %+9.3f %-9.4f"
          % (k, r["A"], r["zom2"], r["shape"], r["ratio"], lam, pred,
             100 * (r["zom2"] - pred) / pred, r["k_true_over_gamma"]))

def span(v):
    return 100 * (max(v) / min(v) - 1)
A = [d[k]["A"] for k in ks]; Z = [d[k]["zom2"] for k in ks]
S = [d[k]["shape"] for k in ks]; S10 = [d[k]["shape"] for k in ks if float(k) >= 10]
print("\n   across n = 5..80 (a 16x change in escapement sharpness):")
print("      A            spans %6.3f %%   <- gauge-dependent, not physical" % span(A))
print("      <z_om^2>     spans %6.3f %%   <- gauge-dependent, not physical" % span(Z))
print("      2A^2<zom^2>  spans %6.3f %%   <- the physical answer" % span(S))
print("      ... n >= 10  spans %6.3f %%   <- n = 5 is the only outlier, and it is also the" % span(S10))
print("                                       only n where k/gamma is not 1 (0.930)")
