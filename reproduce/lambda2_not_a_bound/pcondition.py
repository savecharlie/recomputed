"""The paper's stated parameter condition, checked.

arXiv:2609.34352v1 writes:

  "in order to satisfy the nonnegativity of the off-diagonal elements, which is
   equivalent to a*sqrt(p) >= (d-1)*eps*(1-p), we require p > (d-1)*eps^2."

The first clause is right.  I re-derive it: the only off-diagonals of Q=Q0+Q1
that can go negative are

      Q[0,1] = u - v        with u = a p/(1-p),  v = (d-1) eps sqrt(p)
      Q[2,0] = a/2 - b      with b = (d-1)(1-p) eps / (2 sqrt(p))
      Q[1,2] = h - v        with h = (a + d eps)/2

and u >= v and b <= a/2 reduce to the SAME inequality a sqrt(p) >= (d-1) eps (1-p),
which is a pleasing feature of the construction.  (h >= v is slack: h = O(1),
v = O(eps^2).)

The second clause does NOT follow.  Squaring gives p >= (d-1)^2 eps^2 (1-p)^2/a^2,
so to leading order the requirement is  p > (d-1)^2 eps^2  --  the square.
"""
import sys, numpy as np
sys.path.insert(0, ".")
import model as M

print(__doc__)
print("Testing p = c * eps^2 at d=3, where the paper's clause allows c > d-1 = 2")
print("but my derivation requires c > (d-1)^2 = 4.\n")
print(f"{'L':>7} {'q':>5} {'c':>6} {'p':>12} {'u-v':>13} {'a/2-b':>13}  verdict")
bad_in_paper_window = 0
for L, q in [(2000, 0.75), (10000, 0.75), (120, 0.55), (100000, 0.8)]:
    for c in (2.01, 3.0, 3.9, 4.0, 4.01, 4.5):
        P = M.rates(L, q=q, d=3.0, kappa=c)
        Q = M.Q_matrix(P)
        uv, ab = Q[0, 1], Q[2, 0]
        neg = min(uv, ab) < 0
        inpaper = c > 2.0            # allowed by the paper's stated clause
        if neg and inpaper:
            bad_in_paper_window += 1
        print(f"{L:7d} {q:5.2f} {c:6.2f} {P['p']:12.4e} {uv:13.5e} {ab:13.5e}"
              f"  {'NEGATIVE RATE' if neg else 'ok'}")
    print()

print(f"settings allowed by the paper's clause that have a negative rate: "
      f"{bad_in_paper_window}")
print()
print("The exact threshold in c, solving a^2 c eps^2 = (d-1)^2 eps^2 (1 - c eps^2)^2:")
for L, q in [(120, 0.55), (2000, 0.75), (100000, 0.8)]:
    P = M.rates(L, q=q, d=3.0)
    print(f"  L={L:6d} q={q}  eps={P['eps']:.5e}   c_min = {P['kappa_min']:.6f}"
          f"   (leading order (d-1)^2/a^2 = {(3-1)**2/P['a']**2:.6f})")
print()
print("So: the clause should read p > (d-1)^2 eps^2.  Everything downstream is")
print("unaffected -- p = O(eps^2) is what the entropy-production argument uses,")
print("and (d-1)^2 eps^2 is still O(eps^2).  It is a misprint in a prerequisite,")
print("not a hole in the result: the counterexample stands.")
