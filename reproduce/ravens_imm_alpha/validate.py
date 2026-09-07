"""Validate the IMM rebuild against the three answers the paper already knows.
Nothing this instrument says about alpha is worth reading until these pass."""
import numpy as np, imm

T, R = 4000, 300
print("known answer 1+2: <n(t)> = [(1+b) q t]^(1/(1+b)) and sigma = sqrt(<n>/(1+2b))")
print("   beta    q     t     <n> sim    <n> theory   ratio      sd sim   sd theory")
for beta in (0.37, 0.72, 1.0):
    for q in (0.6, 1.0):
        N, _ = imm.ensemble(R, T, alpha=1.0, beta=beta, q=q, seed=11)
        for t in (500, 4000):
            sim = N[:, t].mean(); th = ((1 + beta) * q * t) ** (1 / (1 + beta))
            sd = N[:, t].std(); sdth = np.sqrt(sim / (1 + 2 * beta))
            print(f"  {beta:5.2f} {q:5.2f} {t:5d}   {sim:8.2f}   {th:9.2f}   "
                  f"{sim/th:6.4f}   {sd:7.3f}   {sdth:8.3f}")

print("\nknown answer 3: at alpha=1, <m_i(t)> ~ t  (log-log slope -> 1)")
N, M = imm.ensemble(R, T, alpha=1.0, beta=0.72, q=0.6, seed=12)
t = np.arange(T + 1)
lo, hi = 400, T
for i in range(3):
    y = M[:, :, i].mean(0)
    s = np.polyfit(np.log(t[lo:hi]), np.log(y[lo:hi]), 1)[0]
    print(f"  rank {i+1}: slope {s:.4f}   m_i({T}) = {y[T]:.1f}")
