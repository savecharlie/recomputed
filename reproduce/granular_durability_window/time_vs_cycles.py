"""Is the lifetime set by CYCLES or by elapsed TIME?

The drive is a sinusoid; the paper reports amplitude 0.31 mm and quotes 40 and
50 Hz.  With the amplitude held fixed, G = A(2*pi*f)^2/g, so f = sqrt(G g / A)/2pi
and the five G rows of Fig. 3 are five FREQUENCIES.  A cycle-counted lifetime and
a clock-counted lifetime then differ by a factor f ~ sqrt(G) -- which is exactly
the direction that matters at the low-G end, where Eq. 2 fails hardest.
"""
import json, numpy as np
P = json.load(open("fig3_points.json"))
canon = np.array([0.41, 0.89, 1.40, 2.40, 3.07])
al = np.array([p["alpha"] for p in P])
G  = canon[np.argmin(np.abs(np.array([p["G"] for p in P])[:, None] - canon), 1)]
N  = np.array([p["Nf"] for p in P])

A, g = 0.31e-3, 9.81
f = np.sqrt(G * g / A) / (2 * np.pi)
print("implied drive frequency for each G row (amplitude fixed at 0.31 mm):")
for gg in canon:
    print(f"   G = {gg:>4}  ->  f = {np.sqrt(gg*g/A)/(2*np.pi):5.1f} Hz")
print("  (the paper's own 50 Hz / 0.31 mm gives G = "
      f"{A*(2*np.pi*50)**2/g:.2f}, and 40 Hz gives {A*(2*np.pi*40)**2/g:.2f})\n")

T = N / f                                   # seconds to failure
def fit(y, tag):
    L = np.log(y); a2 = al**2; ig = 1/G; pr = a2*ig
    X = np.column_stack([np.ones(len(L)), pr])
    b,*_ = np.linalg.lstsq(X, L, rcond=None); r = L - X@b
    r2 = 1 - (r@r)/((L-L.mean())**2).sum()
    # and the additive model
    X2 = np.column_stack([np.ones(len(L)), a2, ig])
    b2,*_ = np.linalg.lstsq(X2, L, rcond=None); r2b = L - X2@b2
    r2_2 = 1 - (r2b@r2b)/((L-L.mean())**2).sum()
    print(f"{tag}:  product form  C = {b[1]:6.2f}  R2 = {r2:.3f}   |  "
          f"additive  B(a^2) = {b2[1]:6.2f}  C(1/G) = {b2[2]:6.2f}  R2 = {r2_2:.3f}")
fit(N, "lifetime in CYCLES N_f     ")
fit(T, "lifetime in SECONDS N_f / f")

print("\nper-alpha sensitivity to 1/G, which Eq. 2 says must equal C1*alpha^2:")
print(f"{'alpha':>5} {'slope(lnN vs 1/G)':>18} {'slope(lnT vs 1/G)':>18} {'26.8*a^2':>10}")
for a in sorted(set(al)):
    m = al == a
    if len(set(G[m])) < 2: continue
    kn = np.polyfit(1/G[m], np.log(N[m]), 1)[0]
    kt = np.polyfit(1/G[m], np.log(T[m]), 1)[0]
    print(f"{a:>5} {kn:>18.2f} {kt:>18.2f} {26.8*a*a:>10.2f}")

print("\n--- does Eq. 2 with the paper's constants predict its own alpha>=0.9 null result? ---")
p0, C1 = 0.0037, 26.8
for a in (0.9, 1.0):
    for gg in (3.07, 2.40):
        Nf = (1/p0)*np.exp(C1*a*a/gg)
        ff = np.sqrt(gg*g/A)/(2*np.pi)
        print(f"  alpha={a}, G={gg}:  Eq.2 N_f = {Nf:10.3g} cycles = "
              f"{Nf/ff/3600:8.2f} h at {ff:.0f} Hz   "
              f"{'-> should have FAILED inside the 10 h window' if Nf/ff/3600 < 10 else '-> consistent with no failure'}")
