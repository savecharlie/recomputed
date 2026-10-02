"""Anchors. Nothing new is reported until every one of these passes.

Each anchor reproduces something whose value is published, exactly derivable,
or known by construction -- BEFORE the instrument is allowed to say anything
about the conjecture.
"""
import sys, numpy as np
sys.path.insert(0, "/home/ivy/iris-the-maker/reading/lambda2-not-a-bound")
import model as M

FAIL = []
def chk(name, cond, detail=""):
    print(("  PASS " if cond else "  FAIL ") + name + ("   " + detail if detail else ""))
    if not cond:
        FAIL.append(name)

print("ANCHOR 1  generator structure (columns sum to zero, off-diagonals >= 0)")
for L, q in [(2000, 0.75), (4000, 0.6)]:
    P = M.rates(L, q=q)
    R = M.R_full(P)
    chk(f"L={L} columns sum to 0", np.max(np.abs(R.sum(axis=0))) < 1e-9,
        f"max|col sum| = {np.max(np.abs(R.sum(axis=0))):.2e}")
    off = R - np.diag(np.diag(R))
    chk(f"L={L} off-diagonals >= 0", off.min() >= 0, f"min = {off.min():.6e}")

print("\nANCHOR 2  p^ss of Q is exactly (p, (1-p)/2, (1-p)/2)   [paper, eq. after Q]")
for L in (200, 2000):
    P = M.rates(L)
    Q = M.Q_matrix(P)
    pq = np.array([P["p"], (1 - P["p"]) / 2, (1 - P["p"]) / 2])
    chk(f"L={L} Q p = 0", np.max(np.abs(Q @ pq)) < 1e-12,
        f"max|Q p| = {np.max(np.abs(Q @ pq)):.2e}")
    # and separately for Q0 and Q1, which the paper asserts share it
    a, u, h, v, b = P["a"], P["u"], P["h"], P["v"], P["b"]
    Q0 = np.array([[-a, u, u], [a/2, -u-h, h], [a/2, h, -u-h]])
    Q1 = np.array([[0, -v, v], [b, 0, -v], [-b, v, 0]])
    chk(f"L={L} Q0 p = 0 and Q1 p = 0",
        max(np.max(np.abs(Q0 @ pq)), np.max(np.abs(Q1 @ pq))) < 1e-12)

print("\nANCHOR 3  characteristic polynomial of Q matches the paper's closed form")
print("          det(zI-Q) = z[(z + (1-eps)/(1-p) + (d/2)eps)^2 + ((d-1)^2 - d^2/4)eps^2]")
for L in (200, 2000, 20000):
    P = M.rates(L)
    eps, p, d = P["eps"], P["p"], P["d"]
    lam = M.np.linalg.eigvals(M.Q_matrix(P))
    pred0 = 0.0
    predpm = -(1 - eps) / (1 - p) - d / 2 * eps + 1j * np.sqrt((d-1)**2 - d**2/4) * eps
    got = np.array(sorted(lam, key=lambda z: -z.real))
    err = max(abs(got[0] - pred0),
              min(abs(got[1] - predpm), abs(got[1] - predpm.conjugate())))
    chk(f"L={L} eigenvalues of Q match eq.(lmd-R0)", err < 1e-12, f"max err = {err:.2e}")

print("\nANCHOR 4  Fourier block-diagonalisation: spectrum(R) == union_k spectrum(R-hat^k)")
for L, q in [(8, 0.6), (60, 0.6), (120, 0.55), (300, 0.75)]:
    P = M.rates(L, q=q)
    e = M.match_max(np.linalg.eigvals(M.R_full(P)), M.spectrum_blocks(P))
    chk(f"L={L} q={q} spectra agree (matched, no sort key)", e < 1e-7,
        f"max matched |diff| = {e:.2e}")

print("\nANCHOR 5  p^ss of the FULL R is the shift-invariant tile (independent route)")
for L, q in [(120, 0.55), (300, 0.75)]:
    P = M.rates(L, q=q)
    R = M.R_full(P)
    w, V = np.linalg.eig(R)
    v = np.real(V[:, np.argmax(w.real)]); v = v / v.sum()
    chk(f"L={L} numerical null vector == analytic tile",
        np.max(np.abs(v - M.pss(P))) < 1e-9,
        f"max|diff| = {np.max(np.abs(v - M.pss(P))):.2e}")

print("\nANCHOR 6  sigma_ss: dense O(n^2) route == fast shift-invariant route")
# ONLY at parameter points inside the model's definition: outside it r^-_A < 0
# and both routes correctly return inf, which is not a test of agreement.
for L, q in [(70, 0.55), (120, 0.55), (300, 0.6), (2000, 0.75)]:
    P = M.rates(L, q=q)
    ok, why = M.valid(P)
    assert ok, f"L={L} q={q} is outside the model: {why}"
    if L <= 300:
        s1 = M.sigma_ss(P)
    else:
        s1 = None
    s2 = M.sigma_ss_fast(P)
    if s1 is None:
        print(f"  (L={L} dense route skipped, n^2 too slow)  fast = {s2:.6e}")
    else:
        rel = abs(s1 - s2) / max(abs(s1), 1e-300)
        chk(f"L={L} two routes agree", rel < 1e-10, f"{s1:.8e} vs {s2:.8e}  rel {rel:.1e}")

print("\nANCHOR 7  KNOWN ZERO: with no driving anywhere, sigma_ss == 0 and lambda is real")
P = M.rates(2000)
P0 = dict(P); P0["v"] = 0.0; P0["b"] = 0.0          # kill inter-layer driving
P0["rA_p"] = P0["rA_m"] = P0["rBC"]                 # kill ring driving
# R_block now reads ring_rates(P), so the override is visible to the blocks too
s = M.sigma_ss_fast(P0)
chk("undriven sigma_ss == 0", abs(s) < 1e-12, f"sigma = {s:.3e}")
w = M.spectrum_blocks(P0)
chk("undriven spectrum is real", np.max(np.abs(w.imag)) < 1e-9,
    f"max|Im| = {np.max(np.abs(w.imag)):.2e}")

print("\nANCHOR 7b the rate floor: EXACT condition r^-_A > 0  <=>  L^-q > 2 tan(pi/L)")
print("          (the paper says only that q<1 'confirms the nonnegativity of r_-';")
print("           that is asymptotic in L -- at fixed q there is a hard floor on L.)")
bad = []
for q in (0.51, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80):
    got = M.L_min_rate(q)
    # the exact inequality, evaluated at the claimed threshold and one below it
    f_ok = (got ** -q) > 2 * np.tan(np.pi / got)
    f_no = ((got - 1) ** -q) > 2 * np.tan(np.pi / (got - 1))
    if not (f_ok and not f_no):
        bad.append((q, got))
chk("exact inequality is the threshold at 7 values of q", not bad, str(bad))
# and the asymptotic form, reported WITH its error rather than asserted exact
devs = []
for q in (0.51, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80):
    devs.append(M.L_min_rate(q) - int(np.ceil((2 * np.pi) ** (1.0 / (1.0 - q)))))
chk("asymptotic ceil((2pi)^(1/(1-q))) is within +1 of exact, never low",
    all(0 <= dv <= 1 for dv in devs), f"deviations {devs} (exact at {devs.count(0)}/7)")

print("\nANCHOR 8  POSITIVE CONTROL: the OSB bound HOLDS on a plain driven ring")
print("          (OSB 2022's own canonical system -- if my detector flags this, it is broken)")
def ring(L, rp, rm):
    R = np.zeros((L, L))
    for i in range(L):
        R[(i+1) % L, i] += rp; R[(i-1) % L, i] += rm; R[i, i] -= rp + rm
    return R
ok = True
for L in (4, 6, 10, 20, 50):
    for rp, rm in [(1.0, 0.5), (1.0, 0.1), (1.0, 0.01), (2.0, 1.9), (1.0, 0.9)]:
        R = ring(L, rp, rm)
        p = np.ones(L) / L
        sig = L * (rp - rm) * p[0] * np.log(rp / rm)
        w = np.linalg.eigvals(R); w = w[np.argsort(-w.real)]
        lR, lI = -w[1].real, abs(w[1].imag)
        if lI >= lR and lI > 0:
            rhs = lI ** 2 / lR
            if sig < rhs - 1e-12:
                ok = False
                print(f"    VIOLATION at L={L} rp={rp} rm={rm}: {sig:.4f} < {rhs:.4f}")
chk("driven ring never violates OSB (25 settings)", ok)

print("\n" + ("ALL ANCHORS PASS" if not FAIL else f"FAILURES: {FAIL}"))
sys.exit(1 if FAIL else 0)
