"""Shiraishi's counterexample to the spectral dissipation-coherence trade-off.

Source: arXiv:2609.34352v1 (LaTeX e-print, OSB-counterexample-ver0.tex,
timestamped 2026-09-27 22:17 UTC). Conjecture under test: Oberreiter, Seifert
& Barato, Phys. Rev. E 106, 014106 (2022) --

    sigma_ss  >=  lambda_I^2 / lambda_R        when lambda_I >= lambda_R

where lambda_2 = -lambda_R + i*lambda_I is the eigenvalue of the transition
rate matrix with the SECOND LARGEST real part.

CONVENTION, stated because the paper's Fourier sign differs from mine and it
does not matter: R is a column generator, dp/dt = R p, sum_i R_ij = 0, and
R_ij is the rate of the jump j -> i.
"""
import numpy as np

K_DEFAULT = 10.0


def rates(L, q=0.75, d=3.0, kappa=None, K=K_DEFAULT):
    """Return the scalar parameters of the model at a given size.

    eps = L**-q with 1/2 < q < 1.   p = kappa * eps**2.

    The paper writes the non-negativity condition as `p > (d-1) eps^2`.  Both
    binding off-diagonals -- u >= v and b <= a/2 -- reduce to the SAME
    inequality a*sqrt(p) >= (d-1)*eps*(1-p), whose solution to leading order is
    p >= (d-1)^2 eps^2 / a^2, NOT (d-1) eps^2.  `kappa_min` below is the exact
    root; see `anchors.py::anchor_7`.
    """
    eps = L ** -float(q)
    a = 1.0 - eps
    # exact root of a^2 p = (d-1)^2 eps^2 (1-p)^2  in p, taking the branch p<1
    # a sqrt(p) = (d-1) eps (1-p)  ->  (d-1) eps p + a sqrt(p) - (d-1) eps = 0
    c = (d - 1.0) * eps
    s = (-a + np.sqrt(a * a + 4.0 * c * c)) / (2.0 * c)   # = sqrt(p) at equality
    kappa_min = (s * s) / eps ** 2
    if kappa is None:
        kappa = 1.05 * kappa_min          # just inside the allowed region
    p = kappa * eps ** 2
    u = a * p / (1.0 - p)
    h = (a + d * eps) / 2.0
    v = (d - 1.0) * eps * np.sqrt(p)
    b = (d - 1.0) * (1.0 - p) * eps / (2.0 * np.sqrt(p))
    theta = 2.0 * np.pi / L
    rA_p = eps / (2.0 * (1.0 - np.cos(theta))) + 1.0 / np.sin(theta)
    rA_m = eps / (2.0 * (1.0 - np.cos(theta))) - 1.0 / np.sin(theta)
    rBC = K / (2.0 * (1.0 - np.cos(theta)))
    return dict(L=L, q=q, d=d, K=K, eps=eps, p=p, kappa=kappa,
                kappa_min=kappa_min, a=a, u=u, h=h, v=v, b=b, theta=theta,
                rA_p=rA_p, rA_m=rA_m, rBC=rBC)


def Q_matrix(P):
    """Inter-layer 3x3 generator Q = Q0 + Q1 (columns sum to zero)."""
    a, u, h, v, b = P["a"], P["u"], P["h"], P["v"], P["b"]
    Q0 = np.array([[-a,      u,       u      ],
                   [ a / 2, -u - h,   h      ],
                   [ a / 2,  h,      -u - h  ]])
    Q1 = np.array([[ 0.0, -v,    v  ],
                   [ b,    0.0, -v  ],
                   [-b,    v,    0.0]])
    return Q0 + Q1


def R_full(P):
    """The full 3L x 3L generator.  State index = 3*i + layer, layer A=0,B=1,C=2."""
    L = P["L"]
    n = 3 * L
    R = np.zeros((n, n))
    for lay, (rp, rm) in enumerate(ring_rates(P)):
        for i in range(L):
            j = 3 * i + lay
            R[3 * ((i + 1) % L) + lay, j] += rp
            R[3 * ((i - 1) % L) + lay, j] += rm
            R[j, j] -= rp + rm
    Q = Q_matrix(P)
    for i in range(L):
        for A in range(3):
            for B in range(3):
                R[3 * i + A, 3 * i + B] += Q[A, B]
    return R


def pss(P):
    """Stationary distribution, exact by shift invariance."""
    L, p = P["L"], P["p"]
    cell = np.array([p, (1 - p) / 2, (1 - p) / 2]) / L
    return np.tile(cell, L)


def ring_rates(P):
    """[(r+, r-)] per layer, A B C.  One place, so an override is seen everywhere."""
    return [(P["rA_p"], P["rA_m"]), (P["rBC"], P["rBC"]), (P["rBC"], P["rBC"])]


def R_block(P, k):
    """Fourier block R-hat^k (3x3), built from the ACTUAL ring rates.

    The ring part of mode k contributes, per layer,
        -(r+ + r-)(1 - cos phi)  -  i (r+ - r-) sin phi,    phi = 2 pi k / L.
    With the paper's parametrisation (r+ + r-) = eps/(1-cos theta) in layer A
    and K/(1-cos theta) in B,C, this is exactly the paper's
        -F_k diag(eps, K, K) + 2i G_k diag(1,0,0)
    up to the sign of the imaginary part, which is the k <-> L-k Fourier
    convention and leaves the spectrum of R unchanged.
    """
    phi = 2.0 * np.pi * k / P["L"]
    diag = []
    for rp, rm in ring_rates(P):
        diag.append(-(rp + rm) * (1.0 - np.cos(phi)) - 1j * (rp - rm) * np.sin(phi))
    return Q_matrix(P) + np.diag(diag)


def sigma_ss(P, R=None, p=None):
    """Stationary entropy production rate, sum_ij R_ij p_j ln(R_ij p_j / R_ji p_i)."""
    if R is None:
        R = R_full(P)
    if p is None:
        p = pss(P)
    n = len(p)
    s = 0.0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            f = R[i, j] * p[j]
            g = R[j, i] * p[i]
            if f > 0.0:
                if g <= 0.0:
                    return np.inf          # irreversible transition
                s += f * np.log(f / g)
    return s


def lambda2(P, R=None):
    """(lambda_R, lambda_I) of the second-largest eigenvalue, from the full matrix."""
    if R is None:
        R = R_full(P)
    w = np.linalg.eigvals(R)
    w = w[np.argsort(-w.real)]
    # w[0] is the stationary eigenvalue 0; the second largest real part follows
    l2 = w[1]
    return -l2.real, abs(l2.imag), w


# --- fast, shift-invariant routes (validated against the dense ones) ----------

def sigma_ss_fast(P):
    """sigma_ss by edges, using shift invariance: L identical sites.

    Per site: 3 ring edge-pairs (one per layer, forward only -- the reverse term
    is picked up by the log) and 3 layer edge-pairs.  Rates and p^ss are both
    site-independent, so one site x L is exact.
    """
    L, p = P["L"], P["p"]
    cell = np.array([p, (1 - p) / 2, (1 - p) / 2]) / L
    tot = 0.0
    # ring: the pair (i -> i+1) at rate r+, (i+1 -> i) at rate r-; same p^ss
    for lay, (rp, rm) in enumerate(ring_rates(P)):
        pi = cell[lay]
        if rp <= 0 or rm <= 0:
            if rp != rm:
                return np.inf
            continue
        # forward term + backward term of the SAME undirected edge
        tot += L * ((rp * pi - rm * pi) * np.log((rp * pi) / (rm * pi)))
    Q = Q_matrix(P)
    for A in range(3):
        for B in range(3):
            if A == B:
                continue
            f, g = Q[A, B] * cell[B], Q[B, A] * cell[A]
            if f > 0:
                if g <= 0:
                    return np.inf
                tot += L * f * np.log(f / g)
    return tot


def spectrum_blocks(P):
    """Every eigenvalue of R, as the union over Fourier blocks k=0..L-1."""
    out = []
    for k in range(P["L"]):
        out.extend(np.linalg.eigvals(R_block(P, k)))
    return np.array(out)


def lambda2_blocks(P):
    w = spectrum_blocks(P)
    w = w[np.argsort(-w.real)]
    return -w[1].real, abs(w[1].imag), w


def rates_nonnegative(P):
    """All off-diagonal rates >= 0?  Returns (ok, worst_value, which)."""
    Q = Q_matrix(P)
    checks = [("rA_minus", P["rA_m"]), ("rA_plus", P["rA_p"]), ("rBC", P["rBC"])]
    for A in range(3):
        for B in range(3):
            if A != B:
                checks.append((f"Q[{A},{B}]", Q[A, B]))
    worst = min(checks, key=lambda t: t[1])
    return worst[1] >= 0.0, worst[1], worst[0]


def L_min_rate(q, Lcap=10 ** 7):
    """Smallest integer L >= 3 with r^-_A > 0.

    EXACT condition: eps/(2(1-cos t)) > 1/sin t with t = 2 pi / L, eps = L^-q,
    which rearranges to   L^-q > 2 tan(pi / L).
    The left side falls like L^-q and the right like 2 pi / L, so for q < 1 the
    inequality holds for all large L and the set of solutions is an up-set in L;
    bisection is therefore exact.  (For q >= 1 it never holds: returns None.)
    The asymptotic threshold is L ~ (2 pi)^(1/(1-q)), which is exact to +1.
    """
    def ok(L):
        return (float(L) ** -q) > 2.0 * np.tan(np.pi / L)
    if q >= 1.0 or not ok(Lcap):
        return None      # no solution at or below Lcap; asymptotically ~(2pi)^(1/(1-q))
    lo, hi = 3, Lcap
    while lo < hi:
        mid = (lo + hi) // 2
        if mid >= 3 and ok(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def match_max(a, b):
    """Max distance under a greedy nearest-neighbour matching of two multisets.
    Convention-free: no sort key, so conjugate-pair ordering cannot fool it."""
    b = list(b); worst = 0.0
    for z in a:
        j = min(range(len(b)), key=lambda i: abs(b[i] - z))
        worst = max(worst, abs(b[j] - z)); b.pop(j)
    return worst


def valid(P):
    """Is this parameter point inside the model's definition?  (ok, reason)"""
    ok, worst, which = rates_nonnegative(P)
    if not ok:
        return False, f"negative rate {which} = {worst:.4e}"
    return True, ""


def spectrum_fast(P):
    """All 3L eigenvalues, batched (L stacked 3x3 eigs)."""
    L = P["L"]
    phi = 2.0 * np.pi * np.arange(L) / L
    Q = Q_matrix(P).astype(complex)
    blocks = np.broadcast_to(Q, (L, 3, 3)).copy()
    for lay, (rp, rm) in enumerate(ring_rates(P)):
        blocks[:, lay, lay] += (-(rp + rm) * (1.0 - np.cos(phi))
                                - 1j * (rp - rm) * np.sin(phi))
    return np.linalg.eigvals(blocks).ravel()


def osb(P):
    """Everything the conjecture is about, at one parameter point.

    ratio = (lambda_I^2 / lambda_R) / sigma_ss;  ratio > 1 means VIOLATED.
    """
    ok, why = valid(P)
    w = spectrum_fast(P)
    w = w[np.argsort(-w.real)]
    zero, l2 = w[0], w[1]
    lR, lI = -l2.real, abs(l2.imag)
    sig = sigma_ss_fast(P)
    rhs = lI ** 2 / lR if lR > 0 else np.inf
    return dict(valid=ok, why=why, zero=abs(zero), lam2=l2, lR=lR, lI=lI,
                sigma=sig, rhs=rhs,
                ratio=(rhs / sig if sig > 0 else np.inf),
                precondition=(lI >= lR))


def rates_eps(L, eps, d=3.0, kappa=None, K=K_DEFAULT):
    """Same model, parametrised by eps DIRECTLY rather than by eps = L^-q.

    The paper's eps = L^-q with 1/2 < q < 1 is a convenient family for taking
    L -> infinity; at a FIXED L nothing in the construction refers to q, so the
    free parameters are (eps, d, kappa, K).  Searching over eps is therefore the
    honest version of "how small can the counterexample be".
    """
    eps = float(eps)
    a = 1.0 - eps
    c = (d - 1.0) * eps
    s = (-a + np.sqrt(a * a + 4.0 * c * c)) / (2.0 * c)
    kappa_min = (s * s) / eps ** 2
    if kappa is None:
        kappa = 1.05 * kappa_min
    p = kappa * eps ** 2
    u = a * p / (1.0 - p)
    h = (a + d * eps) / 2.0
    v = (d - 1.0) * eps * np.sqrt(p)
    b = (d - 1.0) * (1.0 - p) * eps / (2.0 * np.sqrt(p))
    theta = 2.0 * np.pi / L
    rA_p = eps / (2.0 * (1.0 - np.cos(theta))) + 1.0 / np.sin(theta)
    rA_m = eps / (2.0 * (1.0 - np.cos(theta))) - 1.0 / np.sin(theta)
    rBC = K / (2.0 * (1.0 - np.cos(theta)))
    q = -np.log(eps) / np.log(L) if L > 1 else np.nan
    return dict(L=L, q=q, d=d, K=K, eps=eps, p=p, kappa=kappa,
                kappa_min=kappa_min, a=a, u=u, h=h, v=v, b=b, theta=theta,
                rA_p=rA_p, rA_m=rA_m, rBC=rBC)
