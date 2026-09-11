"""A fourth ruler for sigma_R, built on a different principle: sample one phase per cycle.

WHY A FOURTH.  Three rulers have now disagreed at gamma = 0.20.  The boxcar said 0.9856,
the Fourier-fit said 0.9565 with 4 harmonics and 0.9317 with 6, and refitting the same data
with mixed terms said 1.0563.  ruler_domain.py showed on synthetic data with a sigma I typed
in that none of that is a small-N_c artifact -- the Fourier ruler is unbiased to 0.5% at
N_c = 0.83 with a sinusoidal wobble of ANY amplitude up to 6 sigma.  What it cannot handle is
harmonic content outside its basis, and the escapement is a sharp switch, so the clock has
plenty.

Every ruler so far shares one assumption: that the within-cycle wobble can be MODELLED and
subtracted.  Boxcar models it as "whatever averages out over one period"; the Fourier fit
models it as a truncated harmonic series.  Both are a model OF THE CONTAMINANT, which is the
thing Kotoku et al. (arXiv:2609.11268, 10 Sep 2026) set out to avoid: their CROP kills noise
by correlating against something it cannot correlate with, rather than by estimating it.

THE MOVE.  The wobble is not noise, it is a deterministic function of the cycle phase.  So
stop modelling it and stand still instead: record R only at a fixed phase of the cycle.  The
wobble is then the SAME VALUE at every sample and drops out of the variance identically, with
no basis, no truncation, and no harmonic content to get wrong.

Four markers, each of which makes R trivial to read:
    theta = 0 going up    -> R = |omega|
    theta = 0 going down  -> R = |omega|
    omega = 0, theta > 0  -> R = |theta|
    omega = 0, theta < 0  -> R = |theta|

THE ONE THING IT COSTS, stated up front.  If R(t) = A(t) * rho(Phi) with A slow, then the
strobe at phase Phi_0 measures sigma_A * rho(Phi_0), while the autocovariance rulers read the
coefficient of e^{-kL}, which is sigma_A * <rho>^2 under the square root -- i.e. sigma_A *
<rho>.  So each strobe needs dividing by rho(Phi_0)/<rho>, computed from the NOISELESS orbit.
That is exact arithmetic on a deterministic curve, not a fitted correction.

AND THE SELF-TEST THE OTHER RULERS CANNOT DO.  Separability is an assumption, so check it:
four phases, four corrections, four answers.  If they agree, R = A*rho(Phi) held on the real
data.  If they fan out, amplitude and phase are coupled and "the" sigma_R needs more care than
any of these rulers has been giving it.

ANCHOR.  gamma = 0.05 is where the boxcar, the Fourier fit and the isostable theory already
agree (0.9792 / 0.9806 / 0.9813).  A new ruler that does not reproduce a number I already know
is wrong, whatever it says at gamma = 0.20.

Iris (Opus 5), fire 266.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius
import amplitude_phase as ap

T_DEF, TR = 5e-4, 0.4
MARKERS = ("th_up", "th_dn", "om_up", "om_dn")


def orbit_shape(esc, gamma, eps, Rstar, dt=2e-4):
    """rho(Phi) on the noiseless limit cycle: <rho> in time, and rho at each of the 4 markers.

    Returns (mean_rho, {marker: rho}, period).  Deterministic; no noise anywhere in here.
    """
    th, om = Rstar, 0.0
    for _ in range(int(round(80 * 2 * math.pi / dt))):      # settle onto the cycle
        om += (-th - gamma * om + eps * esc.f(th, om)) * dt
        th += om * dt
    # one full period, starting from an upward theta crossing
    pth, pom = th, om
    while not (pth < 0 <= th and om > 0):
        pth, pom = th, om
        om += (-th - gamma * om + eps * esc.f(th, om)) * dt
        th += om * dt
    hit = {}
    hit["th_up"] = abs(om)
    Rs, t, n = [], 0.0, 0
    pth, pom = th, om
    while True:
        om += (-th - gamma * om + eps * esc.f(th, om)) * dt
        th += om * dt
        t += dt; n += 1
        Rs.append(math.hypot(th, om))
        if pth >= 0 > th and om < 0 and "th_dn" not in hit:
            hit["th_dn"] = abs(om)
        if pom >= 0 > om and th > 0 and "om_dn" not in hit:
            hit["om_dn"] = abs(th)
        if pom < 0 <= om and th < 0 and "om_up" not in hit:
            hit["om_up"] = abs(th)
        if n > 100 and pth < 0 <= th and om > 0:
            break
        pth, pom = th, om
    return float(np.mean(Rs)), hit, t


def strobe_run(u, gamma, seed, dt=2e-3, n_traj=1024, t_total=3200.0, T=T_DEF, n=20.0):
    """Integrate the SDE; record R at the four phase markers, sub-step interpolated."""
    esc = Grasshopper(TR, n)
    eps = math.pi * gamma * TR * u / 4.0
    Rstar = limit_cycle_radius(esc, gamma, eps)
    p = ap.predict(esc, gamma, eps, T, Rstar)
    rng = np.random.default_rng(seed)
    sig = math.sqrt(2.0 * gamma * T * dt)
    th = np.full(n_traj, Rstar); om = np.zeros(n_traj)
    for _ in range(int(round((30.0 / gamma) / dt))):
        om += (-th - gamma * om + eps * esc.f(th, om)) * dt + sig * rng.standard_normal(n_traj)
        th += om * dt
    buf = {m: [] for m in MARKERS}
    nst = int(round(t_total / dt))
    for _ in range(nst):
        pth, pom = th, om
        om = om + (-th - gamma * om + eps * esc.f(th, om)) * dt + sig * rng.standard_normal(n_traj)
        th = th + om * dt
        # theta crossings -> R = |omega| interpolated to theta = 0
        c = (pth < 0) & (th >= 0) & (om > 0)
        if c.any():
            f = -pth[c] / (th[c] - pth[c]); buf["th_up"].append(np.abs(pom[c] + f * (om[c] - pom[c])))
        c = (pth >= 0) & (th < 0) & (om < 0)
        if c.any():
            f = pth[c] / (pth[c] - th[c]); buf["th_dn"].append(np.abs(pom[c] + f * (om[c] - pom[c])))
        # omega crossings -> R = |theta| interpolated to omega = 0
        c = (pom >= 0) & (om < 0) & (th > 0)
        if c.any():
            f = pom[c] / (pom[c] - om[c]); buf["om_dn"].append(np.abs(pth[c] + f * (th[c] - pth[c])))
        c = (pom < 0) & (om >= 0) & (th < 0)
        if c.any():
            f = -pom[c] / (om[c] - pom[c]); buf["om_up"].append(np.abs(pth[c] + f * (th[c] - pth[c])))
    mean_rho, hit, per = orbit_shape(esc, gamma, eps, Rstar)
    out = dict(u=u, gamma=gamma, seed=seed, sigma_th=p['sigma_R'], period=per,
               mean_rho=mean_rho, rho=hit, n_samp={}, raw={}, corrected={}, ratio={})
    for m in MARKERS:
        v = np.concatenate(buf[m]) if buf[m] else np.array([])
        out["n_samp"][m] = int(v.size)
        s = float(v.std(ddof=1)) if v.size > 2 else float("nan")
        out["raw"][m] = s
        out["corrected"][m] = s * mean_rho / hit[m]
        out["ratio"][m] = out["corrected"][m] / p['sigma_R']
    return out


if __name__ == "__main__":
    print("STROBE -- sigma_R from one phase per cycle.  no wobble model anywhere.\n")
    print("ANCHOR gamma=0.05: boxcar 0.9792, Fourier 0.9806, isostable theory 0.9813")
    print("       gamma=0.20: boxcar 0.9856, Fourier 0.9565 (4h) / 0.9317 (6h), theory 0.9315\n")
    rows = []
    print("%-7s %-5s | %-8s %-8s %-8s %-8s | %-8s" % ("gamma", "seed", *MARKERS, "spread"))
    for gamma in (0.05, 0.10, 0.20):
        acc = {m: [] for m in MARKERS}
        for seed in (5, 6, 7):
            r = strobe_run(3.5, gamma, seed)
            rows.append(r)
            for m in MARKERS:
                acc[m].append(r["ratio"][m])
            vals = [r["ratio"][m] for m in MARKERS]
            print("%-7.2f %-5d | %-8.4f %-8.4f %-8.4f %-8.4f | %-8.4f"
                  % (gamma, seed, *vals, max(vals) - min(vals)), flush=True)
            json.dump(rows, open("strobe.json", "w"), indent=1)
        mv = [float(np.mean(acc[m])) for m in MARKERS]
        print("%-7.2f %-5s | %-8.4f %-8.4f %-8.4f %-8.4f | %-8.4f   <- mean"
              % (gamma, "mean", *mv, max(mv) - min(mv)))
        r0 = rows[-1]
        print("        rho/<rho> = " + "  ".join("%s %.4f" % (m, r0["rho"][m] / r0["mean_rho"])
                                                 for m in MARKERS)
              + "   (%d samples each)\n" % r0["n_samp"]["th_up"])
    json.dump(rows, open("strobe.json", "w"), indent=1)
