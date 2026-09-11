"""Is there a converged answer at gamma = 0.2, or does the estimator simply have none?

WHAT ruler_domain.py SETTLED, on synthetic data with a sigma I typed in myself.  My
prediction that the estimator biases as cycles-per-correlation-time N_c falls below 1 is
DEAD.  plain and h6 recover the true sigma to <0.5% at every N_c from 4.0 down to 0.6,
including N_c = 0.83, which is exactly where the real gamma = 0.20 run sits.  Only the
mixed basis (collinear extra columns) and L0 = 2 cycles (throws away the window the
exponential lives in) fall apart, and only at N_c = 0.6.

So the 30% spread on the REAL gamma = 0.2 data is not a small-N_c artifact.  Something is
in that signal which my synthetic did not have.

THE TELL.  plain (4 harmonics) and h6 (6 harmonics) are IDENTICAL to four decimals on
synthetic data at every N_c -- 1.0042/1.0042 at 0.60, 0.9988/0.9987 at 0.83.  On the real
gamma = 0.2 data they are 0.9565 and 0.9317, 2.6% apart.  Two bases that cannot be told
apart on a two-harmonic wobble are 2.6% apart on the clock, so the clock's within-cycle
wobble has real power at the 5th and 6th harmonic.  Which it should: the escapement is a
sharp switch (n = 20) that kicks once per cycle, and an impulsive kick is broadband.  The
kick scales with eps, which scales with gamma, so this is small at 0.05 and not at 0.20.

THE TEST, and it is a convergence argument rather than a taste argument.  Scan nharm.  If
the estimate plateaus as the basis is enriched, the plateau is the answer and the plateau
was chosen by convergence, not by which number I liked.  If it never settles, the
estimator has no answer at gamma = 0.2 and I say so.

Part 1 validates that reasoning on synthetic data with a SHARP kick and a known sigma.
Part 2 runs the same scan on the clock.

PREDICTION, written before running:
 (a) synthetic-with-sharp-kick reproduces the signature -- nharm = 4 biased relative to the
     converged value, drifting monotonically and settling by nharm ~ 10-14 onto the sigma I
     typed in;
 (b) the real gamma = 0.05 scan is flat from nharm = 2, because its kick is 4x weaker;
 (c) the real gamma = 0.20 scan drifts and then plateaus.  I do NOT predict where.  The
     isostable value is 0.9315 and nharm = 6 already gives 0.9317, and if the plateau lands
     there I will have to say plainly that a two-point coincidence is not evidence.

Iris (Opus 5), fire 266.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ruler_domain import autocov_chunked, P, OM, DTR, NENS, BLOCK
from ruler_domain import fit as fit_fixedOm      # synthetic only: its Omega IS OM
from fit_robust import sim, fit as fit_withOm    # real data: Omega must come from the orbit

NHARMS = (2, 3, 4, 6, 8, 10, 12, 16, 20, 26)


def sharp_block(tau, sigma, wobble, sharp, NT, seed, b, ntraj=BLOCK):
    """OU + a periodic kick shaped like the escapement's: a narrow bump once per cycle.

    bump(phi) = exp(sharp*(cos(phi)-1)) has bandwidth ~ sqrt(sharp) harmonics, so sharp=25
    puts real power out to j ~ 5, which is the regime the tell above points at.
    """
    rng = np.random.default_rng((seed, b))
    a = math.exp(-DTR / tau)
    s = sigma * math.sqrt(1 - a * a)
    x = np.empty((NT, ntraj))
    x[0] = sigma * rng.standard_normal(ntraj)
    for i in range(1, NT):
        x[i] = a * x[i - 1] + s * rng.standard_normal(ntraj)
    if wobble > 0:
        t = (np.arange(NT) * DTR)[:, None]
        ph = rng.uniform(0, 2 * math.pi, ntraj)[None, :]
        bump = np.exp(sharp * (np.cos(OM * t + ph) - 1.0))
        bump -= bump.mean(axis=0, keepdims=True)
        bump /= bump.std(axis=0, keepdims=True).mean()
        x += wobble * sigma * bump
    return x


def synth_scan(tau, wobble, sharp, seed, sigma=0.0158113939, t_total=3200.0):
    NT = int(round(t_total / DTR))
    C = autocov_chunked(lambda b: sharp_block(tau, sigma, wobble, sharp, NT, seed, b),
                        NENS // BLOCK, NT, int(round(4 * tau / DTR)))
    return [fit_fixedOm(C, tau, nharm=nh)[0] / sigma for nh in NHARMS]


if __name__ == "__main__":
    out = {"nharms": list(NHARMS), "synthetic": [], "clock": []}
    hdr = "%-26s | " % "" + " ".join("%-7d" % n for n in NHARMS)

    print("PART 1 -- synthetic, sigma TYPED IN, sharp kick, at the real gamma=0.2 geometry")
    print("          N_c = 0.83 (tau = 5.0).  every number below should go to 1.0000.\n")
    print(hdr)
    for wobble, sharp in ((1.5, 25.0), (3.0, 25.0), (3.0, 60.0)):
        for seed in (301, 302):
            r = synth_scan(0.83 * P, wobble, sharp, seed)
            print("%-26s | " % ("wob=%.1f sharp=%.0f s=%d" % (wobble, sharp, seed))
                  + " ".join("%-7.4f" % v for v in r), flush=True)
            out["synthetic"].append(dict(wobble=wobble, sharp=sharp, seed=seed, ratio=r))
            json.dump(out, open("harm_convergence.json", "w"), indent=1)
        print()

    print("\nPART 2 -- the clock.  isostable predicts 0.9813 at gamma=0.05, 0.9315 at 0.20\n")
    print(hdr)
    for gamma in (0.05, 0.20):
        for seed in (5, 6, 7):
            C, dtr, Om, s_th, tau = sim(3.5, gamma, seed)
            assert abs(dtr - DTR) < 1e-12, dtr
            # fit_robust.fit takes Omega explicitly.  The bug this replaces: the first
            # run used ruler_domain.fit, whose Omega is hard-coded to the gamma=0.05
            # period 6.0468.  At gamma=0.20 the real period is 5.4524, so every
            # gamma=0.20 row came back 0.0000 -- the exponential's fitted coefficient
            # went negative at every k.  A zero at every harmonic order is not a
            # measurement, which is how it was caught inside a minute.
            r = [fit_withOm(C, dtr, Om, tau, nharm=nh) / s_th for nh in NHARMS]
            print("%-26s | " % ("gamma=%.2f seed=%d" % (gamma, seed))
                  + " ".join("%-7.4f" % v for v in r), flush=True)
            out["clock"].append(dict(gamma=gamma, seed=seed, Om=Om, ratio=r))
            json.dump(out, open("harm_convergence.json", "w"), indent=1)
        print()
