"""Reanalysis of the open compression data from Abdellahi, Rakowska, Treder & Lewis,
"Targeted memory reactivation elicits temporally compressed reactivation linked to
spindles," Imaging Neuroscience 4:IMAG.a.1123 (3 Feb 2026), doi:10.1162/IMAG.a.1123.

Their data (CC BY 4.0), straight from the paper's own GitHub:
  new_compressionRatio.mat   19 temporal ratios  (sleep trial duration / wake trial)
  new_compression_scores.mat 19 x 48 correct classification rates, chance = 0.25

GATE FIRST: reproduce the paper's own significance pattern.  The abstract and
Discussion both say reactivation is "3 to 20 times faster" than wake, tested by
Wilcoxon signed-rank against chance at each ratio.  If that set falls out of their
numbers with my code, the instrument is trustworthy and may then say something new.

THEN THE NEW THING.  Their Discussion argues the compressed reactivation could be
carried by hippocampal ripples, "characterized by 50 to 100 ms of high-frequency
activity," on the grounds that 3-20x compression of a 1.15 s wake trial gives
57-383 ms.  But 1150/100 = 11.5: only compressions FASTER than ~11.5x are short
enough to fit inside one ripple.  383 ms is about four ripples.  So their range
spans two different physiological scales, and the question their own data can
answer is which one the evidence actually sits at.

Iris (Opus 5), Sep 9 2026.
"""
import numpy as np
import scipy.io as sio
from scipy.stats import wilcoxon

WAKE_MS = 1150.0          # paper: "the whole length of wake trial (1.15 second)"
RIPPLE_MS = (50.0, 100.0)  # paper's own bracket, citing Ylinen et al. 1995
CHANCE = 0.25

ratio = sio.loadmat('new_compressionRatio.mat')['new_compressionRatio'].squeeze()
S = sio.loadmat('new_compression_scores.mat')['new_compression_scores']   # 19 x 48
assert S.shape == (ratio.size, 48), S.shape
speed = 1.0 / ratio                 # x faster than wake
dur = ratio * WAKE_MS               # ms that the reactivation occupies

print("=" * 78)
print("GATE  --  reproduce the paper's own Wilcoxon test against chance (0.25)")
print("=" * 78)
print("%9s %8s %9s | %8s %8s %10s %s"
      % ("ratio", "xfaster", "dur(ms)", "mean", "median", "p", "sig"))
sig = []
for i, r in enumerate(ratio):
    x = S[i]
    stat, p = wilcoxon(x - CHANCE)
    s = p < 0.05
    sig.append(s)
    print("%9.4f %8.2f %9.1f | %8.4f %8.4f %10.2g %s"
          % (r, speed[i], dur[i], x.mean(), np.median(x), p, "*" if s else ""))
sig = np.array(sig)

fast_sig = speed[sig]
print("\nsignificant ratios span %.2fx to %.2fx faster than wake"
      % (fast_sig.min(), fast_sig.max()))
print("paper says: \"3 to 20 times faster\"")
lo, hi = round(fast_sig.min()), round(fast_sig.max())
print("GATE %s  (rounded: %dx to %dx)"
      % ("GREEN" if (lo, hi) == (3, 21) or (lo, hi) == (3, 20) else "CHECK", lo, hi))

# ---------------------------------------------------------------- the new part
print()
print("=" * 78)
print("NEW  --  ripple-scale or spindle-scale?  (the paper's Discussion claims ripple)")
print("=" * 78)
ripple_ok = dur <= RIPPLE_MS[1]
print("ratios whose reactivation fits inside ONE 50-100 ms ripple (dur <= 100 ms):")
print("   " + ", ".join("%.3f (%.0f ms, %.0fx)" % (ratio[i], dur[i], speed[i])
                        for i in np.where(ripple_ok)[0]))
print("ratios that are significant but TOO LONG for one ripple:")
print("   " + ", ".join("%.3f (%.0f ms, %.0fx)" % (ratio[i], dur[i], speed[i])
                        for i in np.where(sig & ~ripple_ok)[0]))

# where is the evidence strongest?  per-participant argmax over the significant band
print("\nEffect size by ratio (mean CCR - chance, and Cohen's d over 48 participants):")
for i in np.where(sig)[0]:
    x = S[i] - CHANCE
    print("   %.4f  %5.1fx  %6.0f ms   delta=%.4f   d=%.2f   %s"
          % (ratio[i], speed[i], dur[i], x.mean(), x.mean() / x.std(ddof=1),
             "ripple-scale" if ripple_ok[i] else "spindle-scale"))

best = np.argmax(S.mean(axis=1))
print("\nstrongest mean CCR at ratio %.4f  =  %.1fx faster  =  %.0f ms  (%s)"
      % (ratio[best], speed[best], dur[best],
         "fits one ripple" if dur[best] <= RIPPLE_MS[1] else "needs a spindle"))

# per-participant peak: does the population split into two scales?
pk = np.argmax(S, axis=0)
print("\nper-participant argmax ratio (n=48):")
vals, cnts = np.unique(pk, return_counts=True)
for v, c in zip(vals, cnts):
    print("   %.4f  %5.1fx %6.0f ms : %2d participants %s"
          % (ratio[v], speed[v], dur[v], c, "#" * c))
nrip = int(np.sum(dur[pk] <= RIPPLE_MS[1]))
print("\n%d/48 participants peak at a ripple-compatible ratio (<=100 ms);"
      " %d peak longer." % (nrip, 48 - nrip))
