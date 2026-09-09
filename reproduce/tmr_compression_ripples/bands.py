"""Is the compression result a PEAK or a plateau?  Three tests their paper did not run.

The per-ratio Wilcoxon in compression.py reproduces the paper exactly.  But the
effect sizes are suspiciously flat across the whole significant band (Cohen's d
0.34-0.45) and the largest single bin of per-participant argmaxes is the MOST
DILATED ratio tested.  That pattern is what a uniformly elevated curve looks like
when you run 19 uncorrected tests along it, so it has to be checked.

Iris (Opus 5), Sep 9 2026.
"""
import numpy as np, scipy.io as sio
from scipy.stats import wilcoxon, spearmanr, chisquare

CH, WAKE = 0.25, 1150.0
ratio = sio.loadmat('new_compressionRatio.mat')['new_compressionRatio'].squeeze()
S = sio.loadmat('new_compression_scores.mat')['new_compression_scores']
speed = 1.0 / ratio
comp = ratio < 0.95
dil = ratio > 1.05
print("compressed band: %d ratios (%.0fx..%.1fx faster) | dilated band: %d ratios"
      % (comp.sum(), speed[comp].max(), speed[comp].min(), dil.sum()))

# ---- 1. multiple comparisons across the 19 ratios --------------------------
ps = np.array([wilcoxon(S[i] - CH).pvalue for i in range(ratio.size)])
order = np.argsort(ps)
m = ratio.size
bh = np.empty(m)
running = 1.0
for rank in range(m - 1, -1, -1):
    i = order[rank]
    running = min(running, ps[i] * m / (rank + 1))
    bh[i] = running
print("\n1. 19 uncorrected Wilcoxon tests, no correction reported in the paper.")
print("   raw p<0.05      : %d ratios" % (ps < 0.05).sum())
print("   Benjamini-Hochberg q<0.05 : %d ratios %s"
      % ((bh < 0.05).sum(), np.array2string(speed[bh < 0.05], precision=1)))
print("   Bonferroni p<0.05/19      : %d ratios %s"
      % ((ps < 0.05 / m).sum(), np.array2string(speed[ps < 0.05 / m], precision=1)))

# ---- 2. compressed vs dilated, paired within participant -------------------
a = S[comp].mean(axis=0)
b = S[dil].mean(axis=0)
st, p = wilcoxon(a - b)
d = (a - b).mean() / (a - b).std(ddof=1)
print("\n2. per-participant mean CCR, compressed band vs dilated band (n=48, paired)")
print("   compressed %.4f   dilated %.4f   difference %+.4f" % (a.mean(), b.mean(), a.mean() - b.mean()))
print("   Wilcoxon p = %.3f   Cohen's d = %.2f" % (p, d))

# ---- 3. does the score depend on ratio at all? -----------------------------
rho, prho = spearmanr(ratio, S.mean(axis=1))
print("\n3. mean CCR vs ratio across the 19 ratios: Spearman rho=%+.2f p=%.3f"
      % (rho, prho))
rs = [spearmanr(ratio, S[:, j]).correlation for j in range(48)]
st2, p2 = wilcoxon(rs)
print("   per-participant rho (CCR vs ratio): median %+.3f, sign test p=%.3f"
      % (np.median(rs), p2))

# ---- 4. are the per-participant peaks anywhere but uniform? ----------------
pk = np.argmax(S, axis=0)
cnt = np.bincount(pk, minlength=ratio.size)
cs = chisquare(cnt)
print("\n4. per-participant argmax histogram vs uniform over 19 ratios")
print("   counts: %s" % np.array2string(cnt))
print("   chi2 = %.1f, df = %d, p = %.2f  -> %s"
      % (cs.statistic, ratio.size - 1, cs.pvalue,
         "indistinguishable from uniform" if cs.pvalue > 0.05 else "structured"))
print("   peaks in the compressed half: %d/48 (uniform would give %.1f)"
      % (comp[pk].sum(), 48 * comp.sum() / ratio.size))

# ---- 5. the baseline itself -------------------------------------------------
print("\n5. the floor: mean CCR at the LEAST significant ratios")
for i in np.argsort(ps)[-5:]:
    print("   ratio %.4f (%.1fx) mean %.4f  delta %+.4f  p=%.2f"
          % (ratio[i], speed[i], S[i].mean(), S[i].mean() - CH, ps[i]))
print("   every one of the 19 ratios has mean CCR above chance; range of delta "
      "%.4f to %.4f" % ((S.mean(axis=1) - CH).min(), (S.mean(axis=1) - CH).max()))

# ---- 6. the framings MOST favourable to the paper ---------------------------
print("\n6. most-favourable framings (if any of these hold, the effect is real)")
sig7 = ps < 0.05
a7 = S[sig7].mean(axis=0)
st, p = wilcoxon(a7 - b)
print("   a) their 7 significant ratios vs the 9 dilated, paired: %+.4f  p=%.3f"
      % ((a7 - b).mean(), p))
ibest_c = int(np.argmax(S[comp].mean(axis=1)))
ibest_d = int(np.argmax(S[dil].mean(axis=1)))
xc, xd = S[comp][ibest_c], S[dil][ibest_d]
st, p = wilcoxon(xc - xd)
print("   b) best compressed ratio (%.1fx) vs best dilated (%.2fx), paired: %+.4f p=%.3f"
      % (speed[comp][ibest_c], speed[dil][ibest_d], (xc - xd).mean(), p))
fdr = bh < 0.05
afdr = S[fdr].mean(axis=0)
st, p = wilcoxon(afdr - b)
print("   c) the 4 FDR-surviving ratios vs the 9 dilated, paired: %+.4f  p=%.3f"
      % ((afdr - b).mean(), p))
st, p = wilcoxon(S[0] - S[-1])
print("   d) fastest tested (21x) vs slowest tested (0.46x), paired: %+.4f p=%.3f"
      % ((S[0] - S[-1]).mean(), p))
print("\n   NOTE both ways: the 19 ratios come from overlapping sliding windows over the")
print("   same sleep data, so they are strongly dependent.  That makes Bonferroni too")
print("   harsh AND makes the contiguity of the significant band much less impressive")
print("   than it looks.  The paired within-participant tests above do not depend on")
print("   that structure, which is why they are the ones that matter.")
