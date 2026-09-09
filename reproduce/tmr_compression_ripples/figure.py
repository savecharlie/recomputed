"""Three panels on the open compression data of Abdellahi et al. 2026.

A: mean CCR at every temporal ratio they tested.  The curve is above chance
   EVERYWHERE, including at dilations no replay hypothesis predicts.
B: per-participant argmax, against the uniform expectation.
C: the test the paper did not run -- compressed band vs their own dilated band,
   paired within participant.

Iris (Opus 5), Sep 9 2026.
"""
import numpy as np, scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, chisquare

CH, WAKE, RIP = 0.25, 1150.0, 100.0
ratio = sio.loadmat('new_compressionRatio.mat')['new_compressionRatio'].squeeze()
S = sio.loadmat('new_compression_scores.mat')['new_compression_scores']
speed, dur = 1.0 / ratio, ratio * WAKE
ps = np.array([wilcoxon(S[i] - CH).pvalue for i in range(ratio.size)])
order = np.argsort(ps); bh = np.empty(ratio.size); run = 1.0
for r in range(ratio.size - 1, -1, -1):
    i = order[r]; run = min(run, ps[i] * ratio.size / (r + 1)); bh[i] = run
comp, dil = ratio < 0.95, ratio > 1.05

fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.2))
fig.suptitle("Abdellahi, Rakowska, Treder & Lewis 2026 (doi:10.1162/IMAG.a.1123), "
             "reanalysed from their own posted data  ·  n=48, 4-class, chance 0.25",
             fontsize=11)

# ---- A
m, se = S.mean(axis=1), S.std(axis=1, ddof=1) / np.sqrt(48)
ax[0].axhline(CH, color='k', lw=1, ls='--', label='chance (0.25)')
ax[0].fill_between(dur, CH, m, where=m > CH, color='#d8c8e8', alpha=.5)
ax[0].errorbar(dur, m, yerr=se, marker='o', ms=4.5, lw=1.4, color='#5b2d8e',
               label='mean CCR ± SEM')
ax[0].plot(dur[ps < .05], m[ps < .05], 'o', ms=9, mfc='none', mec='#c0392b', mew=1.8,
           label='raw p<0.05 (the paper: "3 to 20×")')
ax[0].plot(dur[bh < .05], m[bh < .05], 'o', ms=14, mfc='none', mec='#16a085', mew=1.8,
           label='survives BH FDR q<0.05 (5.6–15×)')
ax[0].axvline(RIP, color='#e67e22', lw=1.6)
ax[0].text(RIP * 0.93, CH + 0.0035, "one ripple\n(50–100 ms)", color='#e67e22',
           fontsize=9, va='bottom', ha='right')
ax[0].axvline(WAKE, color='gray', lw=1, ls=':')
ax[0].text(WAKE * 1.06, CH + .0015, "no compression", color='gray', fontsize=8,
           rotation=90, va='bottom')
ax[0].set_xscale('log')
ax[0].set_xlabel("duration of the reactivation (ms)   ←  more compressed")
ax[0].set_ylabel("correct classification rate")
ax[0].set_title("A.  above chance at every ratio tested", loc='left', fontsize=10)
ax[0].legend(fontsize=8, loc='lower right')
ax[0].grid(alpha=.25)

# ---- B
pk = np.argmax(S, axis=0)
cnt = np.bincount(pk, minlength=ratio.size)
cs = chisquare(cnt)
ax[1].bar(np.arange(ratio.size), cnt, color=['#5b2d8e' if c else '#b0b0b0' for c in comp])
ax[1].axhline(48 / ratio.size, color='#c0392b', lw=1.6,
              label='uniform expectation (2.5)')
ax[1].set_xticks(np.arange(ratio.size))
ax[1].set_xticklabels(["%.2g" % s if s >= 1 else "%.2g" % s for s in speed],
                      rotation=90, fontsize=7)
ax[1].set_xlabel("× faster than wake  (purple = compressed, grey = dilated)")
ax[1].set_ylabel("participants peaking here")
ax[1].set_title("B.  per-participant best ratio: $\\chi^2$=%.1f, df=18, p=%.2f"
                % (cs.statistic, cs.pvalue), loc='left', fontsize=10)
ax[1].legend(fontsize=8)
ax[1].grid(alpha=.25, axis='y')

# ---- C
a, b = S[comp].mean(axis=0), S[dil].mean(axis=0)
st, p = wilcoxon(a - b)
for j in range(48):
    ax[2].plot([0, 1], [b[j], a[j]], color='#999', lw=.7, alpha=.6)
ax[2].plot([0, 1], [b.mean(), a.mean()], color='#c0392b', lw=3, marker='o', ms=8,
           label='mean  %+.4f' % (a.mean() - b.mean()))
ax[2].axhline(CH, color='k', lw=1, ls='--')
ax[2].set_xticks([0, 1])
ax[2].set_xticklabels(["dilated band\n(9 ratios, 0.46–0.85×)",
                       "compressed band\n(9 ratios, 1.4–21×)"], fontsize=9)
ax[2].set_ylabel("participant mean CCR")
ax[2].set_title("C.  the test not run: paired, n=48.  Wilcoxon p = %.2f, d = %.2f"
                % (p, (a - b).mean() / (a - b).std(ddof=1)), loc='left', fontsize=10)
ax[2].legend(fontsize=9)
ax[2].grid(alpha=.25, axis='y')

plt.tight_layout(rect=[0, 0, 1, .95])
plt.savefig('compression_plateau.png', dpi=135)
print("wrote compression_plateau.png")
