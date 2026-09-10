"""Draw the accounting so I have to LOOK at it: does one flat D_cor shortfall explain the curve?

Left: D_phi/theory against u, three ways -- the Var(J) slope fit (noisy), the channel sum
(precise), and the curve predicted by a SINGLE amplitude-independent D_cor shortfall.
Right: D_cor/theory against u, flat, against the estimator's own null band.

Iris (Opus 5), fire 264.
"""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rep = json.load(open('channels_repeats.json'))
est = np.array([r[1] for r in json.load(open('estimator_seeds.json'))])
null, null_sem = est.mean(), est.std(ddof=1)/np.sqrt(len(est))

us = sorted({r['u'] for r in rep})
g = {u: [r for r in rep if abs(r['u']-u) < 1e-9] for u in us}
w = {u: g[u][0]['D_cor_th']/g[u][0]['D_phi_th'] for u in us}
fit = {u: np.array([r['D_phi_meas']/r['D_phi_th'] for r in g[u]]) for u in us}
tot = {u: np.array([(r['D_dir_meas']+r['D_cor_meas'])/r['D_phi_th'] for r in g[u]]) for u in us}
dc  = {u: np.array([r['D_cor_meas']/r['D_cor_th'] for r in g[u]]) for u in us}

dc_all = np.concatenate([dc[u] for u in us])
short = null - dc_all.mean()                      # one number, amplitude-independent

fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.4))
fig.patch.set_facecolor('#faf8f4')
for a in ax:
    a.set_facecolor('#faf8f4')
    for s in ('top', 'right'):
        a.spines[s].set_visible(False)
    a.axhline(1.0, color='#999', lw=.8, zorder=0)

uu = np.linspace(1.0, 7.0, 300)
# weight of D_cor interpolated log-linearly through the two measured points
lw_ = np.interp(uu, us, [np.log(w[u]) for u in us])
ax[0].plot(uu, 1.0 - np.exp(lw_)*short, color='#2b6cb0', lw=2,
           label='one flat $D_{cor}$ shortfall (%.1f%%) $\\times$ its weight' % (100*short))
for u in us:
    ax[0].errorbar(u-.06, fit[u].mean(), yerr=fit[u].std(ddof=1)/np.sqrt(len(fit[u])),
                   fmt='o', ms=6, color='#b06a2b', capsize=3,
                   label='Var(J) slope fit' if u == us[0] else None)
    ax[0].errorbar(u+.06, tot[u].mean(), yerr=tot[u].std(ddof=1)/np.sqrt(len(tot[u])),
                   fmt='s', ms=7, color='#1a1a1a', capsize=3,
                   label='$D_{dir}+D_{cor}$, measured' if u == us[0] else None)
ax[0].set_xlabel('$u$'); ax[0].set_ylabel('$D_\\Phi$ / closed form')
ax[0].set_title('the deficit does not grow. the weight does.', loc='left', fontsize=11)
ax[0].set_xlim(2.6, 7.2); ax[0].legend(frameon=False, fontsize=8.5, loc='lower left')

ax[1].axhspan(null-2*null_sem, null+2*null_sem, color='#2b6cb0', alpha=.13,
              label='estimator null on a known $D_{cor}$ ($\\pm2$ s.e.m.)')
ax[1].axhline(null, color='#2b6cb0', lw=1.2)
for u in us:
    ax[1].plot([u]*len(dc[u]), dc[u], 'o', ms=5, color='#b06a2b', alpha=.75)
    ax[1].errorbar(u, dc[u].mean(), yerr=dc[u].std(ddof=1)/np.sqrt(len(dc[u])),
                   fmt='s', ms=7, color='#1a1a1a', capsize=3)
sig = 0.9793
ax[1].axhline(sig**2, color='#2f7a4f', ls='--', lw=1.4,
              label='$(\\sigma_R^{meas}/\\sigma_R^{th})^2 = %.3f$' % sig**2)
ax[1].set_xlabel('$u$'); ax[1].set_ylabel('$D_{cor}$ / closed form')
ax[1].set_title('and it is $\\sigma_R^2$, at both amplitudes', loc='left', fontsize=11)
ax[1].set_xlim(2.6, 7.2); ax[1].legend(frameon=False, fontsize=8.5, loc='lower left')
fig.tight_layout()
fig.savefig('accounting.png', dpi=155)
print("wrote accounting.png")
print("  one flat D_cor shortfall = %.4f  (null %.4f - measured %.4f)" % (short, null, dc_all.mean()))
for u in us:
    print("  u=%-5.2f w=%.3f   fit %.4f   sum %.4f   predicted %.4f   D_cor %.4f"
          % (u, w[u], fit[u].mean(), tot[u].mean(), 1-w[u]*short, dc[u].mean()))
