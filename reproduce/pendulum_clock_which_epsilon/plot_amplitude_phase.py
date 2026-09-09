"""Two panels: the isochronous dip, and the diffusion floor.

A: h = D_Phi / D_dir against u = (R*/theta_r)^2.  The closed form
   h = 1 + (2-u)^2/(4(u-1)) has a MINIMUM OF EXACTLY ONE at u = 2, where the
   escapement is isochronous at its own running amplitude.  Nothing in the paper
   or in finding #6 predicts a minimum anywhere.

B: hold eps, drop gamma.  The paper's D_Phi = D/(2R*^2) says phase diffusion should
   fall as gamma^2 over this sweep (256x).  Adding amplitude-to-phase conversion
   says it falls only as gamma, toward the floor gamma T/(8 theta_r^2).

Iris (Opus 5), Sep 9 2026.
"""
import json, sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import amplitude_phase as ap

dip = json.load(open('isochronous_dip.json'))
sat = json.load(open('saturation.json'))
TR, T = 0.4, 5e-4

fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.4))

# ---------------------------------------------------------------- A: the dip
uu = np.linspace(1.02, 6.0, 800)
ax[0].plot(uu, ap.gh_h(uu), color='#5b2d8e', lw=2,
           label=r'closed form  $1+\dfrac{(2-u)^2}{4(u-1)}$  (sharp escapement)')
ax[0].axhline(1.0, color='k', ls='--', lw=1, label='the paper: $h=1$, no amplitude term')
ax[0].axvline(2.0, color='#e67e22', lw=1.6)
ax[0].text(2.15, 1.58,
           "isochronous amplitude\n$R^*=\\sqrt{2}\\,\\theta_r$,  $d\\Omega/dR=0$",
           color='#e67e22', fontsize=9, va='top')
u_m = np.array([r['u'] for r in dip])
h_m = np.array([r['h_meas'] for r in dip])
h_p = np.array([r['h_pred'] for r in dip])
ax[0].plot(u_m, h_p, 's', ms=8, mfc='none', mec='#16a085', mew=1.8,
           label='prediction on the SMOOTHED force (Gaussian-OU readout)')
h_o = np.array([r['D_meas_origin'] / r['D_dir'] for r in dip])
ax[0].plot(u_m, h_o, 'o', ms=10, mfc='none', mec='#c0392b', mew=1.5,
           label='measured, through-origin slope fit')
ax[0].plot(u_m, h_m, 'o', ms=9, color='#c0392b',
           label='measured, intercept fit (4000 traj/point)')
try:
    rep = json.load(open('dip_repeats.json'))
    for uu_ in sorted({r['u'] for r in rep}):
        hs = [r['h_meas'] for r in rep if r['u'] == uu_]
        ax[0].errorbar([uu_], [np.mean(hs)], yerr=[np.std(hs, ddof=1) if len(hs) > 1 else 0],
                       fmt='D', ms=6, color='#7f8c8d', capsize=4,
                       label='repeat seeds, spread' if uu_ == sorted({r['u'] for r in rep})[0] else None)
except Exception:
    pass
for r in dip:
    ax[0].annotate("%.3f" % r['h_meas'], (r['u'], r['h_meas']), textcoords='offset points',
                   xytext=(0, -16), ha='center', fontsize=8, color='#c0392b')
ax[0].set_xlabel(r'$u=(R^*/\theta_r)^2$   (amplitude$^2$ in units of the critical angle)')
ax[0].set_ylabel(r'$h = D_\Phi\,/\,[\gamma T/(2R^{*2})]$')
ax[0].set_title("A.  the uncertainty penalty vanishes at the isochronous amplitude",
                loc='left', fontsize=11)
ax[0].legend(fontsize=7.6, loc='upper right')
ax[0].grid(alpha=.25)
ax[0].set_ylim(0.9, max(2.0, h_m.max() * 1.25))

# ---------------------------------------------------------------- B: the floor
g = np.array([r['gamma'] for r in sat])
Dm = np.array([r['D_meas'] for r in sat])
Dd = np.array([r['D_dir'] for r in sat])
Dp = np.array([r['D_pred'] for r in sat])
floor = T / (8 * TR * TR) * g
ax[1].loglog(g, Dm, 'o-', ms=9, color='#c0392b', lw=1.6, label=r'measured $D_\Phi$')
ax[1].loglog(g, Dp, 's--', ms=8, mfc='none', mec='#16a085', mew=1.8,
             label='prediction with amplitude-to-phase')
ax[1].loglog(g, Dd, '^:', ms=8, color='#5b2d8e',
             label=r"the paper's $\gamma T/(2R^{*2})$ alone")
ax[1].loglog(g, floor, color='#e67e22', lw=1.4, ls='-.',
             label=r'floor $\gamma T/(8\theta_r^2)$  (amplitude cancels)')
ax[1].set_xlabel(r'$\gamma$   (escapement strength held at $\epsilon=0.05$)')
ax[1].set_ylabel(r'$D_\Phi$')
ax[1].invert_xaxis()
ax[1].set_title("B.  swinging harder stops helping: $D_\\Phi$ falls as $\\gamma$, not $\\gamma^2$",
                loc='left', fontsize=11)
ax[1].legend(fontsize=8.5)
ax[1].grid(alpha=.25, which='both')
for r in sat:
    ax[1].annotate("$R^*$=%.2f" % r['Rstar'], (r['gamma'], r['D_meas']),
                   textcoords='offset points', xytext=(-4, 12), fontsize=8, ha='right')

plt.tight_layout()
plt.savefig('amplitude_to_phase.png', dpi=135)
print("wrote amplitude_to_phase.png")
sp = np.polyfit(np.log(g), np.log(Dm), 1)[0]
print("measured d log D_Phi / d log gamma = %.3f   (floor says 1, the paper says 2)" % sp)
print("D_meas / D_pred: " + "  ".join("%.3f" % (a / b) for a, b in zip(Dm, Dp)))
print("D_meas / D_dir : " + "  ".join("%.2f" % (a / b) for a, b in zip(Dm, Dd)))
