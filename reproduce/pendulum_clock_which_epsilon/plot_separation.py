import json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = json.load(open("separate_gamma_eps.json"))
A = [r for r in rows if r["tag"] == "sweep-eps"]
B = [r for r in rows if r["tag"] == "sweep-gam"]
try:
    V = json.load(open("vdp_control.json"))
except FileNotFoundError:
    V = []

fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.4))

# ---- A: gamma fixed, eps swept
e = np.array([r["eps"] for r in A]); qa = np.array([r["Q"] for r in A])
ra = np.array([r["Rstar"] for r in A]); g0 = A[0]["gamma"]
ax[0].loglog(e, qa, "o-", color="#1b6ca8", ms=7, label=r"measured $Q_\infty$")
ax[0].axhline(g0**2/2, color="#c1121f", ls="--", lw=2,
              label=r"$\gamma^2/2$   (damping reading)")
ax[0].loglog(e, e**2/2, ":", color="#555", lw=2,
             label=r"$\varepsilon^2/2$   (nonlinearity reading)")
ax[0].set_xlabel(r"escapement strength  $\varepsilon$   (16$\times$)")
ax[0].set_ylabel(r"$Q_\infty$")
ax[0].set_title("A.  hold the damping, sweep the escapement\n"
                r"$\gamma=%.3g$,  $R^*$ %.2f$\to$%.2f" % (g0, ra.min(), ra.max()),
                fontsize=10)
ax[0].legend(fontsize=8, loc="upper left"); ax[0].grid(alpha=.25, which="both")

# ---- B: eps fixed, gamma swept
g = np.array([r["gamma"] for r in B]); qb = np.array([r["Q"] for r in B])
rb = np.array([r["Rstar"] for r in B]); e0 = B[0]["eps"]
ax[1].loglog(g, qb, "o-", color="#1b6ca8", ms=7, label=r"measured $Q_\infty$")
ax[1].loglog(g, g**2/2, "--", color="#c1121f", lw=2, label=r"$\gamma^2/2$")
ax[1].axhline(e0**2/2, color="#555", ls=":", lw=2, label=r"$\varepsilon^2/2$")
ax[1].set_xlabel(r"damping  $\gamma$   (16$\times$)")
ax[1].set_title("B.  hold the escapement, sweep the damping\n"
                r"$\varepsilon=%.3g$,  $R^*$ %.2f$\to$%.2f" % (e0, rb.min(), rb.max()),
                fontsize=10)
ax[1].legend(fontsize=8, loc="upper left"); ax[1].grid(alpha=.25, which="both")

# ---- C: the residual collapses on eps/gamma, and vanishes for van der Pol
for S, mk, c, lab in ((A, "o", "#1b6ca8", "grasshopper, sweep $\\varepsilon$"),
                      (B, "s", "#e07a1f", "grasshopper, sweep $\\gamma$")):
    x = np.array([r["eps"]/r["gamma"] for r in S])
    y = np.array([r["Q"]/(r["gamma"]**2/2) for r in S])
    o = np.argsort(x)
    ax[2].semilogx(x[o], y[o], mk+"-", color=c, ms=7, label=lab)
if V:
    x = np.array([r["eps"]/r["gamma"] for r in V])
    y = np.array([r["Q"]/(r["gamma"]**2/2) for r in V])
    o = np.argsort(x)
    ax[2].semilogx(x[o], y[o], "^-", color="#2a7f3e", ms=7,
                   label=r"van der Pol ($\bar f_\Phi\equiv 0$)")
ax[2].axhline(1.0, color="#c1121f", ls="--", lw=2)
ax[2].set_xlabel(r"$\varepsilon/\gamma$")
ax[2].set_ylabel(r"$Q_\infty \,/\, (\gamma^2/2)$")
ax[2].set_title("C.  what is left over is not the nonlinearity:\n"
                "it collapses on $\\varepsilon/\\gamma$, and it is the\n"
                "$R\\to R^*$ approximation, not the escapement", fontsize=10)
ax[2].legend(fontsize=8); ax[2].grid(alpha=.25, which="both")
ax[2].set_ylim(0.8, 3.2)

fig.suptitle(r"Which $\varepsilon$ is the pendulum-clock uncertainty law about?"
             "   (arXiv:2609.04957,  $T=0.005$, 2000 trajectories/point)", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig("which_epsilon.png", dpi=150)
print("wrote which_epsilon.png")
