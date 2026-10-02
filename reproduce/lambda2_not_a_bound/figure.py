"""Three panels, each plotting a CLAIM rather than the quantity it came from.

(a) the violation: sigma_ss falls to zero while lambda_I^2/lambda_R sits at 4.
    Compensated inset: sigma_ss / (4 kappa eps) -> 1, which is the claim
    "sigma_ss = O(eps)" made into a horizontal line.
(b) the two amplitudes: A_rel -> 1 (the coherent mode is the WHOLE correlation
    function of the best observable) while A_abs -> 0 like the stationary
    weight p.  Compensated: A_abs / 3p -> 1.
(c) the 72-state witness: the actual spectrum of the actual generator.
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ.setdefault(v, "2")
import sys, json, numpy as np
sys.path.insert(0, ".")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import model as M, amplitude as A, repair as Rp

Ls = [100, 200, 400, 800, 1600, 3200, 6400, 12800]
eps, sig, rhs, arel, aabs, pp, kap = [], [], [], [], [], [], []
for L in Ls:
    P = M.rates(L, q=0.55)
    r = M.osb(P)
    eps.append(P["eps"]); sig.append(r["sigma"]); rhs.append(r["rhs"])
    arel.append(A.best_A(P, n=20)[0]); aabs.append(Rp.A_abs(P)[0])
    pp.append(P["p"]); kap.append(P["kappa"])
eps, sig, rhs = map(np.array, (eps, sig, rhs))
arel, aabs, pp, kap = map(np.array, (arel, aabs, pp, kap))

fig = plt.figure(figsize=(13.2, 4.3))
fig.patch.set_facecolor("white")
INK, ACC, WARN = "#1b1b1b", "#1f6fb2", "#c0392b"

# ---- (a) ------------------------------------------------------------------
ax = fig.add_subplot(1, 3, 1)
ax.plot(eps, rhs, "o-", color=ACC, ms=4.5, lw=1.6,
        label=r"$\lambda_I^2/\lambda_R$  (what OSB demands)")
ax.plot(eps, sig, "s-", color=WARN, ms=4.5, lw=1.6,
        label=r"$\sigma^{\rm ss}$  (what the system pays)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.invert_xaxis()
ax.set_xlabel(r"$\epsilon = L^{-0.55}$   (left = larger system)")
ax.set_ylabel("rate")
ax.set_title("(a)  the violation", loc="left", fontsize=11, color=INK)
ax.legend(fontsize=8, frameon=False, loc="lower left")
ax.grid(alpha=0.18, lw=0.6)
axi = ax.inset_axes([0.52, 0.52, 0.44, 0.30])
axi.plot(eps, sig / (4 * kap * eps), "s-", color=WARN, ms=3, lw=1.2)
axi.axhline(1.0, color=INK, lw=0.8, ls="--")
axi.set_xscale("log"); axi.invert_xaxis()
axi.set_ylim(0.9, 1.45)
axi.tick_params(labelsize=6)
axi.set_title(r"$\sigma^{\rm ss}/4\kappa\epsilon \to 1$", fontsize=7)
axi.grid(alpha=0.15, lw=0.4)

# ---- (b) ------------------------------------------------------------------
ax = fig.add_subplot(1, 3, 2)
ax.plot(eps, 3 * pp, "-", color="#9ecae1", lw=5.0, solid_capstyle="round",
        zorder=1, label=r"$3p$   (stationary weight of layer A)")
ax.plot(eps, arel, "o-", color=ACC, ms=4.5, lw=1.6, zorder=3,
        label=r"$A_{\rm rel}$ — relative to the observable's own variance")
ax.plot(eps, aabs, "s-", color=WARN, ms=4.5, lw=1.6, zorder=4,
        label=r"$A_{\rm abs}$ — at fixed size in the uniform measure")
ax.axhline(1.0, color=INK, lw=0.7, ls="--", alpha=0.5, zorder=0)
ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
ax.set_xlabel(r"$\epsilon = L^{-0.55}$")
ax.set_ylabel("amplitude of the coherent mode")
ax.set_title("(b)  two things called 'the amplitude'", loc="left",
             fontsize=11, color=INK)
ax.legend(fontsize=7.5, frameon=False, loc="lower left",
          bbox_to_anchor=(0.0, 0.52))
ax.set_ylim(1.5e-4, 3.0)
ax.grid(alpha=0.18, lw=0.6)
axi = ax.inset_axes([0.60, 0.07, 0.37, 0.25])
axi.plot(eps, aabs / (3 * pp), "s-", color=WARN, ms=3, lw=1.2)
axi.axhline(1.0, color=INK, lw=0.8, ls="--")
axi.set_xscale("log"); axi.invert_xaxis(); axi.set_ylim(0.995, 1.012)
axi.tick_params(labelsize=6)
axi.set_title(r"$A_{\rm abs}/3p \to 1$", fontsize=7)
axi.grid(alpha=0.15, lw=0.4)

# ---- (c) ------------------------------------------------------------------
ax = fig.add_subplot(1, 3, 3)
hit = next(r for r in json.load(open("minsize.json")) if r["ratio"] > 1)
P = M.rates_eps(hit["L"], hit["eps"], d=hit["d"], kappa=hit["kappa"], K=hit["K"])
w = np.linalg.eigvals(M.R_full(P))
r = M.osb(P)
ax.scatter(w.real, w.imag, s=22, color=INK, alpha=0.55, zorder=3)
l0 = M.np.linalg.eigvals(M.Q_matrix(P))          # the k=0 block: Q itself
ax.scatter(l0.real, l0.imag, s=70, marker="x", color="#2a7f3e", lw=1.6,
           zorder=5)
ax.scatter([-r["lR"], -r["lR"]], [r["lI"], -r["lI"]], s=190,
           facecolors="none", edgecolors=WARN, lw=1.8, zorder=4)
ax.axvline(0, color=INK, lw=0.6, alpha=0.5)
ax.axvline(-r["lR"], color=WARN, lw=0.7, ls=":", alpha=0.8, zorder=1)
ax.set_xlim(-2.15, 0.45); ax.set_ylim(-3.6, 3.6)
ax.set_xlabel(r"${\rm Re}\,\lambda$   (zoom: the fast modes run to $-15$)")
ax.set_ylabel(r"${\rm Im}\,\lambda$")
from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([], [], ls="", marker="o", mfc="none", mec=WARN, mew=1.6, ms=9,
           label=r"$\lambda_2 = %.4f \pm %.4f\,i$   (the $k{=}1$ mode in layer A)"
                 % (-r["lR"], r["lI"])),
    Line2D([], [], ls="", marker="x", color="#2a7f3e", mew=1.6, ms=7,
           label=r"eigenvalues of $Q$ — the DISSIPATIVE mode, at $-1.123$"),
    Line2D([], [], ls="", marker="o", color=INK, alpha=0.55, ms=5,
           label="the other 69 eigenvalues"),
], fontsize=7.0, frameon=False, loc="lower left", handletextpad=0.6,
   borderaxespad=0.3)
ax.grid(alpha=0.18, lw=0.6)
ax.annotate("the oscillation is here,\nand costs almost nothing",
            xy=(-r["lR"] + 0.04, r["lI"] + 0.12), xytext=(-0.62, 3.05),
            fontsize=7, color=WARN, ha="left",
            arrowprops=dict(arrowstyle="->", color=WARN, lw=0.8,
                            connectionstyle="arc3,rad=0.2"))
ax.annotate("the dissipation is here,\nand is NOT $\\lambda_2$",
            xy=(float(sorted(l0.real)[1]), 0.12), xytext=(-2.08, 2.65),
            fontsize=7, color="#2a7f3e", ha="left",
            arrowprops=dict(arrowstyle="->", color="#2a7f3e", lw=0.8,
                            connectionstyle="arc3,rad=0.25"))
ax.set_title("(c)  the smallest witness: 72 states\n"
             r"$\Delta S = %.3f \;<\; 4\pi^2\mathcal{N} = %.3f$ — the conjecture is false by %.1f%%"
             % (2 * np.pi * r["sigma"] / r["lI"],
                4 * np.pi ** 2 * r["lI"] / (2 * np.pi * r["lR"]),
                100 * (r["ratio"] - 1)),
             loc="left", fontsize=11, color=INK)

fig.suptitle("arXiv:2609.34352v1 rebuilt — what the second eigenvalue does and "
             "does not bound", fontsize=11.5, color=INK, y=1.005)
fig.tight_layout()
fig.savefig("lambda2_not_a_bound.png", dpi=155, bbox_inches="tight",
            facecolor="white")
print("wrote lambda2_not_a_bound.png")
print(f"panel (a) eps range {eps.min():.2e}..{eps.max():.2e}, "
      f"sigma {sig.min():.3e}..{sig.max():.3e}, rhs {rhs.min():.4f}..{rhs.max():.4f}")
print(f"panel (b) A_rel {arel.min():.6f}..{arel.max():.6f}, "
      f"A_abs {aabs.min():.3e}..{aabs.max():.3e}")
print(f"panel (c) {3*P['L']} states, ratio {r['ratio']:.6f}")
