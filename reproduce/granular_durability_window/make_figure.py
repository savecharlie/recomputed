import json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

P = json.load(open("fig3_points.json"))
canon = np.array([0.41, 0.89, 1.40, 2.40, 3.07])
al = np.array([p["alpha"] for p in P])
G  = canon[np.argmin(np.abs(np.array([p["G"] for p in P])[:, None] - canon), 1)]
N  = np.array([p["Nf"] for p in P])
x  = al**2 / G

COL = {0.4:"#e8467c", 0.5:"#eb9a3c", 0.6:"#2fa14b", 0.7:"#1f8fc4", 0.8:"#6b4fbb"}
fig, ax = plt.subplots(1, 2, figsize=(13.2, 5.4))

# ---- left: everything on the collapse axis --------------------------------
a0 = ax[0]
a0.axvspan(0.037, 0.205, color="#c9c9c9", alpha=.45, zorder=0)
a0.text(0.121, 1.6e3, "the window Fig. 4 actually uses\n"r"$G = 2.4-4.1$",
        ha="center", va="top", fontsize=9.5, color="#555")
xs = np.linspace(0.01, 0.42, 50)
a0.plot(xs, (1/0.0037)*np.exp(26.8*xs), "k-", lw=1.6,
        label=r"Eq. 2 as published:  $N_f=p_0^{-1}e^{\,26.8\,\alpha^2/G}$")
hi = G >= 2.4
for a in sorted(COL):
    m = (al == a)
    a0.scatter(x[m & hi], N[m & hi], s=52, color=COL[a], edgecolor="k",
               lw=.7, zorder=3, label=rf"$\alpha={a}$")
    a0.scatter(x[m & ~hi], N[m & ~hi], s=62, facecolor="none",
               edgecolor=COL[a], lw=2.0, zorder=3, marker="s")
a0.set_yscale("log"); a0.set_xlim(0, 0.43); a0.set_ylim(2e2, 2e7)
a0.set_xlabel(r"$\alpha^2/G$", fontsize=13)
a0.set_ylabel(r"$N_f$  (cycles to failure)", fontsize=13)
a0.set_title("Fig. 3's data on Fig. 4's axis\n"
             r"filled = $G\geq2.4$ (in the collapse)   open squares = $G\leq1.4$ (not)",
             fontsize=10.5)
a0.legend(fontsize=8.4, loc="upper left", framealpha=.95)
a0.grid(alpha=.25, which="both", lw=.4)

# ---- right: sensitivity to 1/G, family by family --------------------------
a1 = ax[1]
req, got = [], []
for a in sorted(COL):
    m = al == a
    a1.scatter(1/G[m], N[m], s=46, color=COL[a], edgecolor="k", lw=.6, zorder=3)
    k, b = np.polyfit(1/G[m], np.log(N[m]), 1)
    t = np.linspace((1/G[m]).min(), (1/G[m]).max(), 20)
    a1.plot(t, np.exp(k*t+b), color=COL[a], lw=2.0, zorder=2)
    ref = np.exp(26.8*a*a*t + (np.log(N[m]).mean() - 26.8*a*a*(1/G[m]).mean()))
    a1.plot(t, ref, color=COL[a], lw=1.1, ls=":", zorder=1)
    req.append(26.8*a*a); got.append(k)
    a1.annotate(rf"$\alpha={a}$", (t[-1], np.exp(k*t[-1]+b)), fontsize=9,
                color=COL[a], xytext=(4, -2), textcoords="offset points")
a1.set_yscale("log"); a1.set_xlabel(r"$1/G$", fontsize=13)
a1.set_ylabel(r"$N_f$", fontsize=13); a1.set_ylim(2e2, 1e8)
a1.set_title("solid = measured slope inside each $\\alpha$;  dotted = the slope\n"
             r"Eq. 2 requires ($C_1\alpha^2$, same intercept). Every family is too flat.",
             fontsize=10.5)
a1.grid(alpha=.25, which="both", lw=.4)
for a, r, g in zip(sorted(COL), req, got):
    print(f"alpha={a}: required {r:6.2f}   measured {g:6.2f}   ratio {r/g:5.1f}")

fig.suptitle("arXiv:2609.00587 — the collapse is fitted on $G\\geq2.4$; the paper's own "
             "Fig. 3 goes down to $G=0.41$", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig("collapse_window.png", dpi=140)
print("wrote collapse_window.png")
