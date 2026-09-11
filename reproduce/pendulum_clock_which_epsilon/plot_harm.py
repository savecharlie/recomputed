"""The picture for the harmonic-convergence result: a truncation with a plateau.

Left panel is synthetic, where sigma was typed in by hand, so the correct answer is the
line at 1.0.  Right panel is the clock, where nothing is known in advance.  The point of
putting them side by side is that the STOPPING RULE -- take the plateau -- is chosen on the
left, where it can be checked, and only then applied on the right.

Iris (Opus 5), fire 266.
"""
import sys, os, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(here, "harm_convergence.json")))
NH = np.array(d["nharms"], float)

fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))

# ---- left: synthetic, answer known
for r in d["synthetic"]:
    if r["seed"] != 301:
        continue
    ax[0].plot(NH, r["ratio"], "o-", ms=4.5, lw=1.4,
               label="wobble %.1f$\\sigma$, sharpness %d" % (r["wobble"], r["sharp"]))
ax[0].axhline(1.0, color="k", lw=1.0, ls="--")
ax[0].text(13, 1.045, "the $\\sigma$ I typed in", fontsize=8.5)
ax[0].set_xlim(1, 28)
ax[0].set_yscale("log")
ax[0].set_title("synthetic: the answer is known", fontsize=11)
ax[0].set_xlabel("harmonics in the fit basis")
ax[0].set_ylabel("estimated $\\sigma$ / true $\\sigma$")
ax[0].legend(fontsize=8, frameon=False, loc="upper right")

# ---- right: the clock
cols = {0.05: "#1f77b4", 0.20: "#d62728"}
seen = set()
for r in d["clock"]:
    g = r["gamma"]
    lbl = None if g in seen else "$\\gamma$ = %.2f" % g
    seen.add(g)
    ax[1].plot(NH, r["ratio"], "o-", ms=4, lw=1.2, color=cols[g], alpha=0.85, label=lbl)
for g, pred in ((0.05, 0.9813), (0.20, 0.9315)):
    ax[1].axhline(pred, color=cols[g], lw=1.0, ls=":")
    ax[1].text(1.4, pred + 0.0016, "isostable prediction", fontsize=8, color=cols[g])
ax[1].annotate("sigma_nofilter.py\nstopped here",
               xy=(4, 0.9565), xytext=(6.4, 0.899), fontsize=8.5, color="#d62728",
               arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.0))
ax[1].annotate("...and 6 harmonics crosses the prediction\non the way past it, which is not agreement",
               xy=(6, 0.9317), xytext=(9.0, 0.9455), fontsize=8.5, color="#555555",
               arrowprops=dict(arrowstyle="->", color="#555555", lw=0.9))
ax[1].set_xlim(1, 30.5)
ax[1].set_title("the clock: nothing known in advance", fontsize=11)
ax[1].set_xlabel("harmonics in the fit basis")
ax[1].set_ylabel("measured $\\sigma_R$ / theory")
ax[1].legend(fontsize=9, frameon=False)

for a in ax:
    a.set_xticks([2, 4, 6, 8, 10, 12, 16, 20, 26])
    a.grid(alpha=0.22, lw=0.6)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)

fig.suptitle("A model of the nuisance always fits something: where the truncation stops mattering",
             fontsize=12)
fig.tight_layout()
out = os.path.join(here, "harm_convergence.png")
fig.savefig(out, dpi=150)
print("wrote", out)
