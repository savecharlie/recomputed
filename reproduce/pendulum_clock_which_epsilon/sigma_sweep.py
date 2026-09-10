"""sigma_R reads 0.979 x sqrt(T/2).  Which dials move it?  (And how noisy is the reading?)

WHY THIS FILE EXISTS AS A SWEEP AND NOT AS A CONCLUSION: dt_check.py printed 0.9699 at
dt=2e-3 and 0.9850 at dt=1e-3 and I wrote "clear O(dt) convergence" in my own notes before
the third point arrived at 0.9806 and killed it.  Two points make a line whatever they are.
So: four seeds at the reference setting to get the reading's own scatter, then one dial at a
time -- timestep, noise temperature, ensemble size -- against that scatter.

Iris (Opus 5), fire 264.
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from dt_check import run

REF = dict(u=3.5, dt=2e-3)
rows = []


def rec(tag, **kw):
    a = dict(REF); a.update(kw)
    c, raw, dtr = run(**a)
    rows.append(dict(tag=tag, ratio=c, raw=raw, **a))
    print("  %-22s %s  ->  sigma_R/th = %.4f" % (tag, {k: v for k, v in kw.items()}, c), flush=True)
    return c


print("reference scatter (u=3.5, dt=2e-3, T=5e-4):")
base = [rec("seed", seed=s) for s in (5, 6, 7, 8)]
b = np.array(base)
print("  mean %.4f  s.d. %.4f (%.2f%%)  s.e.m %.4f\n"
      % (b.mean(), b.std(ddof=1), 100*b.std(ddof=1)/b.mean(), b.std(ddof=1)/2))

print("dials:")
for T in (1.25e-4, 2e-3):
    rec("temperature", T=T)
for dt in (1e-3, 5e-4):
    rec("timestep", dt=dt)
rec("ensemble", n_traj=2048)
rec("record length", t_total=3600.0)
rec("high amplitude", u=6.37)
json.dump(rows, open('sigma_sweep.json', 'w'), indent=1)
