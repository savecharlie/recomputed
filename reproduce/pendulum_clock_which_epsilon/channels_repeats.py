"""Three seeds at u = 3.5 and u = 6.37 through the channel split.

CHANNELS_FINDINGS.md ends on a statistics problem, not a theory problem: one seed per u,
and dip_repeats.json puts the seed-to-seed scatter at u = 3.5 at 4.1% (s.d., 3 seeds), which
is the same size as the residuals being chased.  This is the run that decides whether the
D_phi deficit and its apparent growth with u are real.  Nothing is to be claimed about a
trend before this finishes.

Iris (Opus 5), Sep 9 2026.  Detached with nohup so it survives the fire that started it.
"""
import sys, os, time, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper
import channels as ch

T, GAM, TR = 5e-4, 0.05, 0.4
esc = Grasshopper(TR, 20.0)
rows = []
for u in (3.5, 6.37):
    eps = math.pi * GAM * TR * u / 4.0
    for sd in (2001, 2002, 2003):
        t0 = time.time()
        r = ch.run_channels(esc, GAM, eps, T, seed=sd)
        o = ch.analyse(r, esc)
        o['wall_s'] = time.time() - t0
        rows.append(o)
        json.dump(rows, open('channels_repeats.json', 'w'), indent=1)
        print("  [%.0fs]\n" % o['wall_s'], flush=True)
for u in (3.5, 6.37):
    v = [r['D_phi_int'] / r['D_phi_th'] for r in rows if abs(r['u'] - u) < 1e-6]
    if len(v) > 1:
        print("u=%.2f  D_phi/th = %s  mean %.4f  s.d. %.4f (%.1f%%)  s.e.m %.4f"
              % (u, " ".join("%.4f" % x for x in v), np.mean(v), np.std(v, ddof=1),
                 100 * np.std(v, ddof=1) / np.mean(v), np.std(v, ddof=1) / np.sqrt(len(v))),
              flush=True)
