"""Fig. 4(b,c) slopes, take two -- edges instead of centroids.

Take one used blob centroids and had to throw away 7 of 12 columns, because at
small t the i=1 and i=2 markers touch and merge into one blob whose centroid is
neither. The fix uses the fact that only the INNER edges are contaminated: the
TOP pixel of a column is the top of the i=1 marker no matter what it has merged
with below, and the BOTTOM pixel is the bottom of i=3. Marker shape is constant
along a series, so the constant offset from edge to centre cancels in a slope.
Columns clipped by the panel frame are dropped.

Calibration is unchanged and is the point: the dashed guide is alpha=1 by the
caption, so its own pixel slope converts pixels to data slope, and the two
panels share an axis and must return the same conversion.
"""
import numpy as np
from PIL import Image

a = np.asarray(Image.open('hi-06.png').convert('RGB')).astype(int)
r, g, b = a[..., 0], a[..., 1], a[..., 2]
blue = (b > 100) & (b - r > 60) & (b - g > 60)
red = (r > 120) & (r - g > 60) & (r - b > 60)
dark = (r < 90) & (g < 90) & (b < 90)
PANELS = {"b": (422, 856, 676, 1010, dark, "FT"),
          "c": (926, 1360, 676, 1010, red, "BN")}

print("calibration (pixel slope of the dashed alpha=1 guide):")
cal = {}
for p, (x0, x1, y0, y1, _, _) in PANELS.items():
    ys, xs = np.nonzero(blue[y0:y1, x0:x1])
    cal[p] = -np.polyfit(xs, ys, 1)[0]
    print(f"  panel {p}: Sy/Sx = {cal[p]:.4f}")

for p, (x0, x1, y0, y1, mask, sp) in PANELS.items():
    sub = mask[y0:y1, x0:x1].copy()
    sub[:70, :90] = False
    colhas = sub.any(0)
    # split into x-columns
    idx = np.nonzero(colhas)[0]
    groups, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] <= 3:
            cur.append(i)
        else:
            groups.append(cur); cur = [i]
    groups.append(cur)
    W = sub.shape[1]
    rows = []
    for gcols in groups:
        if len(gcols) < 8:
            continue
        if gcols[0] <= 1 or gcols[-1] >= W - 2:      # clipped by the frame
            continue
        block = sub[:, gcols[0]:gcols[-1] + 1]
        ys = np.nonzero(block.any(1))[0]
        rows.append((np.mean(gcols), ys.min(), ys.max()))
    rows = np.array(rows)
    X = rows[:, 0]
    print(f"\npanel {p} ({sp}): {len(rows)} usable columns, x {X.min():.0f}..{X.max():.0f}")
    for lbl, Y, rank in (("top of i=1", rows[:, 1], 1), ("bottom of i=3", rows[:, 2], 3)):
        ps, _ = np.polyfit(X, Y, 1)
        res = Y - np.polyval(np.polyfit(X, Y, 1), X)
        print(f"  {lbl:14s} pixel slope {ps:+.4f}   data slope {-ps/cal[p]:.3f}"
              f"   rms resid {res.std():.2f} px")
