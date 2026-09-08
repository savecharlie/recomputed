"""Which (alpha, G) conditions from Fig. 3 actually appear in the Fig. 4 collapse?

Fig. 4 markers are FILLED from a colour bar running alpha = 0.4 (bottom) to 0.8
(top), so a fill colour decodes to alpha and the x position gives alpha^2/G;
the two together give the G at which each collapse point was taken.
Calibration is certified by the paper's own printed constants: reading the drawn
fit line back out returns p0 = 0.00375 (printed 0.0037) and C1 = 27.7 (26.8).
"""
import numpy as np, collections
from PIL import Image
from scipy import ndimage as ndi

a = np.asarray(Image.open("Figure4.png").convert("RGB")).astype(int)
X0, PX = 205.0, 179.33
to_x = lambda px: (px - X0) / PX * 0.05

cb_x = 1008
s_ = a[:, cb_x, :]; v = (s_.max(-1) - s_.min(-1))
nz = np.nonzero(v > 60)[0]
ys = np.arange(nz.min() + 5, nz.max() - 4)
lut = a[ys, cb_x, :]
alph = 0.8 - (ys - ys.min()) / (ys.max() - ys.min()) * 0.4

sub = a[16:664, 200:990]                      # plot interior only
flat = sub.reshape(-1, 3)
d = np.abs(flat[:, None, :] - lut[None, :, :]).sum(-1)
best = d.min(1); who = d.argmin(1)
m = (best < 40).reshape(sub.shape[:2])
who = who.reshape(sub.shape[:2])
lab, n = ndi.label(m)
sizes = ndi.sum(m, lab, range(1, n + 1))
cents = ndi.center_of_mass(m, lab, range(1, n + 1))

by = collections.defaultdict(list)
kept = 0
for i, (s, (cy, cx)) in enumerate(zip(sizes, cents), 1):
    if s < 120:
        continue
    al = float(np.median(alph[who[lab == i]]))
    al = round(al * 10) / 10
    by[al].append(to_x(cx + 200)); kept += 1
print(f"blobs {n}, markers kept {kept}\n")
print(f"{'alpha':>5} {'n':>3}   alpha^2/G span        implied G span")
for al in sorted(by):
    xs = np.array(by[al])
    print(f"{al:>5} {len(xs):>3}   {xs.min():.3f} .. {xs.max():.3f}    "
          f"G = {al*al/xs.max():.2f} .. {al*al/xs.min():.2f}")
