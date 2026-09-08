"""Digitize Fig. 3 of arXiv:2609.00587v1 (Yokota & Kurita) -> (alpha, G, N_f) triples.

Markers are OPEN rings sharing their hue with the guide curves; the two are
separable by exact RGB (the curve uses a second shade that never reaches the
legend box).  Overlapping markers merge under connected components, so centres
come from MATCH FILTERING with a per-family ring template.

Axis calibration, checked two independent ways before anything else:
  x  major ticks at 148, 415, 683, 950  -> decade widths 267/268/267 (0.4% spread)
  y  major ticks 616..26 in 8 steps     -> 84/85/84/84/85/84/84 px per 0.5 G
"""
import numpy as np, json
from PIL import Image
from scipy import ndimage as ndi

a = np.asarray(Image.open("Figure3.png").convert("RGB")).astype(int)
H, W, _ = a.shape

X0, DEC = 148.0, 267.333          # x of 1e2, px per decade
Y0, PGY = 616.0, 168.571          # y of G=0, px per unit G
to_N = lambda x: 10 ** (2 + (x - X0) / DEC)
to_G = lambda y: (Y0 - y) / PGY

FAM = {0.4: (255, 0, 0), 0.5: (0, 58, 255), 0.6: (254, 0, 205),
       0.7: (13, 205, 0), 0.8: (0, 169, 255)}
LEGEND = (150, 260, 330, 590)     # x0,x1,y0,y1 of the legend box (generous)

def mask_of(col, tol=90):
    return (np.abs(a - np.array(col)).sum(-1) < tol)

def template_from_legend(m):
    """The legend marker is the one clean, isolated instance of each glyph."""
    sub = np.zeros_like(m)
    x0, x1, y0, y1 = LEGEND
    sub[y0:y1, x0:x1] = m[y0:y1, x0:x1]
    lab, n = ndi.label(sub)
    assert n == 1, f"legend gave {n} components, expected 1"
    sl = ndi.find_objects(lab)[0]
    t = (lab[sl] == 1).astype(float)
    return t, sl

out = []
report = {}
for al, col in FAM.items():
    m = mask_of(col)
    tmpl, sl = template_from_legend(m)
    th, tw = tmpl.shape
    # correlate the OUTLINE with the outline mask; normalise by template mass
    corr = ndi.correlate(m.astype(float), tmpl, mode="constant")
    corr /= tmpl.sum()
    # blank the legend so it is not reported as data
    x0, x1, y0, y1 = LEGEND
    corr[y0:y1, x0:x1] = 0
    # non-maximum suppression on a footprint slightly smaller than the glyph
    fp = np.ones((max(3, th - 6), max(3, tw - 6)), bool)
    peak = (corr == ndi.maximum_filter(corr, footprint=fp)) & (corr > 0.62)
    lab, n = ndi.label(peak)
    cents = ndi.center_of_mass(peak, lab, range(1, n + 1))
    pts = [(cx, cy) for cy, cx in cents]
    report[al] = dict(template=[th, tw], legend_slice=str(sl), found=len(pts))
    for cx, cy in pts:
        out.append(dict(alpha=al, x=round(cx, 1), y=round(cy, 1),
                        G=round(to_G(cy), 3), Nf=round(to_N(cx), 1)))

for al in FAM:
    pts = [p for p in out if p["alpha"] == al]
    gs = sorted({round(p["G"], 2) for p in pts})
    print(f"alpha={al}: {len(pts):2d} markers  G rows ~ {gs}")
json.dump(out, open("fig3_points.json", "w"), indent=1)
print("total", len(out))
