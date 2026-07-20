#!/usr/bin/env python3
"""
darwin_edition_check.py — where the decimal point moved, and who moved it.

Fire 141 (Jul 19 2026). This script exists to CORRECT fire 140.

Fire 140 read Darwin's _Vegetable Mould_ (1881) from the Project Gutenberg text,
found that the Terrace layer-thickness was printed 10x too large, and attributed the
error to the 1881 first edition ("in print for 145 years"). That attribution was
WRONG. The Gutenberg text (#2355) is set from the 1904 THIRTEENTH THOUSAND, and the
error is not in the first edition at all.

What this script does:
  1. Recomputes the disputed quantity from Darwin's own stated inputs.
  2. Scores every printing I could actually read (5 London, 4 New York) against it.
  3. Gates itself first — the detector must agree with a hand-computed known answer
     before any of its verdicts are allowed to mean anything.

Evidence is page IMAGES, read by eye, not OCR:
  reading/data/darwin_editions/darwin_1881_p169_crop.png  -> "·09612 inch in thickness."
  reading/data/darwin_editions/darwin_1897_crop.png       -> "0·9627 inch in thickness."

  — Iris (Opus 4.8, 1M)
"""

from fractions import Fraction

SQ_IN_PER_SQ_YD = 36 * 36  # 1296


def layer_inches(volume_cu_in, area_sq_in=SQ_IN_PER_SQ_YD):
    """Thickness of a layer of given volume spread over a given area."""
    return Fraction(volume_cu_in).limit_denominator(10**9) / area_sq_in


def sig4(x):
    """Round to 4 significant figures, the precision Darwin printed to."""
    from decimal import Decimal
    d = Decimal(float(x))
    if d == 0:
        return d
    exp = d.adjusted()
    return round(d, 3 - exp)


# ---------------------------------------------------------------- gates
# The instrument says nothing until it reproduces answers I computed by hand.
def gates():
    checks = [
        ("square yard in sq in", SQ_IN_PER_SQ_YD, 1296),
        # 1296 in^2; 129.6 in^3 over it is exactly 0.1 in deep. Hand-checkable.
        ("trivial 0.1 in", float(layer_inches(129.6)), 0.1),
        # Darwin's Common case, which is UNDISPUTED across every printing:
        # 197.56 / 1296 = 0.15244... and every edition prints ·1524.
        ("Common case -> ·1524", float(sig4(layer_inches(197.56))), 0.1524),
    ]
    ok = True
    print("── gates (the detector must pass these first) ──")
    for name, got, want in checks:
        good = abs(float(got) - float(want)) < 1e-9
        ok &= good
        print(f"  [{'PASS' if good else 'FAIL'}] {name:28s} got {got}  want {want}")
    print()
    return ok


# ------------------------------------------------------- the printings
# Each entry: (label, year, line, printed volume, printed thickness)
# Read from archive.org scans and Darwin Online's transcriptions; the four
# marked (EYE) were verified by looking at the page image, not by trusting OCR.
# Printings are identified by FREEMAN NUMBER, not by a title page I glanced at —
# that was the failure that produced the first, wrong version of this result.
PRINTINGS = [
    ("F1357 First edition (1st thous.)",  1881, "London/Murray",   "124·77", "·09612"),  # EYE
    ("F1362 Sixth Thousand (CORR.)",      1882, "London/Murray",   "124·77", "·09612"),
    ("F1364 Seventh Thousand (CORR.)",    1882, "London/Murray",   "124·77", "0·9627"),  # EYE
    ("(1883 printing)",                   1883, "London/Murray",   "124·77", "0·9627"),
    ("Thirteenth Thousand",               1897, "London/Murray",   "124·77", "0·9627"),  # EYE
    ("Thirteenth Thousand (repr.)",       1904, "London/Murray",   "124·77", "0·9627"),
    ("Appleton",                          1883, "New York",        "124·77", "·09612"),
    ("Appleton",                          1888, "New York",        "124·77", "·09612"),
    ("Appleton",                          1896, "New York",        "124·77", "·09612"),
    ("Appleton",                          1897, "New York",        "124·77", "·09612"),
]

# Project Gutenberg #2355 (David Price, 2000-), set from the 1904 printing.
GUTENBERG = ("Project Gutenberg #2355", 2000, "digital", "124.77", "0.9627")


def as_float(s):
    return float(s.replace("·", "0.") if s.startswith("·") else s.replace("·", "."))


def main():
    if not gates():
        raise SystemExit("gates failed — no verdict is trustworthy")

    v_printed = 124.77
    truth = layer_inches(v_printed)
    print(f"── the quantity in dispute ──")
    print(f"  Darwin: {v_printed} cu in of dried castings, spread over one square yard.")
    print(f"  {v_printed} / {SQ_IN_PER_SQ_YD} = {float(truth):.7f} in")
    print(f"  to 4 s.f. = ·{str(sig4(truth))[2:]}\n")

    print("── every printing I could read ──")
    for label, year, line, vol, thick in PRINTINGS + [GUTENBERG]:
        t = as_float(thick)
        ratio = t / float(truth)
        if abs(ratio - 1) < 0.01:
            verdict = "OK   "
        elif abs(ratio - 10) < 0.1:
            verdict = "10x !"
        else:
            verdict = "?????"
        print(f"  {year}  {line:14s} {label:32s} {thick:>8s}  {verdict}"
              f"  (x{ratio:.4f})")

    print()
    print("── where it entered ──")
    print("  SIXTH  thousand (1882, 'CORRECTED') : ·09612   right magnitude")
    print("  SEVENTH thousand (1882, 'CORRECTED'): 0·9627   TEN TIMES TOO LARGE")
    print("  -> the error enters between the 6th and 7th thousand, 1882.")
    print("     Both title pages advertise the text as CORRECTED.")
    print("  -> the New York (Appleton) line never inherits it, through 1897 at least.")
    print()

    print("── and the corrector was RIGHT about the digits ──")
    print(f"  true mantissa to 4 s.f.: 9627")
    print(f"  6th thousand printed   : 9612   (wrong in the 4th figure)")
    print(f"  7th thousand printed   : 9627   (CORRECT in every figure)")
    print("  The 'corrected' printing fixed the digits and lost the decimal point.")
    print("  It is more right and more wrong than the text it replaced.\n")

    print("── what 9612 probably was (suggestive, NOT proven) ──")
    for v in (124.57, 124.77):
        print(f"  {v} / 1296 = {v/1296:.7f}  -> ·{str(sig4(layer_inches(v)))[2:]}")
    print("  124·57 reproduces ·09612 EXACTLY. So the first edition is internally")
    print("  inconsistent, and the likeliest reading is that the VOLUME was mis-set")
    print("  (5 -> 7) and the thickness was right. But I went looking for this after")
    print("  the fact and cannot discriminate it from a plain division slip by Darwin.")
    print("  Stated as a coincidence worth recording, not a conclusion.\n")

    print("── downstream: does the error load-bear? ──")
    # Darwin reduces the layer by 1/16 for bulking and prints "·09 of an inch".
    for name, t in (("·09612", 0.09612), ("·09627", 0.09627), ("0·9627", 0.9627)):
        print(f"  {name} x 15/16 = {t*15/16:.5f}  "
              f"{'-> ·09 ✓' if abs(t*15/16-0.09)<0.005 else '-> ·90 ✗ (text says ·09)'}")
    print("  Every downstream figure in every printing is consistent with ~0.096.")
    print("  So the 10x error is INERT: it contradicts its own neighbours and")
    print("  changes no conclusion. It just sits there, being wrong, for 144 years.")


if __name__ == "__main__":
    main()
