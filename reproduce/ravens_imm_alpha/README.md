# ravens — arXiv:2609.01858, rebuilt

The PDF (11.7 MB) is not kept; `paper.txt` is `pdftotext -layout` of it and
`hi-06.png` is page 6 at 340 dpi, which is what `digitize2.py` reads. Re-fetch
with `https://arxiv.org/pdf/2609.01858v1` if a figure other than 4 is needed.

    imm.py          the model, written from Eqs. (1)-(2)
    validate.py     the three known answers it must return first
    alpha_scan.py   p1 vs alpha, and q fitted from the discovery curve
    sensitivity.py  alpha under +/-25% in n(200d) and +/-0.10 in beta
    convention.py   both readings of "rank-i site", to show the slope barely cares
    digitize2.py    edges not centroids; self-calibrated on the dashed guide
    axis_calib.txt  independent check of that calibration off the tick labels
    curvature.py    local slopes at the figure's own sample times

Findings and caveats: the finding note in `../../findings/2026-09_ravens-imm-alpha.md`.
