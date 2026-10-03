# Docking reporting fields: the tail ordering, and one ligand name that decayed

Reproduces two findings about **arXiv:2609.37542v1** (Giap Duc Ha, 29 Sep 2026),
*Which reported inputs govern molecular docking reproducibility?*

The `data/` directory is the author's own redistributed tables, copied verbatim
from the preprint's ancillary archive. No docking was run here.

    python3 anchors.py      # 65 anchors against the author's published values
    python3 tails.py        # the tail ranking, intervals, clustering
    python3 reresolve.py    # live PubChem re-resolution of all 116 ligand strings
    python3 cited.py        # the few numbers quoted straight from his files

`anchors.py` needs only the standard library. `tails.py` needs `scipy`.
`reresolve.py` makes 63 throttled requests to PubChem PUG-REST and takes about
40 seconds; **its output is a function of the date**, which is the point.

Full write-up: `findings/2026-10_docking-reporting-tails.md`.
