# lambda2_not_a_bound

Rebuild of **arXiv:2609.34352v1** (Naoto Shiraishi, 28 Sep 2026), read from the
LaTeX e-print rather than from extracted PDF text.

Run `python3 anchors.py` first. Eight anchors, all of which must pass before any
other script in this directory is allowed to say anything:

1. the generator's columns sum to zero and its off-diagonals are non-negative
2. the paper's stated stationary distribution of `Q`, for `Q0`, `Q1` and `Q`
3. the paper's closed form for the eigenvalues of `Q`, to 1e-12 at three sizes
4. the Fourier block decomposition against a dense diagonalisation of the full
   matrix, matched without a sort key so conjugate-pair ordering cannot fool it
5. the stationary distribution of the full `R` by two independent routes
6. `sigma_ss` by a dense O(n^2) double sum and by a shift-invariant edge sum
7. a known zero: no driving anywhere gives `sigma_ss = 0` and a real spectrum
8. **a positive control** — the OSB inequality holds on 25 plain driven rings.
   Without this, a violation detector that simply always fires looks identical
   to a correct one.

Then: `pcondition.py` (the misprint), `minsize.py` and `robust.py` (how small),
`verify.py` / `verify2.py` / `twolayer_verify.py` (60-digit witnesses),
`amplitude.py` and `repair.py` and `kscale.py` (the two amplitudes),
`candidate.py` and `falsify.py` (a candidate bound, attacked),
`twolayer.py` / `twolayer_min2.py` (two layers are enough at finite size),
`figure.py`.

Everything is seeded. Nothing here needs a network.
