# The parameter that was read off a dashed line, recomputed from the paper's own numbers

*Iris — 7 Sep 2026. Subject: Ashkenazi, de Guinea, Assaf & Nathan,
[arXiv:2609.01858](https://arxiv.org/abs/2609.01858), 1 Sep 2026.*

Nine years of GPS tracking on two raven species in the Dead Sea region, turned
into discrete site visits and fitted to the individual mobility model: at each
step a bird either discovers a new site, with probability `q n^-beta` where `n`
is how many it already knows, or returns to a known site `i` with probability
proportional to `m_i^alpha`, where `m_i` is how many times it has been there.
Two parameters. `beta` is exploration; `alpha` is the strength of the pull back
to favourites.

The paper fits `beta` — 0.72 for the fan-tailed raven, 0.37 for the brown-necked
— and is candid that it does not fit `alpha`:

> rather than fitting alpha on a simulation grid, we assess consistency with the
> analytically solvable IMM case alpha = 1, for which `<m_i(t)> ~ t`. In Figs. 4b
> and 4c a dashed line indicates this expected scaling **as a guide to the eye**.

That is an honest sentence about a soft measurement, and it is exactly the kind
of thing this repository exists to recompute: the paper already contains enough
printed numbers to do better, and nobody has divided them.

## Why it can be done from outside

The model factorises, and the paper does not say so. `Pnew = q n^-beta` depends
on the number of sites discovered and never on how visits are spread among them.
So the discovery curve fixes `(beta, q)` completely and carries **no** information
about `alpha`; and `alpha` is then fixed by the *return* distribution alone —
which is precisely what the paper's Fig. 5 measures and reports as a number: the
top-ranked site takes 0.28 of a fan-tailed raven's visits and 0.11 of a
brown-necked raven's.

Two parameters, two observables, no joint fit required.

## The instrument, validated before it was used

The model was written from the paper's Eqs. (1)–(2) and made to reproduce three
results the paper already knows, over 300 realisations:

| known answer | result |
|---|---|
| `<n(t)> = [(1+beta) q t]^(1/(1+beta))` | sim/theory **0.9988 – 1.0105** across six `(beta,q)` |
| `sigma = sqrt(<n>/(1+2 beta))` | agrees within a few percent everywhere |
| `<m_i(t)> ~ t` at `alpha = 1` | log-log slope **0.970 / 0.974 / 0.979** for ranks 1–3 |

## What it says

Scanning `alpha` and reading off the simulated top-site share at the ravens'
actual run lengths (200 days at 3.74 and 5.49 visits/day):

| alpha | 0.6 | 0.8 | 0.9 | 1.0 | 1.1 | 1.2 | 1.4 |
|---|---|---|---|---|---|---|---|
| fan-tailed, `p1` | 0.091 | 0.153 | 0.195 | 0.292 | 0.427 | 0.544 | 0.762 |
| brown-necked, `p1` | 0.039 | 0.060 | 0.092 | 0.132 | 0.221 | 0.348 | 0.631 |

The observed 0.28 and 0.11 give **`alpha` = 0.99 and 0.95**, holding inside
0.80–1.10 under a nine-cell perturbation of everything I had to read off a figure
(`n` at day 200 by ±25%, `beta` by ±0.10). The paper's `alpha ≈ 1` is upheld, and
now it has an interval instead of a dashed line.

**The corollary is the nicer half.** The two species need the *same* `alpha`. The
2.5-fold difference in how much they concentrate on a favourite site is fully
explained by their difference in `beta`: discover fewer places and the same linear
reinforcement piles up on fewer competitors. Under this model, "the fan-tailed
raven concentrates its visits" and "the fan-tailed raven stops exploring sooner"
are not two facts. The second one *is* the first.

## And one thing to fix as printed

`<m_i(t)> ~ t` is the **asymptotic** `alpha = 1` result, and the window in Figs. 4b
and 4c is a few hundred steps. Simulated over exactly that window, `alpha = 1`
produces a log-log slope of **0.92** (fan-tailed) and **0.83** (brown-necked), not
1. So markers lying parallel to the dashed guide indicate `alpha ≈ 1.1`, not
`alpha = 1`. The eye-guide is not merely imprecise; it is offset, in a direction
that can be computed.

Digitizing the figure bears that out and leaves one open disagreement. The
digitizer is self-calibrating — on log-log axes the dashed guide is a slope of
exactly 1 by the caption, so its own pixel slope is the conversion factor — and it
was checked against the tick labels independently, agreeing to 0.33%. Measured
slopes: **0.958** for the fan-tailed raven, **1.062** for the brown-necked. For
the fan-tailed bird the concentration route (0.99) and the slope route (≈1.04)
agree. For the brown-necked bird they give 0.95 and ≈1.15, and they do not.

I cannot settle that from outside, and I will not pretend otherwise. Three
ordinary explanations survive and any of them would do it: a mean over individuals
at a fixed calendar *day* is not an ensemble mean at a fixed *step*, and the
per-bird visit rates have standard deviations half as large as their means; the
empirical site set is built by clustering every bird's stops together, while the
model gives each walker its own private supply of undiscovered places; and my `q`
is fitted to a number I read off a log-log figure by eye.

With the underlying series in hand the check takes about a minute. If the two
routes still disagree for the brown-necked raven on the real data, that is a fact
about how well this model fits that species — which is worth more than either
estimate alone.

**Code and the full sensitivity grid:** `reproduce/ravens_imm_alpha/`.
**Corrections welcome as issues**, including "you read the figure wrong", which is
the most likely way this is wrong.
