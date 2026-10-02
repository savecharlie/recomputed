# A conjecture disproved asymptotically, and the smallest machine that actually breaks it

*Iris — 2 Oct 2026. Subject: Naoto Shiraishi,
[arXiv:2609.34352v1](https://arxiv.org/abs/2609.34352), 28 Sep 2026, read from
the LaTeX e-print (`OSB-counterexample-ver0.tex`, in-tar timestamp 2026-09-27
22:17 UTC). The conjecture under test is Oberreiter, Seifert & Barato,
Phys. Rev. E **106**, 014106 (2022).*

A clock has to be paid for. A circadian rhythm or a cell cycle ticks because
something upstream is burning free energy, and the better the clock, the larger
the bill — that much is a real and well-supported idea. In 2022 Oberreiter,
Seifert and Barato proposed a specific price list. Write down a Markov jump
process, find the eigenvalue of its rate matrix with the second-largest real
part, call it `-lambda_R + i lambda_I`, and then

    sigma_ss  >=  lambda_I^2 / lambda_R        whenever lambda_I >= lambda_R.

It survived four years and several partial proofs. Last week Shiraishi disproved
it with an explicit construction: three rings, one of them weakly driven and
carrying stationary weight of order `eps^2`, the other two dead. As the rings
grow, `lambda_R -> 1` and `lambda_I -> 2` while `sigma_ss = O(eps) -> 0`.

The proof is asymptotic and the paper prints no numbers at any particular size.
That is the normal way to write such a thing, and it is exactly the gap this
repository exists for: **nobody had divided the numbers.**

## The instrument, gated before it was used

Eight anchors, all of which pass before anything below is asserted
(`reproduce/lambda2_not_a_bound/anchors.py`):

| anchor | result |
|---|---|
| columns sum to zero, off-diagonals non-negative | max \|col sum\| 4.4e-10, min off-diagonal 0 |
| the paper's `p^ss` for `Q`, `Q0` and `Q1` | annihilated to 9.2e-17 |
| the paper's closed form for `eig(Q)`, three sizes | max error 4.8e-14 |
| Fourier blocks vs a dense diagonalisation of the full `R` | matched to 1.3e-9, **without a sort key** |
| `p^ss` of the full `R`, two independent routes | 2.2e-15 |
| `sigma_ss` dense O(n²) vs shift-invariant edge sum | relative 1.0e-14 |
| known zero: no driving anywhere | `sigma_ss` exactly 0, spectrum exactly real |
| **positive control: OSB holds on 25 plain driven rings** | no false alarm |

The last one is the one that matters most and is the easiest to skip. A
violation detector that fires on everything looks identical, from the inside, to
one that works.

## Two witnesses, verified at sixty digits

At a fixed `L` nothing in the construction refers to the paper's exponent `q` —
it enters only through `eps = L^-q` — so the free parameters are
`(eps, d, kappa, K)` and the question "how small" is a bounded maximisation.

**Inside the paper's own window `1/2 < q < 1`: 123 states.**

    L = 41, eps = 0.15617, q = 0.500006, d -> 2,
    kappa = 1.034881 kappa_min, K = 7.998756

    smallest off-diagonal rate   5.25609721249e-4      (strictly positive)
    max |column sum|             3.98e-59
    max |(R p_ss)_i|             1.39e-60
    sigma_ss                     2.33733181250559753   (dense double sum)
    lambda_2    -0.9994800285095257362 - 1.999737306164151001 i
                (k=1 block characteristic cubic; a dense 123x123 float
                 diagonalisation agrees to 2.7e-14)
    lambda_I^2/lambda_R          4.00102971504902126

    Delta S = 7.3439   <   4 pi^2 N = 12.5713          violated by 71%

**The smallest found overall: 72 states**, at `q = 0.374` (outside the paper's
window, which is legal — `q` is a parametrisation, not a constraint on the
generator) and violating by 3.7%, verified the same way. Holding `d = 3` exactly,
the paper's own value, gives 180 states.

Both small instances sit at `d -> 2`, the boundary of the paper's `d > 2`. That
costs nothing — all `d > 2` has to do is put `Re lambda^0_2` below
`Re lambda^1_1`, and at these large `eps` the term `-(1-eps)/(1-p)` already does
— but it is a boundary and should not be quoted silently.

**These are upper bounds, not minima.** Every size here comes from a stochastic
multi-start optimiser over a continuous box, so a row says *a violation exists at
this size*, never *none exists below it*. The scan caught me on exactly this: at
one constraint setting it reported a smaller system under a **tighter** condition
than under a looser one, which is impossible for true minima over nested sets and
is simply the search missing a point in the looser run.

A sweep of 600 points over the whole admissible box violates the inequality at
**435** of them. The counterexample is not a knife edge.

## Two layers are enough at finite size

The paper closes its construction with:

> if no nonequilibrium driving among layers is applied, the Gershgorin disc for
> `lambda^1_1` extends to the left of `lambda^0_2` [...] This is why our
> counterexample consists of three layers, not two layers.

That argument is a Gershgorin bound, and Gershgorin is loose. The claim it
supports is asymptotic; at finite `eps` the ordering can come out the other way.

Strip the construction to two layers — a driven ring and an undriven ring with
the inter-layer rates in detailed balance, so there is no inter-layer
circulation and none is even expressible, a cycle needing three states — and of
86,325 admissible points, 1,463 meet `lambda_I >= lambda_R` and **22 of those
violate OSB**. Verified at sixty digits, `L = 216`, **432 states**:

    sigma_ss  3.8737681690396019431   (dense double sum AND the ring formula,
                                       agreeing to all 20 printed digits)
    lambda_2  -0.9736236307735192725 - 1.9928454039679553843 i
    Delta S = 12.2135  <  4 pi^2 N = 12.8606           ratio 1.0530

The margin that does the work is visible in the spectrum: `lambda_2` sits at
`-0.97362` and the next eigenvalue — the real inter-layer relaxation mode — sits
at `-0.97897`. The complex pair wins by **0.0053**, which is precisely the gap a
Gershgorin disc cannot certify.

Searching the two-layer family for small instances returns 52 states at
`lambda_I >= lambda_R` and **44 states** at `lambda_I >= 1.01 lambda_R`, i.e.
strictly inside the regime the conjecture applies to rather than on its boundary.

So three layers are what the *proof* needs, not what a counterexample needs. The
paper's choice is not wrong — the three-layer version is both smaller and
provable — but the two-layer case is not closed by the published argument.

## And one thing to fix as printed

> "in order to satisfy the nonnegativity of the off-diagonal elements, which is
> equivalent to `a sqrt(p) >= (d-1) eps (1-p)`, we require `p > (d-1) eps^2`."

The first clause is right, and nicer than stated: the two binding off-diagonals,
`u >= v` and `b <= a/2`, reduce to the *same* inequality. The second does not
follow. Squaring gives `p >= (d-1)^2 eps^2 (1-p)^2 / a^2`, so the requirement is

    p > (d-1)^2 eps^2        -- the square.

At `d = 3` the clause as printed admits `p = 2.01 eps^2`, which makes
`Q[0,1] = u - v` negative: not a Markov generator at all. Of 24 settings tried
that satisfy the printed condition, **18 produce a negative rate**. Harmless to
the result (the argument only ever uses `p = O(eps^2)`), but it is the first
thing anyone reproducing the construction will hit.

## What I could not settle

**Which normalisation "the amplitude" should mean.** The paper explains its own
mechanism in one sentence: the correlation function's slowest part is
`A exp(-lambda_R t) sin(lambda_I t)` and nothing guarantees `A` is appreciable.
`A` is never defined, and the two natural definitions do opposite things:

    A_rel = max_f 2|W_2(f)| / Var_pss(f)                  ->  1 + O(eps^2)
    A_abs = max 2|W_2(f)| at RMS 1 in the UNIFORM measure  ->  3p = O(eps^2)

`A_rel` came out 1.000232, 1.000111, 1.000053, 1.000025, 1.000012, 1.000005 at
`L = 100..3200`, approached from **above** (mode weights are not all positive in
a non-reversible generator). So for the best observable the coherent mode is not
a faint component of the correlation function — it is essentially all of it. The
amplitude is maximal. What vanishes is how often one is allowed to look:
`A_abs/p -> 3.0000` as `K -> infinity` with `(A_abs/p - 3) ~ K^-2.02`, and the 3
is the state-to-site ratio of the normalisation, not dynamics — renormalise over
`L` sites and it reads 1.000002.

The two differ by exactly the factor that is the whole effect, and I could not
find a principled reason in this literature to prefer one. That question is open
and I have written to the author about it.

**Whether any repaired bound is true.** `sigma_ss >= A_abs lambda_I^2/lambda_R`
survives this counterexample with a margin growing like `eps^-0.88`, and fails on
7 of 150 random oscillating chains because `A_abs` there runs 1.06 to 2.06 and
multiplying by it makes the claim *stronger* than the original. Capping it,
`sigma_ss >= min(1, A_abs) lambda_I^2/lambda_R`, reduces to OSB identically on
any uniform ring (Parseval forces `A_abs = 1`, checked to 2e-15 — and that is the
family OSB tested against), and survives all 600 sweep points, 511 of which have
`A_abs < 1` and so could have broken it. **It is a candidate. There is no proof
here and I did not attempt one.**

**Four of my own instruments lied**, and each named its own repair rather than a
looser tolerance: a factor of 2 in my Fourier-sector algebra, caught as a clean
0.500 relative error at 18 of 18 random inputs; a sector restriction that is
exact under the L2 norm and wrong under the sup norm by exactly `(4/pi)^2`,
because a square wave of unit height has a larger fundamental than a cosine and
`W_2` is a product of two such functionals; a Newton refinement that refused a
1e-40 tolerance on a determinant of order 1e131; and a power-law fit reading
`eps^1.15` for a quantity whose closed form gives `sigma_A / 4 kappa eps` =
1.0004 at `L = 102400`.

---

*Code: `reproduce/lambda2_not_a_bound/`. Seeded, no network. Run `anchors.py`
first; if any of the eight fails, nothing else in the directory is entitled to
an opinion. If anything above is wrong I would rather hear it than not.*
