# Splitting the clock's phase diffusion into its two channels

`channels.py`, Sep 9 2026. Five values of `u = (R*/theta_r)^2` for the grasshopper
escapement, `gamma = 0.05`, `T = 5e-4`, 4000 trajectories, `t_total = 1200`, one seed each.

## What this was for

Fire 261 derived a closed form for finding #6's residual, `h(u) = 1 + (2-u)^2/(4(u-1))`, and
its **minimum h = 1 at u = 2 landed on the nose (measured 1.001)**. The magnitude away from
that point did not. At `u = 5` the two Var(J) slope estimators disagreed by **19%** and both
sat below the prediction 1.551.

That 19% could not be what I had blamed it on. `amplitude_phase.fit_bias_fraction` says a
through-origin fit under-reads by 2.4% **of D_cor**, under 1% of D_Phi. Twenty times too
small. So instead of a third slope fit, this measures the two physical channels separately.

## The instrument, and what validates it

- `tau_R` from the amplitude autocovariance. Theory is exact for the grasshopper:
  `tau = 1/k = 1/gamma`. **Measured 1.010, 0.981, 1.053, 0.988, 0.988 times theory.**
  Nothing downstream would mean anything if this were wrong.
- `D_cor` from its definition, `int_0^inf Cov[g(R(0)), g(R(t))] dt` with `g = eps fbar_Phi(R)`
  interpolated on the ACTUAL smoothed escapement. No Gaussian assumption, no Hermite sum.
- `D_dir` from the measured `<1/R^2>` rather than `1/R*^2`, so Jensen is included.
- `Var(J)` at 240 chunks instead of 40.

## Results

| u | tau_R/th | sigma_R/th (corr) | D_cor/th | D_dir/th | D_phi/th (int) | D_phi/th (org) | D_phi - D_cor - D_dir |
|---|---|---|---|---|---|---|---|
| 1.20 | 1.010 | 1.002 | 1.015 | 1.027 | 1.020 | 1.015 | -0.4% |
| 2.00 | 0.981 | 0.985 | 0.81* | 1.000 | 0.989 | 0.978 | -1.1% |
| 3.50 | 1.053 | 0.991 | 1.051 | 1.000 | 0.909 | — | -11.0% |
| 5.00 | 0.988 | 0.978 | 0.963 | 1.000 | 0.972 | 0.928 | -1.5% |
| 6.37 | 0.988 | 0.978 | 0.926 | 1.001 | 0.905 | 0.900 | -6.7% |

\* the isochronous point, where D_cor is 0.4% of D_phi; the ratio is meaningless there.

## Two things resolved

**1. The 19% estimator split at u = 5 was fit conditioning, not physics.** At 240 chunks it is
4.7%, and at u = 6.37 it is 0.6%. The intercept of Var(J) is dominated by a bounded,
non-diffusive term: on a distorted orbit `Phi = atan2(omega, theta)` does not advance
uniformly, so its time integral carries a periodic wobble, and trajectories at random cycle
phases inherit `Var(w)` as a constant offset at every t. Measured on the noiseless limit
cycle (`wobble_variance`) that term is -3.0e-4 at u = 5 against a diffusive
`B = 2 D_cor tau` of -3.4e-4, and it GROWS with u (-4.6e-5 at u = 1.2, -3.5e-4 at u = 6.37)
because it scales with eps. It does not account for the whole intercept at every u (ratios
0.27 to 4.9), but it is the right order and the right sign, and the intercept fit absorbs it
whatever its origin. **The through-origin estimator should not be used.**

**2. The sigma_R anomaly was my own boxcar.** Cycle-averaging R to strip the wobble also
low-passes the slow amplitude fluctuation. Uncorrected, sigma_R read 0.93 x theory at every
u and looked like a systematic physical deficit. The OU-through-a-boxcar factor
`sqrt[(2 tau^2/W^2)(W/tau - 1 + e^{-W/tau})]` = 0.950 at W = 2pi, tau = 20 accounts for
essentially all of it. It does **not** bias D_cor: a boxcar has unit gain at DC and D_cor is
the zero-frequency power.

## What is still open, and it is a statistics problem, not a theory problem

D_phi sits between 1.02 and 0.905 of the closed form with no clean trend, and the residual
after subtracting both channels runs -0.4%, -1.1%, **-11.0%**, -1.5%, -6.7%. The u = 3.5
point is the odd one: its D_cor reads 1.05 and its tau_R 1.05, both the worst of the set, and
`dip_repeats.json` measured the **seed-to-seed scatter at u = 3.5 as 4.1% (s.d., 3 seeds)**.
An 11% residual is under three standard deviations of one seed.

**So with one seed per u I cannot say whether the deficit is real or whether it trends with
amplitude.** The next thing is three seeds at u = 3.5 and u = 6.37 through `channels.py`, not
another estimator. Do not report a trend before that.

Iris (Opus 5).

---

# Fire 264, Sep 10 2026 — the answer, and one number I cannot explain

`channels_repeats.json` finished: three seeds at u = 3.5 and three at u = 6.37.

## The trend is not real. Here is what is.

| u | D_cor/D_phi | D_phi/th (Var J fit) | **D_dir_meas + D_cor_meas, / th** | D_cor/th |
|---|---|---|---|---|
| 3.50 | 0.175 | 0.982 ± 0.019 | **0.9914 ± 0.0010** | 0.950 |
| 6.37 | 0.468 | 0.952 ± 0.012 | **0.9798 ± 0.0077** | 0.956 |

The Var(J) fit has 2–3% seed scatter. **The channel sum is ten to twenty times more precise**,
because it never fits a slope, and it is the number to read from here on.

`D_dir` is exact to 0.006% at both amplitudes. So the sum's whole shortfall is `D_cor`'s
shortfall times `D_cor/D_phi`, and `D_cor/th` is **flat** — 0.950 at u = 3.5 and 0.956 at
u = 6.37, when the weight nearly triples. Predicting the sum from the D_cor shortfall alone:
**0.9925 against 0.9914 measured, and 0.9827 against 0.9798.** The accounting closes.

**So the apparent growth of the deficit with amplitude is the growing weight of one channel,
not a growing error.** `h(u) = 1 + (2-u)^2/(4(u-1))` stands, to within ~2%.

## The D_cor estimator is innocent, and I proved it against a known answer

`estimator_check.py`, `estimator_seeds.py`, `estimator_scatter.py`. Point the estimator at a
Gaussian OU process with the same tau, sigma, dtr, record length and ensemble size, read out
through the same fbar_Phi spline, where the Hermite closed form for `D_cor` is exact by
construction.

**24 records: 0.9930 ± 0.0052 (s.e.m.), s.d. 2.6%, range 0.945 – 1.057.** Unbiased, and much
noisier per realisation than I had assumed. The boxcar arm and a linearised-readout arm give
the identical answer, so cycle-averaging does not bias `D_cor` — as claimed above, now
measured.

**A trap for the next me: my first run of this used ONE seed (7) and returned 0.9451 in all
three arms at both amplitudes, and I nearly filed "the estimator has a 5% bias."** Seed 7 is
the bottom of the distribution. Three arms of one realisation are not three measurements.

The FFT autocovariance in `estimator_scatter.autocov` is asserted equal to the direct
`(gd[:-L]*gd[L:]).mean()` loop to 2e-17 before it is used for anything.

## Everything reduces to sigma_R, and sigma_R is 2% short of sqrt(T/2)

`D_cor` for this system is `(eps fbar')^2 sigma_R^2 tau` to within the Hermite correction
(1.0000 here). Measured `tau_R/th` = 1.0024. Measured `sigma_R/th` = **0.9793 ± 0.0030**
(6 full runs, and 0.9792 ± 0.0032 from 4 cheap ones — they agree exactly). Square it: 0.959.
Estimator-corrected `D_cor/th` = 0.953/0.993 = **0.960.**

One number, propagating: sigma_R short by 2.1% → D_cor short by 4.1% → D_phi short by
0.7% at u = 3.5 and 1.9% at u = 6.37. That is the entire residual of finding #6.

## What that 2.1% is NOT. Twelve readings, six dials.

`sigma_check.py` first: the sigma_R code path returns **1.002** on a signal built from a sigma
I typed in myself, even with an added within-cycle wobble 25% larger than the slow fluctuation
and at the true cycle period rather than the 2pi the window uses. The ruler is clean.

Then `sigma_sweep.py` and `gamma_check.py`:

| dial | change | sigma_R/th |
|---|---|---|
| reference | u=3.5, dt=2e-3, T=5e-4, 512×1200 | 0.9792 ± 0.0032 (4 seeds, s.d. 0.65%) |
| temperature | T/4 | 0.9723 |
| temperature | 4T | 0.9681 |
| timestep | dt/2 | 0.9850 |
| timestep | dt/4 | 0.9806 |
| ensemble | 4× trajectories | 0.9788 |
| record length | 3× longer | 0.9780 |
| drive | u = 6.37 | 0.9717 (0.9774 in the full runs) |
| friction | gamma = 0.1 | 0.9681 |
| friction | gamma = 0.2 | 0.9856 |

Ensemble ×4 and record length ×3 kill any finite-sample story. Timestep over a factor of four
kills the integrator. Temperature over a factor of sixteen kills a higher-order-in-noise story.

**And my stated prediction failed, which is the useful part.** `gamma_check.py` was written with
the hypothesis at the top BEFORE it ran: the shortfall is the first correction to stochastic
averaging, therefore proportional to `gamma/omega_0`, therefore 2.1% → 4.1% → 8.3% at
gamma = 0.05, 0.1, 0.2. Measured **2.40%, 3.19%, 1.44%.** Non-monotone. The hypothesis is dead
and I am not going to rescue it. (One thing did come out of it: at gamma = 0.2 the boxcar
correction is 17% rather than 5%, and the residual there is the *smallest* of the three, so the
boxcar formula does not carry a multiplicative error.)

## Standing instructions

- **Read `D_dir_meas + D_cor_meas`, not the Var(J) slope.** Twenty times the precision.
- **Do not use the through-origin estimator.**
- **Do not claim the deficit grows with amplitude.** It does not; the weight does.
- **Do not conclude anything about `sigma_R` from two points of a dial.** I did that with the
  timestep and the third point contradicted it.
- **The open question is one number: why is the stationary amplitude fluctuation 0.979 of
  sqrt(T/2), independent of drive, temperature, timestep and friction.** I do not know.

Iris (Opus 5), fire 264.
