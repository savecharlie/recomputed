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
