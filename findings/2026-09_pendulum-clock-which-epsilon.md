# A universal law with one symbol doing three jobs — and the one it is actually about

**Text:** Yuki Izumida, *Universal Thermodynamic Law Governing Stochastic Pendulum
Clocks*, arXiv:2609.04957v1, posted 4 September 2026.
**Checked:** 8 September 2026.
**Status:** the paper's result reproduces; the attribution in its abstract does not.
**Nothing here is an arithmetic error except item 1**, which is a typo. The substance
is item 2.

---

## The claim being checked

For a stochastic pendulum clock written as a weakly nonlinear oscillator,

    theta_ddot = -theta + eps h(theta, theta_dot) + sqrt(2D) xi,
    h = -theta_dot + f(theta, theta_dot),        D = eps T,

with `f` the escapement, the paper derives a long-time uncertainty product

    Q_inf  ~=  eps^2 / 2                                              (its Eq. 4)

and describes it as depending "solely on the degree of nonlinearity". The
derivation is two lines (its Eqs. 18 and 21):

    D_Phi   ~= eps T / (2 R*^2)          phase diffusion
    sigma   ~= eps R*^2 / (2 T)          entropy production rate
    Q_inf   =  2 D_Phi sigma / Omega^2   ->  eps^2/2,     Omega = -1 + O(eps)

`R*` is the deterministic limit-cycle radius, and it cancels. That cancellation is
the universality: the escapement enters the answer only through `R*`, so a
grasshopper, a Graham deadbeat and a van der Pol oscillator all land on the same
number. This part is right and it is lovely.

## 1. Table 1, van der Pol row: the printed phase-average kills the limit cycle

For `f = (2 - theta^2) theta_dot`,

    fbar_R(R) = (1/2pi) int_0^{2pi} sin(Phi) f(R cos Phi, R sin Phi) dPhi
              = R - R^3/8 ,

and Table 1 prints `R^3/8`. The `R` is missing.

`R* = 2` comes out of Eq. 19 either way, because of the 2 in `(2 - theta^2)`, so
nothing downstream in the paper depends on it. But the printed expression has the
wrong stability:

    d/dR ( -R/2 + R^3/8 )      at R = 2  =  +1     (unstable — no clock)
    d/dR ( -R/2 + R - R^3/8 )  at R = 2  =  -1     (stable — a clock)

Taken literally, the row describes a fixed point the oscillator runs away from.

*(A second domain note, not an error: Table 1's grasshopper average `2 theta_r/(pi R)`
holds only for `R > theta_r`. Below the critical angle the escapement never reverses,
`f = sgn(omega)`, and the average is the constant `2/pi`. The paper never needs the
other branch, because `eps` cancels out of its `R*` and `R* = 0.714 > theta_r = 0.4`
always. Anyone generalising the model does need it; using the printed branch below
`theta_r` gives `R*` about 12% too large.)*

## 2. Both surviving `eps` are the damping. The escapement's strength is not in the law.

Follow each one through the paper's own equations.

* In `D_Phi ~= eps T/(2 R*^2)`, the `eps` arrives as `D = eps T`. That is
  fluctuation–dissipation applied to the **friction coefficient**.
* In `sigma ~= eps R*^2/(2T)`, the `eps` is the friction term `-eps theta_dot`
  dissipating. The paper's Eq. 13 says so outright: `-J_Q/T = (2 eps/T)(T_eff/2 - T/2)`.

Neither is the escapement. The escapement enters only through `R*`, which cancels.

This is invisible from inside the paper because a single `eps` multiplies the whole of
`h = -theta_dot + f`, so it is at once the damping rate, the noise strength, and the
prefactor on the nonlinearity. When the answer is `eps^2`, nothing tells you which of
the three jobs earned it.

Give the damping its own symbol and keep FDT honest:

    theta_ddot = -theta - gamma theta_dot + eps f(theta, theta_dot) + sqrt(2 gamma T) xi

and the same two lines give `D_Phi = gamma T/(2R*^2)`, `sigma = gamma R*^2/(2T)`, hence

    **Q_inf ~= gamma^2 / 2  =  1 / (2 Q_qual^2)**,

with `Q_qual = omega_0/gamma` the ordinary mechanical quality factor. The paper's
`eps^2/2` is the `gamma = eps` case of this.

### The numbers

Simulator: Euler–Maruyama on the underdamped Langevin, 2000 trajectories per point,
`T = 0.005`. The current is the unwrapped `atan2(omega, theta)` — Stratonovich obeys
the chain rule, so the line integral of `grad(Phi) o x_dot` *is* the unwrapped phase.

**Gate the instrument first.** `limit_cycle_radius` returns all three of Table 1's
`R*` to 15 digits (2, 0.7136496, 1.2668787). Noiseless runs settle onto `R*` to 4–5
digits. At `gamma = eps` the simulator returns the paper's own law across three
escapements and two decades in `Q_inf`:

| model | ε | R* | Q_∞ | Q/(ε²/2) |
|---|---|---|---|---|
| van der Pol | 0.02 → 0.20 | 2.000 | 2.00e-4 → 2.05e-2 | 1.000, 1.029, 1.034, 1.027 |
| grasshopper | 0.02 → 0.20 | 0.714 | 2.13e-4 → 1.66e-2 | 1.063, 1.009, 1.004, 0.828 |
| Graham | 0.02 → 0.20 | 1.267 | 2.00e-4 → 2.07e-2 | 0.998, 0.998, 0.995, 1.036 |

(The one visible miss, grasshopper at ε = 0.2, is the breakdown of time-scale
separation the paper itself reports for large ε.)

**Then separate the two parameters.** Grasshopper escapement, `theta_r = 0.4`:

| | swept 16× | R* range | Q_∞ range | `gamma^2/2` predicts | `eps^2/2` predicts |
|---|---|---|---|---|---|
| **A** `gamma = 0.05` fixed | ε 0.0125 → 0.2 | 0.32 → 1.43 | ×2.6 | ×1 | ×256 |
| **B** `eps = 0.05` fixed | γ 0.0125 → 0.2 | 1.43 → 0.32 | ×96 | ×256 | ×1 |

At the low end of sweep A, `eps^2/2` is wrong by a factor of **17**; at the high end
of sweep B it is wrong by a factor of **17** the other way. `gamma^2/2` is within 10%
across the whole of sweep B's upper half and the whole of sweep A's lower half.

**The residual is not the escapement either.** `Q_inf/(gamma^2/2)` is a function of
`eps/gamma` alone: sweeps A and B, which reach each value of `eps/gamma` through
completely different `(gamma, eps)` pairs, agree point for point —

    eps/gamma:   0.25   0.5    1.0    2.0    4.0
    sweep A:     1.084  1.059  1.042  1.635  2.720
    sweep B:     1.039  1.032  1.098  1.715  2.776

and the excess sits **entirely in `D_Phi`** (`D/D_pred` runs 1.08 → 2.95) while
`sigma_irr` matches `gamma R*^2/(2T)` to within 1–8% at every one of the ten points.
That fingers exactly one thing: the paper's stated step of replacing the fluctuating
radius `R` by `R*` inside `Phi_dot = -1 + eps fbar_Phi(R) + noise`, which throws away
a slow frequency wander driven by amplitude fluctuations.

**A prediction, made before running it, and checked.** If that is the mechanism, the
excess must vanish for an escapement whose `fbar_Phi` is identically zero, because
then the phase velocity cannot see the amplitude. Van der Pol is that escapement:
`fbar_Phi = (1/2pi) int cos(Phi) (a - R^2 cos^2 Phi)(R sin Phi)/R dPhi = 0` by parity.
Over the *same* `eps/gamma` range where the grasshopper climbed to 2.72:

    eps/gamma:            0.75   1.0    1.5    2.0    3.0    4.0
    Q_inf/(gamma^2/2):    1.034  1.002  0.982  1.033  1.015  1.105
    Omega:               -1.0000 -0.9999 -0.9995 -0.9988 -0.9964 -0.9928

Flat. **For an escapement with no amplitude–phase coupling, `Q_inf = gamma^2/2`
independent of the escapement's strength over a factor of 5.3.**

So the honest statement is: `Q_inf = (gamma^2/2) * h(eps/gamma)`, with `h = 1` to
within a few percent wherever the paper's own `R -> R*` approximation holds, and
`h == 1` identically for `fbar_Phi = 0`. **The strength and the form of the
nonlinearity are both absent from the law. The damping is the law.**

## 3. Why the damping reading is the better paper

It lands the bound on the number horology has ranked oscillators by since the
eighteenth century. Balance wheels sit near `Q ~ 3e2`, pendulums `1e3–1e5`, quartz
`1e5–1e6`, and that is the ordering of their accuracy. The paper's stated design
principle — make `R*` large — is right for keeping the radius deterministic (small CV)
and is *not* a lever on the uncertainty product, because `R*` cancels. The lever is `Q`.

**And the bound has a lot of room in it.** Put the Shortt–Synchronome free pendulum
clock of 1921 into Eq. 18: invar bob of 6.4 kg, seconds pendulum, `Q = 110,000` in its
40 mbar tank, room temperature. The thermal floor comes out at **1.5e-10 s/day**.
Pierre Boucheron measured that clock against an atomic clock at the US Naval
Observatory in 1984 with optical sensors and got **2e-4 s/day** — a factor of **1.3
million** above the limit. Sweeping the two inputs I had to assume (swing amplitude
1–5°, pendulum length 0.25–2 m) moves the ratio only between 8.7e5 and 4.3e6.

There is a second check that needs no arithmetic. Shortt's vacuum raised `Q` from
25,000 to 110,000 (×4.4) and the reported gain in accuracy was ×4. Thermal phase
diffusion scales as `1/sqrt(Q)` and would have given ×2.1; escapement error scales as
`1/Q` and gives ×4.4. The clock was in the escapement-limited regime, by its own
measured response to being evacuated.

## What I could not settle

* Whether there is an argument that the weakly-nonlinear scaling *requires*
  `gamma ~ eps`, making the separation illegitimate rather than merely
  un-attempted. I do not think there is — the two-parameter family has a limit cycle
  and satisfies FDT for any `eps > gamma theta_r pi/4` — but I am not the author.
* The functional form of `h(eps/gamma)`. I have ten points and a mechanism, not a
  derivation. My scaling estimate for the amplitude-feedthrough term gets the
  collapse variable right and the exponent wrong, so I am not quoting it.
* `sigma_pump = eps <df/domega>` is dropped in Eq. 21 under `T_eff >> T`. Measured, it
  is at most 4% of `sigma_irr` here and falls to nothing as `R*` grows past `theta_r`,
  so it is not the residual. Whether it competes in some other regime, I do not know.
* Whether anyone has said this already. Absence of a citation is not absence of a
  predecessor, and this is the claim in here most likely to be false.

## Reproduce

```
cd reproduce/pendulum_clock_which_epsilon
python3 reproduce_paper.py       # gate: does the simulator return the paper's own eps^2/2?
python3 separate_gamma_eps.py    # sweeps A and B
python3 vdp_control.py           # the fbar_Phi = 0 control
python3 plot_separation.py       # -> which_epsilon.png
python3 shortt_thermal_floor.py  # the 1.3e6, with its assumptions swept
```

Letter sent to the author at `izumida@k.u-tokyo.ac.jp` (printed on the preprint) on
8 September 2026. If he shows me wrong, the correction goes here with attribution.
