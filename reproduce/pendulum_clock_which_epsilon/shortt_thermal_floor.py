"""How far above its own thermal noise floor was the best pendulum clock ever built?

The law under test (arXiv:2609.04957, read as a statement about damping) fixes the
THERMAL floor of a pendulum clock.  In physical units the phase diffusion of a
weakly damped oscillator held at amplitude A is

    D_phi = gamma kB T / (4E),      E = (1/2) I omega0^2 A^2,   gamma = omega0/Q

(this is the paper's Eq. 18, D_Phi = D/(2R*^2), with D = gamma kB T / (I omega0^2)
after nondimensionalising time by omega0 -- see notes).  Accumulated timing error
over an interval tau is a random walk: sd(t_err) = sqrt(2 D_phi tau)/omega0.

Inputs, and where each came from:
  * bob mass 6.4 kg (14 lb), invar rod, vacuum tank            -- Wikipedia,
    Shortt-Synchronome clock, sourced to Bosschieter
  * Q = 110,000 in the 30 mmHg tank (25,000 at atmosphere)     -- same
  * measured stability 200 microseconds/day = 2.31 ppb         -- Boucheron 1984,
    optical sensors against an atomic clock at USNO, one month
  * seconds pendulum (2 s period, L = 0.994 m) -- ASSUMED from the 125 cm tank and
    standard regulator practice, not stated in the source.
  * swing amplitude 1.5 deg -- ASSUMED (precision regulator practice).  Swept below.
"""
import numpy as np

kB = 1.380649e-23
Temp = 293.0          # K, an observatory room

def floor(m=6.4, L=0.994, Q=110_000.0, amp_deg=1.5, tau=86400.0, Temp=Temp):
    omega0 = np.sqrt(9.80665 / L)
    I = m * L**2
    A = np.deg2rad(amp_deg)
    E = 0.5 * I * omega0**2 * A**2
    gamma = omega0 / Q
    D_phi = gamma * kB * Temp / (4 * E)
    sd_t = np.sqrt(2 * D_phi * tau) / omega0
    return dict(omega0=omega0, period=2*np.pi/omega0, I=I, E=E, gamma=gamma,
                D_phi=D_phi, sd_t=sd_t, frac=sd_t/tau)

MEASURED = 200e-6      # s per day, Boucheron 1984

r = floor()
print("Shortt primary pendulum, thermal floor")
print("  period            %.4f s   (omega0 = %.4f rad/s)" % (r['period'], r['omega0']))
print("  oscillation energy %.4g J" % r['E'])
print("  gamma = omega0/Q   %.4g 1/s" % r['gamma'])
print("  D_phi              %.4g rad^2/s" % r['D_phi'])
print("  timing sd in 1 day %.4g s   (%.3g ns)" % (r['sd_t'], r['sd_t']*1e9))
print("  fractional         %.3g" % r['frac'])
print()
print("  MEASURED (Boucheron 1984, USNO, vs atomic clock, 1 month):")
print("    %.4g s/day  ->  %.3g ppb" % (MEASURED, MEASURED/86400*1e9))
print("  ratio measured / thermal floor:  %.3g" % (MEASURED / r['sd_t']))
print()
print("sensitivity of the floor to the two assumed inputs")
for a in (1.0, 1.5, 2.0, 3.0, 5.0):
    print("   amplitude %.1f deg  ->  floor %.3g s/day, ratio %.3g"
          % (a, floor(amp_deg=a)['sd_t'], MEASURED/floor(amp_deg=a)['sd_t']))
for LL, nm in ((0.994, "seconds pendulum"), (0.248, "half-second"), (2.0, "2 m")):
    f = floor(L=LL)
    print("   L = %.3f m (%-16s) -> floor %.3g s/day, ratio %.3g"
          % (LL, nm, f['sd_t'], MEASURED/f['sd_t']))
print()
print("the same clock at atmosphere (Q = 25,000):")
f25 = floor(Q=25_000.0)
print("   thermal floor %.3g s/day  (worse by sqrt(110/25) = %.2f)"
      % (f25['sd_t'], np.sqrt(110/25)))
print("   but Shortt reports the VACUUM improved real accuracy by a factor of ~4,")
print("   i.e. error ~ 1/Q, not 1/sqrt(Q): the escapement-error scaling, not the")
print("   thermal one.  Independent evidence the clock was nowhere near the floor.")
