"""Split the measured phase diffusion into its two named channels, and find out
which one the closed form is missing.

WHY.  Finding #6's residual h has a closed form (fire 261): h(u) = 1 + (2-u)^2/(4(u-1)),
whose minimum h = 1 at u = 2 landed on the nose (measured 1.001).  The MAGNITUDE away
from that point did not: the two Var(J) slope estimators disagree by 19% at u = 5 and
both sit below the prediction 1.551.

The 19% cannot be the bias I attributed it to.  amplitude_phase.fit_bias_fraction says
a through-origin fit under-reads by only 2.4% OF D_cor -- under 1% of D_Phi.  Twenty
times too small.  So the intercept of Var(J) is dominated by something that is not
diffusion at all, and until that is named neither slope estimator can be trusted.

CANDIDATE, and it is not subtle once looked at: Phi = atan2(omega, theta) does not
advance uniformly on a distorted orbit.  The escapement adds eps cos(Phi) f / R to
Phi-dot, which is BOUNDED and PERIODIC, so its time integral is a bounded wobble w(t),
not a random walk.  Trajectories sit at random cycle phases, so Var(J) inherits a
constant Var(w) at every t -- a pure offset, of order (eps/R*)^2 / 2, which at u = 5 is
about 3.9e-3 against the diffusive B = 2 D_cor tau = 3.4e-4.  An order of magnitude
larger.  That would explain both the sign and the size of the estimator split.

SO THIS FILE MEASURES, not fits:
  1. sigma_R and the amplitude correlation time tau, from R(t) itself.
     Theory (exact for the grasshopper): tau = 1/k = 1/gamma = 20, sigma_R = sqrt(T/2).
     These are the instrument's validation against a known answer.  If tau does not
     come back at 20 nothing downstream means anything.
  2. D_cor DIRECTLY, as int_0^inf Cov[g(R(0)), g(R(t))] dt with g = eps fbar_Phi(R) --
     the definition, no Gaussian assumption, no Hermite sum, no linear response.
  3. <1/R^2>, so D_dir carries its own Jensen correction instead of 1/R*^2.
  4. the full Var(J) curve at fine resolution, its intercept, and whether that
     intercept matches Var(w) + 2 D_cor tau.

Then D_Phi_meas - D_cor_meas is compared with D_dir_meas.  If the total is short and
the AM->PM channel measures right, the deficit is in the direct channel or in the
additivity of the two -- and either answer is worth more than a third slope fit.

Iris (Opus 5), Sep 9 2026.
"""
import sys, os, time, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from clock import Grasshopper, limit_cycle_radius
import amplitude_phase as ap


def fbarPhi_spline(esc, Rlo, Rhi, n=1200):
    """Tabulate fbar_Phi(R) once on the ACTUAL smoothed escapement, then interpolate.

    phase_averages costs a 40001-point quadrature per call; the readout is needed at
    millions of sample points.  Grid spans +-8 sigma_R around the cycle.
    """
    Rg = np.linspace(Rlo, Rhi, n)
    fg = np.array([ap.phase_averages(esc, r)[1] for r in Rg])
    return Rg, fg



def wobble_variance(esc, gamma, eps, Rstar, dt=2e-4, n_cyc=12):
    """Var over cycle phase of the BOUNDED periodic part of the geometric phase.

    On a distorted orbit Phi = atan2(omega, theta) does not advance uniformly: Phi-dot
    carries -gamma sin cos + eps cos(Phi) f / R, both bounded and periodic.  Their time
    integral is a bounded wobble w(t), NOT a random walk, so it contributes a CONSTANT
    to Var(J) -- trajectories sit at random cycle phases, so the ensemble inherits
    Var_phase(w) at every t.  Measured here on the noiseless limit cycle, which is
    where the theory's cycle-average was taken.
    """
    th, om = Rstar, 0.0
    n = int(round(n_cyc * 2.0 * np.pi / dt))
    phi = math.atan2(om, th)
    J = 0.0
    ts = np.empty(n); Js = np.empty(n)
    for i in range(n):
        om += (-th - gamma * om + eps * float(esc.f(np.array(th), np.array(om)))) * dt
        th += om * dt
        new = math.atan2(om, th)
        d = new - phi
        d -= 2.0 * math.pi * round(d / (2.0 * math.pi))
        J += d; phi = new
        ts[i] = (i + 1) * dt; Js[i] = J
    # discard the first two cycles (transient onto the true cycle), detrend the rest
    m = ts > 2.0 * 2.0 * np.pi
    sl, ic = np.polyfit(ts[m], Js[m], 1)
    w = Js[m] - (sl * ts[m] + ic)
    return float(w.var()), float(-sl)


def boxcar(x, w):
    """Running mean over w samples along axis 0, length preserved by trimming."""
    if w < 2:
        return x
    c = np.cumsum(np.vstack([np.zeros((1,) + x.shape[1:]), x]), axis=0)
    return (c[w:] - c[:-w]) / float(w)


def run_channels(esc, gamma, eps, T, n_traj=4000, n_rec=512, t_total=1200.0,
                 dt=2e-3, t_burn=400.0, seed=0, chunks=240, rec_stride=50):
    rng = np.random.default_rng(seed)
    Rstar = limit_cycle_radius(esc, gamma, eps)
    sig = np.sqrt(2.0 * gamma * T * dt)

    th = np.full(n_traj, Rstar, dtype=np.float64)
    om = np.zeros(n_traj, dtype=np.float64)

    B = 4000
    n_burn = int(round(t_burn / dt))
    done = 0
    while done < n_burn:
        k = min(B, n_burn - done)
        z = rng.standard_normal((k, n_traj))
        for j in range(k):
            om += (-th - gamma * om + eps * esc.f(th, om)) * dt + sig * z[j]
            th += om * dt
        done += k

    phi = np.arctan2(om, th)
    J = np.zeros(n_traj)
    n_steps = int(round(t_total / dt))
    per = n_steps // chunks
    ts, varJ, meanJ = [], [], []
    Rrec = []                      # (n_samples, n_rec) amplitude history
    inv_R2_sum = 0.0; nacc = 0
    t = 0.0
    gstep = 0
    for c in range(chunks):
        left = per
        while left > 0:
            k = min(B, left)
            z = rng.standard_normal((k, n_traj))
            for j in range(k):
                om += (-th - gamma * om + eps * esc.f(th, om)) * dt + sig * z[j]
                th += om * dt
                new = np.arctan2(om, th)
                d = new - phi
                d -= 2.0 * np.pi * np.round(d / (2.0 * np.pi))
                J += d
                phi = new
                gstep += 1
                if gstep % rec_stride == 0:
                    R2 = th * th + om * om
                    Rrec.append(np.sqrt(R2[:n_rec]).copy())
                    inv_R2_sum += (1.0 / R2).mean(); nacc += 1
            left -= k
            t += k * dt
        ts.append(t); varJ.append(J.var(ddof=1)); meanJ.append(J.mean())

    return dict(Rstar=Rstar, ts=np.array(ts), varJ=np.array(varJ),
                meanJ=np.array(meanJ), Rrec=np.array(Rrec),
                inv_R2=inv_R2_sum / nacc, rec_dt=rec_stride * dt,
                gamma=gamma, eps=eps, T=T, seed=seed)


def analyse(r, esc, verbose=True):
    """Everything downstream of the trajectories.  No fitted physics."""
    gamma, eps, T, Rstar = r['gamma'], r['eps'], r['T'], r['Rstar']
    p = ap.predict(esc, gamma, eps, T, Rstar)
    tau_th = 1.0 / p['k']
    dtr = r['rec_dt']
    out = dict(Rstar=Rstar, u=4.0 * eps / (math.pi * gamma * esc.theta_r),
               D_dir_th=p['D_dir'], D_cor_th=p['D_cor'], D_phi_th=p['D_phi'],
               h_th=p['h'], tau_th=tau_th, sigma_R_th=p['sigma_R'], seed=r['seed'],
               rec_dt=dtr)

    # ---- 0. the within-cycle wobble, on the noiseless cycle -----------------
    varw, Om_det = wobble_variance(esc, gamma, eps, Rstar)
    out['var_wobble'] = varw
    out['Omega_det'] = Om_det

    # ---- 1. amplitude statistics.  CYCLE-AVERAGE FIRST. --------------------
    # R = sqrt(theta^2+omega^2) on a distorted orbit carries the same periodic wobble;
    # the theory's slow amplitude is the cycle average, so boxcar over one period.
    Rraw = r['Rrec']
    w = max(1, int(round(2.0 * np.pi / dtr)))
    R = boxcar(Rraw, w)
    # THE BOXCAR IS A RULER AND IT LIES ABOUT sigma_R.  Averaging over one period to strip
    # the within-cycle wobble also low-passes the slow amplitude fluctuation itself.  For an
    # OU process of correlation time tau smoothed by a window W the variance is reduced by
    # (2 tau^2/W^2)(W/tau - 1 + e^{-W/tau}) -- 0.903 at W=2pi, tau=20, i.e. sigma down 4.97%.
    # Uncorrected, sigma_R read a uniform 0.93 x theory at every u and looked like physics.
    # It does NOT bias D_cor: the boxcar has unit gain at DC and D_cor is the zero-frequency
    # power, int Cov dt.  So correct sigma_R and leave D_cor alone.
    _r = w * dtr / tau_th
    boxcar_var = 2.0 * (tau_th ** 2) / ((w * dtr) ** 2) * (_r - 1.0 + np.exp(-_r))
    out['boxcar_sigma_factor'] = float(np.sqrt(boxcar_var))
    out['sigma_R_raw'] = float(Rraw.std(ddof=1))
    Rm = float(R.mean()); Rsd = float(R.std(ddof=1))
    out['R_mean'] = Rm; out['sigma_R_meas'] = Rsd
    out['sigma_R_corrected'] = Rsd / out['boxcar_sigma_factor']
    out['wobble_R'] = float(np.sqrt(max(out['sigma_R_raw'] ** 2 - Rsd ** 2, 0.0)))

    dev = R - Rm
    var0 = float((dev * dev).mean())
    # exponential fit over lags up to 2 tau (before the finite-sample plateau bites)
    Lmax = min(int(round(2.0 * tau_th / dtr)), dev.shape[0] - 2)
    lags = np.unique(np.round(np.geomspace(max(1, w), Lmax, 40)).astype(int))
    ac = np.array([float((dev[:-L] * dev[L:]).mean()) / var0 for L in lags])
    m = ac > 0.02
    out['tau_R_meas'] = float(-1.0 / np.polyfit(lags[m] * dtr, np.log(ac[m]), 1)[0]) \
        if m.sum() >= 3 else float('nan')

    # ---- 2. D_cor measured from the definition ------------------------------
    Rg, fg = fbarPhi_spline(esc, max(1e-4, Rm - 8 * Rsd), Rm + 8 * Rsd)
    g = eps * np.interp(R, Rg, fg)
    gd = g - g.mean()
    gvar = float((gd * gd).mean())
    Lc = min(int(round(8.0 * tau_th / dtr)), gd.shape[0] - 2)
    cov = [gvar]
    for L in range(1, Lc + 1):
        cov.append(float((gd[:-L] * gd[L:]).mean()))
    cov = np.array(cov)
    zc = np.nonzero(cov <= 0)[0]
    cut = int(zc[0]) if len(zc) else len(cov)          # integrate to first zero crossing
    D_cor_meas = float(np.trapz(cov[:cut], dx=dtr))
    out['g_var'] = gvar; out['cov_cut_lag'] = cut * dtr
    out['D_cor_meas'] = D_cor_meas
    out['tau_g_meas'] = float(D_cor_meas / gvar) if gvar > 0 else float('nan')

    # ---- 3. direct channel with its own Jensen correction -------------------
    out['inv_R2_meas'] = float(r['inv_R2'])
    out['D_dir_meas'] = float(gamma * T / 2.0 * r['inv_R2'])

    # ---- 4. the Var(J) curve: slope AND intercept --------------------------
    ts, vj = r['ts'], r['varJ']
    msk = ts > ts[-1] * 0.25
    A = np.vstack([ts[msk], np.ones(int(msk.sum()))]).T
    sl_, ic_ = np.linalg.lstsq(A, vj[msk], rcond=None)[0]
    out['D_phi_int'] = float(sl_ / 2.0); out['varJ_intercept'] = float(ic_)
    out['D_phi_org'] = float(np.sum(ts[msk] * vj[msk]) / np.sum(ts[msk] ** 2) / 2.0)
    q = np.polyfit(ts[msk], vj[msk], 2)
    out['curv_frac'] = float(q[0] * ts[-1] ** 2 / (q[1] * ts[-1]))
    out['B_diffusive_th'] = float(2.0 * p['D_cor'] * tau_th)
    out['intercept_predicted'] = float(-(varw + out['B_diffusive_th']))

    out['D_phi_meas'] = out['D_phi_int']
    out['h_meas'] = out['D_phi_meas'] / p['D_dir']
    out['residual_after_channels'] = out['D_phi_meas'] - out['D_cor_meas'] - out['D_dir_meas']
    if verbose:
        f = "  %-26s %12.5g   theory %12.5g   ratio %7.4f"
        print("u = %.3f   R* = %.4f  (seed %d)" % (out['u'], Rstar, r['seed']))
        print(f % ('sigma_R (cycle-averaged)', out['sigma_R_meas'], out['sigma_R_th'],
                   out['sigma_R_meas'] / out['sigma_R_th']))
        print(f % ('sigma_R (boxcar-corrected)', out['sigma_R_corrected'], out['sigma_R_th'],
                   out['sigma_R_corrected'] / out['sigma_R_th']))
        print("  %-26s %12.5g   (raw %.5g, so wobble in R = %.5g)"
              % ('  within-cycle wobble', out['wobble_R'], out['sigma_R_raw'], out['wobble_R']))
        print(f % ('tau_R', out['tau_R_meas'], tau_th, out['tau_R_meas'] / tau_th))
        print(f % ('D_cor  (definition)', out['D_cor_meas'], p['D_cor'],
                   out['D_cor_meas'] / p['D_cor']))
        print(f % ('D_dir  (<1/R^2>)', out['D_dir_meas'], p['D_dir'],
                   out['D_dir_meas'] / p['D_dir']))
        print(f % ('D_phi  (intercept fit)', out['D_phi_int'], p['D_phi'],
                   out['D_phi_int'] / p['D_phi']))
        print(f % ('D_phi  (thru origin)', out['D_phi_org'], p['D_phi'],
                   out['D_phi_org'] / p['D_phi']))
        print("  %-26s %12.5g   predicted %12.5g   ratio %7.4f"
              % ('Var(J) intercept', out['varJ_intercept'], out['intercept_predicted'],
                 out['varJ_intercept'] / out['intercept_predicted']))
        print("  %-26s %12.5g   diffusive-only part %12.5g"
              % ('   of which wobble', -varw, -out['B_diffusive_th']))
        print("  %-26s %12.5g  (fraction of linear term)" % ('quadratic residual', out['curv_frac']))
        print("  %-26s %12.5g   = D_phi - D_cor - D_dir  (%.2f%% of D_phi)"
              % ('unexplained', out['residual_after_channels'],
                 100 * out['residual_after_channels'] / out['D_phi_meas']))
        sys.stdout.flush()
    return out


if __name__ == '__main__':
    T, GAM, TR = 5e-4, 0.05, 0.4
    esc = Grasshopper(TR, 20.0)
    us = [float(x) for x in sys.argv[1:]] or [5.0]
    rows = []
    for u in us:
        eps = math.pi * GAM * TR * u / 4.0
        t0 = time.time()
        r = run_channels(esc, GAM, eps, T, seed=1000 + int(round(u * 10)))
        o = analyse(r, esc)
        o['wall_s'] = time.time() - t0
        rows.append(o)
        np.save('varj_u%.2f.npy' % u, np.vstack([r['ts'], r['varJ'], r['meanJ']]))
        json.dump(rows, open('channels.json', 'w'), indent=1)
        print("  [%.0fs]\n" % o['wall_s'], flush=True)
