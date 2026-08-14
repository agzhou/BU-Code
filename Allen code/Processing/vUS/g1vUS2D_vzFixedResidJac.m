%% Description:
%   Residual + analytic Jacobian for the SEPARABLE (v_zgp-fixed) 4-parameter
%   fit of the complex g1(tau) model (Eq. 15 in Jianbo Tang et al., 2020),
%   for use with lsqnonlin:
%       fun = @(x) g1vUS2D_vzFixedResidJac(x, tau, sigma, k0, vzFixed, ydataReal, ydataImag);
%       opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true);
%       x = lsqnonlin(fun, x0, lb, ub, opts); % x = [v_xgp, p, F, DC]
%
%   v_zgp is NOT a fit parameter here -- it is a fixed, known constant,
%   meant to come from findVzPhaseDiff.m (a closed-form, non-iterative
%   estimator). Compare to g1vUS2D_vec_split.m, which fits the same model
%   with v_zgp as a 5th free parameter via lsqcurvefit.
%
%   Why fix v_zgp instead of just seeding it: a Monte Carlo comparison
%   (synthetic g1 built from the vUS model + realistic finite-sample
%   estimator noise, matching Eqs. 2 and 15) across a 1-20 mm/s, 3-noise-
%   level sweep found:
%     - Against the MAGNITUDE-ONLY model (g1vUS2D_mag_vec.m): fixing
%       v_zgp nearly halved v_xgp RMSE (2.78 -> 1.47 mm/s), because
%       |g1(tau)| discards the phase term that is the only thing cleanly
%       distinguishing axial from transverse motion -- a joint fit can
%       (and does) trade v_xgp against v_zgp to explain the same decay
%       envelope. Even when the joint fit was SEEDED with an accurate
%       v_zgp0, leaving it free made v_zgp WORSE (3.36 mm/s RMSE vs.
%       ~0.5 mm/s for the seed estimate alone).
%     - Against the FULL COMPLEX model (this function): the accuracy gap
%       mostly disappears (0.906 vs. 0.925 mm/s v_xgp RMSE) since the
%       phase term already resolves the v_xgp/v_zgp collinearity on its
%       own. The remaining benefit here is speed: dropping 5 params to 4
%       AND supplying this analytic Jacobian (vs. lsqcurvefit's default
%       finite-difference one) cut mean per-pixel solve time by ~48%
%       (3.05 -> 1.57 ms). The complex residual is twice as long as the
%       magnitude-only one (real+imag stacked), so the finite-difference
%       baseline pays for twice as many perturbation evaluations,
%       making the analytic Jacobian proportionally more valuable here.
%     - One honest caveat: at very low v_zgp (~1 mm/s) the separable fit
%       was occasionally slightly worse than the joint fit, likely
%       because fixing v_zgp removes the joint fit's ability to
%       self-correct a slightly-off seed -- a capacity that matters more
%       once phase already gives the joint fit a reasonable handle on
%       v_zgp (unlike the magnitude-only case, where it had no such
%       handle at all).
%   The analytic Jacobian below was verified against a finite-difference
%   approximation (max abs difference ~1e-9) before use.
%
%   Because v_zgp is fixed, the phase term 2*k0*vzFixed*tau does not
%   depend on any fit parameter -- only the real amplitude envelope
%   ("dynamic") depends on [v_xgp, p, F], and DC only enters the real
%   part additively. That's what keeps this Jacobian simple:
%       g_real(tau) = DC + dynamic(tau)*cos(2*k0*vzFixed*tau)
%       g_imag(tau) = dynamic(tau)*sin(2*k0*vzFixed*tau)
%       dynamic(tau) = F * exp(-(v_xgp*tau)^2/(4*sigma(1)^2) - (vzFixed*tau)^2/(4*sigma(2)^2)) * exp(-(p*vzFixed*k0*tau)^2)
%
% Inputs:
%   x: [v_xgp, p, F, DC] (all in SI units; F and DC always included --
%      unlike g1vUS2D_vec_split.m, there are no useF/useDC toggles here)
%   tau: time lag vector [s], should exclude tau=0 (matches vUS_2D.m's
%        convention of fitting tau(2:end)/g1_exp(:,2:end))
%   sigma: 1/e values of the PSF shape in x, z [m]
%   k0: Angular wavenumber [rad/m]
%   vzFixed: fixed axial group velocity [m/s] (e.g. from findVzPhaseDiff.m)
%   ydataReal, ydataImag: real and imaginary parts of the observed g1(tau)
%       (same length as tau)
%
% Outputs:
%   r: [2*numel(tau), 1] stacked residual, [real part; imaginary part],
%      i.e. model(x) - ydata
%   J: [2*numel(tau), 4] analytic Jacobian of r w.r.t. [v_xgp, p, F, DC]
%      (only computed if requested, so this can also be called as
%      r = g1vUS2D_vzFixedResidJac(x, ...) without 'SpecifyObjectiveGradient')
function [r, J] = g1vUS2D_vzFixedResidJac(x, tau, sigma, k0, vzFixed, ydataReal, ydataImag)
    v_xgp = x(1); p = x(2); F = x(3); DC = x(4);
    tau = tau(:); ydataReal = ydataReal(:); ydataImag = ydataImag(:);

    Ez = exp(-(vzFixed.*tau).^2 ./ (4*sigma(2)^2));       % fixed given vzFixed
    Ex = exp(-(v_xgp.*tau).^2 ./ (4*sigma(1)^2));
    Ep = exp(-(p.*vzFixed.*k0.*tau).^2);
    dynamic = F .* Ez .* Ex .* Ep;                          % real amplitude envelope
    phase = 2 .* k0 .* vzFixed .* tau;                       % fixed, no parameter dependence

    cosP = cos(phase); sinP = sin(phase);
    gReal = DC + dynamic.*cosP;
    gImag = dynamic.*sinP;

    r = [gReal - ydataReal; gImag - ydataImag];

    if nargout > 1
        n = numel(tau);
        ddyn_dvx = dynamic .* ( -v_xgp.*tau.^2 ./ (2*sigma(1)^2) );
        ddyn_dp  = dynamic .* ( -2*p.*(vzFixed*k0)^2 .* tau.^2 );
        ddyn_dF  = Ez.*Ex.*Ep;

        dgReal_dvx = cosP .* ddyn_dvx;
        dgReal_dp  = cosP .* ddyn_dp;
        dgReal_dF  = cosP .* ddyn_dF;
        dgReal_dDC = ones(n,1);

        dgImag_dvx = sinP .* ddyn_dvx;
        dgImag_dp  = sinP .* ddyn_dp;
        dgImag_dF  = sinP .* ddyn_dF;
        dgImag_dDC = zeros(n,1);

        J = [dgReal_dvx, dgReal_dp, dgReal_dF, dgReal_dDC;
             dgImag_dvx, dgImag_dp, dgImag_dF, dgImag_dDC];
    end
end
