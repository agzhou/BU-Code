%% Description:
%   Residual + analytic Jacobian for the FULL 5-parameter joint fit of the
%   complex g1(tau) model (Eq. 15 in Jianbo Tang et al., 2020), for use
%   with lsqnonlin:
%       fun = @(x) g1vUS2D_residJac(x, tau, sigma, k0, ydataReal, ydataImag);
%       opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true);
%       x = lsqnonlin(fun, x0, lb, ub, opts); % x = [v_xgp, v_zgp, p, F, DC]
%
%   Fits the SAME model as g1vUS2D_vec_split.m (which remains the more
%   readable reference and is what lsqcurvefit needs), just packaged as a
%   residual+Jacobian pair for lsqnonlin. v_zgp is free here, unlike
%   g1vUS2D_vzFixedResidJac.m's separable (v_zgp-fixed) 4-parameter
%   version -- use this when you want the standard joint fit, or that one
%   when you have an independent v_zgp estimate (e.g. from
%   findVzPhaseDiff.m) and want to fix it.
%
%   Analytic Jacobian is more involved here than in the v_zgp-fixed case:
%   v_zgp now appears in the phase term too (not just the amplitude
%   envelope), so the real/imaginary parts each pick up an extra term from
%   the chain rule through cos/sin(2*k0*v_zgp*tau). Verified against a
%   finite-difference approximation (max abs difference ~4e-7) before use.
%
%   IMPORTANT CAVEAT, carried over from g1vUS2D_vzFixedResidJac.m: this
%   function does NOT resolve the v_xgp/p identifiability problem found in
%   this model. v_xgp and p only ever enter through the same combined
%   decay coefficient (v_xgp^2/(4*sigma(1)^2) + p^2*(v_zgp*k0)^2), so a
%   SINGLE g1(tau) curve (this function's input) cannot separate them --
%   confirmed by a noise-free test where v_zgp/F/DC recovered exactly but
%   v_xgp and p landed on a different point of the same degenerate ridge
%   (true [v_xgp,p]=[5.0,0.30] fit to [3.12,0.31]) while the combined
%   coefficient matched exactly (14968.39 both times). This function is
%   purely a faster way to run the SAME (still degenerate) joint fit --
%   supplying this analytic Jacobian to lsqnonlin should give a speed
%   benefit over lsqcurvefit's default finite-difference Jacobian, the
%   same way it did for the separable 4-parameter case (~18-48% faster
%   per pixel there), but that has NOT been separately benchmarked for
%   this 5-parameter case.
%
% Inputs:
%   x: [v_xgp, v_zgp, p, F, DC] (all SI units; F and DC always included --
%      no useF/useDC toggles, matching g1vUS2D_vzFixedResidJac.m's style)
%   tau: time lag vector [s], should exclude tau=0 (matches vUS_2D.m's
%        convention of fitting tau(2:end)/g1_exp(:,2:end))
%   sigma: 1/e values of the PSF shape in x, z [m]
%   k0: Angular wavenumber [rad/m]
%   ydataReal, ydataImag: real and imaginary parts of the observed g1(tau)
%       (same length as tau)
%
% Outputs:
%   r: [2*numel(tau), 1] stacked residual, [real part; imaginary part]
%   J: [2*numel(tau), 5] analytic Jacobian of r w.r.t.
%      [v_xgp, v_zgp, p, F, DC] (only computed if requested)
function [r, J] = g1vUS2D_residJac(x, tau, sigma, k0, ydataReal, ydataImag)
    v_xgp = x(1); v_zgp = x(2); p = x(3); F = x(4); DC = x(5);
    tau = tau(:); ydataReal = ydataReal(:); ydataImag = ydataImag(:);

    Ex = exp(-(v_xgp.*tau).^2 ./ (4*sigma(1)^2));
    Ez = exp(-(v_zgp.*tau).^2 ./ (4*sigma(2)^2));
    Ep = exp(-(p.*v_zgp.*k0.*tau).^2);
    dynamic = F .* Ex .* Ez .* Ep;                 % real amplitude envelope
    phase = 2 .* k0 .* v_zgp .* tau;                % now depends on v_zgp

    cosP = cos(phase); sinP = sin(phase);
    gReal = DC + dynamic.*cosP;
    gImag = dynamic.*sinP;

    r = [gReal - ydataReal; gImag - ydataImag];

    if nargout > 1
        n = numel(tau);
        ddyn_dvx = dynamic .* ( -v_xgp.*tau.^2 ./ (2*sigma(1)^2) );
        ddyn_dvz = dynamic .* ( -v_zgp.*tau.^2./(2*sigma(2)^2) - 2*p^2*k0^2.*v_zgp.*tau.^2 );
        ddyn_dp  = dynamic .* ( -2*p.*(v_zgp*k0)^2 .* tau.^2 );
        ddyn_dF  = Ex.*Ez.*Ep;
        dphase_dvz = 2*k0.*tau;

        dgReal_dvx = cosP .* ddyn_dvx;
        dgReal_dvz = cosP .* ddyn_dvz - dynamic .* sinP .* dphase_dvz;
        dgReal_dp  = cosP .* ddyn_dp;
        dgReal_dF  = cosP .* ddyn_dF;
        dgReal_dDC = ones(n,1);

        dgImag_dvx = sinP .* ddyn_dvx;
        dgImag_dvz = sinP .* ddyn_dvz + dynamic .* cosP .* dphase_dvz;
        dgImag_dp  = sinP .* ddyn_dp;
        dgImag_dF  = sinP .* ddyn_dF;
        dgImag_dDC = zeros(n,1);

        J = [dgReal_dvx, dgReal_dvz, dgReal_dp, dgReal_dF, dgReal_dDC;
             dgImag_dvx, dgImag_dvz, dgImag_dp, dgImag_dF, dgImag_dDC];
    end
end
