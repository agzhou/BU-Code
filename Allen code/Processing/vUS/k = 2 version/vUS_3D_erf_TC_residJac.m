function [r, J] = vUS_3D_erf_TC_residJac(x, tau, k0, sigma, ydataReal, ydataImag)
%% Description:
%   TC (Transverse-Combined) analog of vUS_3D_erf_residJac.m. Residual +
%   analytic Jacobian for the new (Poiseuille / uniform-velocity-
%   distribution) g1 model with combined transverse velocity, for direct
%   use with lsqnonlin:
%       fun = @(x) vUS_3D_erf_TC_residJac(x, tau, k0, sigma, ydataReal, ydataImag);
%       opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true);
%       x = lsqnonlin(fun, x0, lb, ub, opts); % x = [v_tgp, v_zgp, F, DC]
%
%   Since the residual is model - data and data doesn't depend on x, the
%   residual's Jacobian equals the model's (from
%   vUS_3D_erf_TC_vec_split_Jac.m / vUS_3D_erf_TC_complex_Jac.m).
%
%   NOTE: v_tgp is a combined magnitude (v_tgp^2 = v_xgp^2+v_ygp^2), so
%   it is non-negative by construction -- pass a lower bound of 0 to
%   lsqnonlin, same as the sign-not-of-interest convention used for
%   v_xgp elsewhere in this codebase, except here it isn't a convention
%   choice: v_tgp has no sign to begin with.
%
% Inputs:
%   x: [v_tgp, v_zgp, F, DC] (DC real-valued)
%   tau: time lag vector [s] (exclude tau=0, matching the rest of the
%       pipeline's convention)
%   k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_y, sigma_z] [m] -- sigma_x and sigma_y must be
%       equal, matching vUS_3D_erf_TC.m
%   ydataReal, ydataImag: real/imag parts of the observed g1(tau), same
%       length as tau
%
% Outputs:
%   r: [2*numel(tau), 1] stacked residual, [real part; imaginary part]
%   J: [2*numel(tau), 4] analytic Jacobian of r w.r.t. [v_tgp, v_zgp, F, DC]

    tau = tau(:); ydataReal = ydataReal(:); ydataImag = ydataImag(:);

    if nargout > 1
        [g1_split, J] = vUS_3D_erf_TC_vec_split_Jac(x, tau, k0, sigma);
    else
        g1_split = vUS_3D_erf_TC_vec_split_Jac(x, tau, k0, sigma);
    end

    r = [g1_split(:,1) - ydataReal; g1_split(:,2) - ydataImag];
end
