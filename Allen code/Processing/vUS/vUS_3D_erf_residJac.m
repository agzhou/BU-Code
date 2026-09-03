function [r, J] = vUS_3D_erf_residJac(x, tau, k0, sigma, ydataReal, ydataImag)
%% Description:
%   3D analog of vUS_2D_erf_residJac.m. Residual + analytic Jacobian for
%   the new (Poiseuille / uniform-velocity-distribution) g1 model, for
%   direct use with lsqnonlin:
%       fun = @(x) vUS_3D_erf_residJac(x, tau, k0, sigma, ydataReal, ydataImag);
%       opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true);
%       x = lsqnonlin(fun, x0, lb, ub, opts); % x = [v_xgp, v_ygp, v_zgp, F, DC]
%
%   Since the residual is model - data and data doesn't depend on x, the
%   residual's Jacobian equals the model's (from
%   vUS_3D_erf_vec_split_Jac.m / vUS_3D_erf_complex_Jac.m).
%
% Inputs:
%   x: [v_xgp, v_ygp, v_zgp, F, DC] (DC real-valued)
%   tau: time lag vector [s] (exclude tau=0, matching the rest of the
%       pipeline's convention)
%   k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_y, sigma_z] [m]
%   ydataReal, ydataImag: real/imag parts of the observed g1(tau), same
%       length as tau
%
% Outputs:
%   r: [2*numel(tau), 1] stacked residual, [real part; imaginary part]
%   J: [2*numel(tau), 5] analytic Jacobian of r w.r.t.
%      [v_xgp, v_ygp, v_zgp, F, DC]

    tau = tau(:); ydataReal = ydataReal(:); ydataImag = ydataImag(:);

    if nargout > 1
        [g1_split, J] = vUS_3D_erf_vec_split_Jac(x, tau, k0, sigma);
    else
        g1_split = vUS_3D_erf_vec_split_Jac(x, tau, k0, sigma);
    end

    r = [g1_split(:,1) - ydataReal; g1_split(:,2) - ydataImag];
end
