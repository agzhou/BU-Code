function [r, J] = vUS_2D_erf_windowed_residJac(x, tau, k0, sigma, f_min, f_max, ydataReal, ydataImag)
%% Description:
%   Residual + analytic Jacobian for the windowed (partial-volume) g1
%   model, for direct use with lsqnonlin:
%       fun = @(x) vUS_2D_erf_windowed_residJac(x, tau, k0, sigma, f_min, f_max, ydataReal, ydataImag);
%       opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true);
%       x = lsqnonlin(fun, x0, lb, ub, opts); % x = [v_xgp, v_zgp, F, DC]
%
% Inputs:
%   x: [v_xgp, v_zgp, F, DC]
%   tau: time lag vector [s] (exclude tau=0)
%   k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_z] [m]
%   f_min, f_max: sampled window bounds within [0, 2] (fixed, not fitted)
%   ydataReal, ydataImag: real/imag parts of the observed g1(tau)
%
% Outputs:
%   r: [2*numel(tau), 1] stacked residual, [real part; imaginary part]
%   J: [2*numel(tau), 4] analytic Jacobian of r w.r.t. [v_xgp,v_zgp,F,DC]

    tau = tau(:); ydataReal = ydataReal(:); ydataImag = ydataImag(:);

    if nargout > 1
        [g1_split, J] = vUS_2D_erf_windowed_vec_split_Jac(x, tau, k0, sigma, f_min, f_max);
    else
        g1_split = vUS_2D_erf_windowed_vec_split_Jac(x, tau, k0, sigma, f_min, f_max);
    end

    r = [g1_split(:,1) - ydataReal; g1_split(:,2) - ydataImag];
end
