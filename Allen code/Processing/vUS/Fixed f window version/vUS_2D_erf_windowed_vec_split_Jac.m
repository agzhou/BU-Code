function [g1_split, J] = vUS_2D_erf_windowed_vec_split_Jac(x, tau, k0, sigma, f_min, f_max)
%% Description:
%   Value + analytic Jacobian of the windowed (partial-volume) g1 model,
%   packaged to match vUS_2D_erf_windowed_vec_split.m's [nTau, 2]
%   real/imaginary-split output convention, for use with lsqnonlin's
%   'SpecifyObjectiveGradient' option. f_min/f_max are fixed inputs, not
%   fitted -- see vUS_2D_erf_windowed_vec_split.m's header.
%
% Inputs:
%   x: [v_xgp, v_zgp, F, DC] (DC real-valued)
%   tau: vector of time lags [s]
%   k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_z] [m]
%   f_min, f_max: sampled window bounds within [0, 2]
%
% Outputs:
%   g1_split: [numel(tau), 2] = [real(g1), imag(g1)]
%   J: [2*numel(tau), 4] Jacobian of g1_split(:) w.r.t. [v_xgp,v_zgp,F,DC]

    tau = tau(:);
    x = x(:); sigma = sigma(:);

    if nargout > 1
        [g1, Jc] = vUS_2D_erf_windowed_complex_Jac(x, tau, k0, sigma, f_min, f_max);
        J = [real(Jc); imag(Jc)];
    else
        g1 = vUS_2D_erf_windowed_complex_Jac(x, tau, k0, sigma, f_min, f_max);
    end

    g1 = g1(:);
    g1_split = [real(g1), imag(g1)];
end
