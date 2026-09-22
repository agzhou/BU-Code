function [g1_split, J] = vUS_3D_erf_vec_split_Jac(x, tau, k0, sigma)
%% Description:
%   3D analog of vUS_2D_erf_vec_split_Jac.m. Value + analytic Jacobian of
%   the new (Poiseuille / uniform-velocity-distribution) g1 model,
%   packaged to match vUS_3D_erf_vec_split.m's [nTau, 2] real/imaginary-
%   split output convention, for use with lsqnonlin's
%   'SpecifyObjectiveGradient' option.
%
%   The actual derivative formulas live in vUS_3D_erf_complex_Jac.m (see
%   generate_vUS_3D_erf_Jac.m to (re)build it symbolically). This file
%   just splits that complex result into the real/imag-stacked form
%   lsqnonlin needs -- differentiation w.r.t. a real parameter commutes
%   with real()/imag(), so no further derivative work happens here.
%
% Inputs:
%   x: [v_xgp, v_ygp, v_zgp, F, DC] (DC real-valued)
%   tau: vector of time lags [s]
%   k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_y, sigma_z] [m]
%
% Outputs:
%   g1_split: [numel(tau), 2] = [real(g1), imag(g1)]
%   J: [2*numel(tau), 5] Jacobian of g1_split(:) w.r.t.
%      [v_xgp,v_ygp,v_zgp,F,DC] (real-part rows stacked on top of
%      imaginary-part rows, matching g1_split(:)'s column-major order).
%      Only computed if requested.

    tau = tau(:);
    % vUS_3D_erf_complex_Jac.m (symbolic-generated) indexes x/sigma as
    % column vectors internally; force column here so this wrapper
    % tolerates either orientation from the caller.
    x = x(:); sigma = sigma(:);

    if nargout > 1
        [g1, Jc] = vUS_3D_erf_complex_Jac(x, tau, k0, sigma);
        J = [real(Jc); imag(Jc)];
    else
        g1 = vUS_3D_erf_complex_Jac(x, tau, k0, sigma);
    end

    g1 = g1(:);
    g1_split = [real(g1), imag(g1)];
end
