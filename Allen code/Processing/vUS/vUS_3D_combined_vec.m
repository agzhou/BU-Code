function g1 = vUS_3D_combined_vec(x, tau, k0)
%% Description:
%   Plain complex-output counterpart of vUS_3D_combined_split.m: Jianbo's simplified model, further
%   reduced to a single combined decay parameter C (see the derivation on the "From Jianbo's simplified
%   model" / "Final model" slides -- C absorbs v_xgp, v_ygp, v_zgp and sigma_vz under the "f1, f2, f3, p
%   invariant in time" assumption):
%
%       g1(tau) = DC + F * exp(-C tau^2) * exp(i 2 k0 v_zgp tau)
%
%   Used as the independent reference for finite-difference-checking vUS_3D_combined_batchModel.m's
%   analytic Jacobian (vUS_3D_combined_complex_Jac.m and vUS_3D_combined_batchModel.m share no code with
%   this file). Returns the same g1 as vUS_3D_combined_split.m, just complex instead of [real, imag].
%
% Inputs:
%   x: [C, v_zgp, F, DC] (C >= 0, DC real-valued)
%       C: combined decay parameter [rad^2/s^2] (equivalently 1/s^2, since e^{-C tau^2} is dimensionless)
%       v_zgp: axial group velocity [m/s]
%       F: dynamic fraction;  DC: static fraction
%   tau: vector of time lags [s]
%   k0: wavenumber [rad/m]
%
% Outputs:
%   g1: [numel(tau), 1] complex

    C = x(1); v_zgp = x(2); F = x(3); DC = x(4);
    tau = tau(:);

    g1 = DC + F .* exp(-C.*tau.^2 + 1i*2*k0.*v_zgp.*tau);
end
