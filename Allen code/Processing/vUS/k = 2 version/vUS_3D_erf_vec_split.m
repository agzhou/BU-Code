function [g1_split] = vUS_3D_erf_vec_split(x, tau, k0, sigma)
%% Description:
%   3D analog of vUS_2D_erf_vec_split.m: the new g1 model derived from
%   the Poiseuille flow model (uniform velocity probability
%   distribution), extended from vUS_3D_erf.m with a multiplicative
%   dynamic fraction F and an additive static offset DC, packaged as a
%   single x-vector input and real/imaginary-split output for use with
%   lsqnonlin (matching vUS_2D_erf_vec_split.m's convention).
%
% Inputs:
%   x: [v_xgp, v_ygp, v_zgp, F, DC] (DC real-valued)
%       v_xgp, v_ygp, v_zgp: group velocity components [m/s]
%       F: dynamic fraction
%       DC: static/offset component (real-valued)
%   tau: vector of time lags [s]
%   k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_y, sigma_z], 1/e PSF widths [m]
%
% Outputs:
%   g1_split: [numel(tau), 2] = [real(g1), imag(g1)]

    v_xgp = x(1); v_ygp = x(2); v_zgp = x(3); F = x(4); DC = x(5);
    tau = tau(:);

    M = v_xgp.^2./sigma(1)^2 + v_ygp.^2./sigma(2)^2 + v_zgp.^2./sigma(3)^2;

    g1 = DC + F .* 1/2 .* sqrt(pi./M)./tau .* exp(-4 .* k0^2 .* v_zgp.^2 ./ M) .* ...
         ( erfz(sqrt(M).*tau - 2.*1i.*k0.*v_zgp ./ sqrt(M)) - ...
         erfz(-2.*1i.*k0.*v_zgp ./ sqrt(M)));

    g1 = g1(:);
    g1_split = [real(g1), imag(g1)];

end
