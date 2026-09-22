function [g1] = vUS_3D_erf_TC_vec(x, tau, k0, sigma)
%% Description:
%   TC (Transverse-Combined) analog of vUS_2D_erf_vec_split.m /
%   vUS_3D_erf_vec_split.m: the new g1 model derived from the Poiseuille
%   flow model (uniform velocity probability distribution), assuming
%   sigma_x = sigma_y and combining the two transverse components into a
%   single v_tgp via v_tgp^2 = v_xgp^2 + v_ygp^2 (see vUS_3D_erf_TC.m).
%   Extended from that bare 2-parameter model with a multiplicative
%   dynamic fraction F and an additive static offset DC, packaged as a
%   single x-vector input and real/imaginary-split output for use with
%   lsqnonlin (matching vUS_2D_erf_vec_split.m / vUS_3D_erf_vec_split.m's
%   convention).
%
% Inputs:
%   x: [v_tgp, v_zgp, F, DC] (DC real-valued)
%       v_tgp: transverse group velocity [m/s], v_tgp^2 = v_xgp^2+v_ygp^2
%       v_zgp: axial group velocity [m/s]
%       F: dynamic fraction
%       DC: static/offset component (real-valued)
%   tau: vector of time lags [s]
%   k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_y, sigma_z], 1/e PSF widths [m] -- sigma_x and
%       sigma_y must be equal (checked below), matching vUS_3D_erf_TC.m
%
% Outputs:
%   g1_split: [numel(tau), 2] = [real(g1), imag(g1)]

    if sigma(1) ~= sigma(2)
        error('vUS_3D_erf_TC_vec_split:sigmaMismatch', ...
            'For this version of the g1 model, sigma_x and sigma_y must be equal.')
    end

    v_tgp = x(1); v_zgp = x(2); F = x(3); DC = x(4);
    tau = tau(:);

    M = v_tgp.^2./sigma(1)^2 + v_zgp.^2./sigma(3)^2;

    g1 = DC + F .* 1/2 .* sqrt(pi./M)./tau .* exp(-4 .* k0^2 .* v_zgp.^2 ./ M) .* ...
         ( erfz(sqrt(M).*tau - 2.*1i.*k0.*v_zgp ./ sqrt(M)) - ...
         erfz(-2.*1i.*k0.*v_zgp ./ sqrt(M)));

end
