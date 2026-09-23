function g1 = vUS_3D_quad_vec(x, tau, k0, sigma, s, w)
%% Description:
%   Corrected Tangelder-profile g1 model, evaluated with fixed Gauss-Legendre
%   quadrature over s = (r/R)^2 (the fraction of the vessel cross-section
%   inside radius r). Hand-written and independent of the symbolic
%   Jacobian code, so it can be used to check vUS_3D_quad_complex_Jac.m.
%
%       g1(tau) = DC + F * int_0^1 exp( -M tau^2 f(s)^2/4 + i 2 k0 v_zgp tau f(s) ) ds
%       f(s)    = fmax * (1 - a^k s^(k/2)),   fmax = (k+2) / (k+2-2 a^k)
%       M       = v_tgp^2/sigma_t^2 + v_zgp^2/sigma_z^2
%
%   g1(0) = DC + F for every a and k (the previous vUS_3D_num_vec.m form gave
%   DC + F/a^2). a = 0 is plug flow (f == 1); a = 1 is the k-only model.
%
% Inputs:
%   x: [v_tgp, v_zgp, F, DC, k, a] (DC real-valued)
%       v_tgp: transverse group velocity [m/s], v_tgp^2 = v_xgp^2 + v_ygp^2
%       v_zgp: axial group velocity [m/s]
%       F: dynamic fraction;  DC: static fraction
%       k: profile bluntness (typically 2 to 3);  a in [0, 1]: wall-speed coefficient
%   tau: time lags [s]
%   k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_y, sigma_z] [m]; sigma_x must equal sigma_y
%   s, w: Gauss-Legendre nodes [1, N] and weights [N, 1] on [0, 1], from
%       gaussLegendre01.m
%
% Outputs:
%   g1: [numel(tau), 1] complex

    if sigma(1) ~= sigma(2)
        error('vUS_3D_quad_vec:sigmaMismatch', ...
            'For the transverse-combined model, sigma_x and sigma_y must be equal.')
    end

    tau = tau(:);
    s = s(:).';
    w = w(:);

    v_tgp = x(1); v_zgp = x(2); F = x(3); DC = x(4); k = x(5); a = x(6);
    M = v_tgp^2/sigma(1)^2 + v_zgp^2/sigma(3)^2;

    fmax = (k + 2) / (k + 2 - 2*a^k);
    f = fmax * (1 - a^k * s.^(k/2));                                   % [1, N]
    E = exp(-M/4 .* tau.^2 .* f.^2 + 2i*k0*v_zgp .* tau .* f);         % [nTau, N]

    g1 = DC + F * (E * w);
end
