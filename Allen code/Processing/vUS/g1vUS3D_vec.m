%% Description:
%   Calculate the g1 model as in the vUS paper (Eq. 15 in Jianbo Tang et al., 2020)
%   This is just a version where parameters are a vector input, which is
%   what fitting functions require

% Inputs:
%   x: v_xgp, v_ygp, v_zgp, p(all in SI units)
%   tau: time lag vector [any units, probably seconds]
%   sigma: 1/e values of the PSF shape in x, y, z [m]
%   k0: Wavenumber [rad/m]
%   useF: boolean value --> consider the F parameter or not

% Outputs:
%   g1(tau)

function [g1] = g1vUS3D_vec(x, tau, sigma, k0, useF)

    v_xgp = x(1);
    v_ygp = x(2);
    v_zgp = x(3);
    p = x(4);
    if useF
        F = x(5);
    else
        F = 1; % No need for F if we start at tau = 0 --> |g1(0)| = 1
    end

    g1 = ( F.*exp(-(v_xgp .* tau).^2 ./ (4 * sigma(1)^2) ...
        - (v_ygp .* tau).^2 ./ (4 * sigma(2)^2) ...
        - (v_zgp .* tau).^2 ./ (4 * sigma(3)^2)) ...
        .* exp(-(p .* v_zgp .* k0 .* tau).^2) ...
        .* exp(2.*1i.*k0.*tau.*v_zgp) );
end