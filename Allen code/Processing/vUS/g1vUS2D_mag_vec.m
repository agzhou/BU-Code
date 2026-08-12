%% Description:
%   Calculate the g1 model as in the vUS paper (Eq. 15 in Jianbo Tang et al., 2020)
%   This is just a version where parameters are a vector input, which is
%   what fitting functions require

% This outputs the magnitude only -- |g1|

% Inputs:
%   x: v_xgp, v_zgp, p, F (all in SI units)
%   tau: time lag vector [any units, probably seconds]
%   sigma: 1/e values of the PSF shape in x, z [m]
%   k0: Wavenumber [rad/m]
%   useF: boolean value --> consider the F parameter or not

% Outputs:
%   |g1(tau)|

function [g1] = g1vUS2D_mag_vec(x, tau, sigma, k0, useF, useDC)

    v_xgp = x(1);
    v_zgp = x(2);
    p = x(3);
    if useF
        F = x(4);
    else
        F = 1; % No need for F if we start at tau = 0 --> |g1(0)| = 1
    end

    if useDC
        DC = x(5);
    else
        DC = 0;
    end

    g1 = ( DC + F.*exp(-(v_xgp .* tau).^2 ./ (4 * sigma(1)^2) ...
        - (v_zgp .* tau).^2 ./ (4 * sigma(2)^2)) ...
        .* exp(-(p .* v_zgp .* k0 .* tau).^2) );
end