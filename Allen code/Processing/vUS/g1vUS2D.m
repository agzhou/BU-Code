%% Description:
%   Calculate the g1 model as in the vUS paper (Eq. 15 in Jianbo Tang et al., 2020)

% Inputs:
%   F: dynamic fraction
%   tau: vector of time lags
%   k0: wavenumber [rad/m]
%   sigma: vector of sigma values (sigma_x, sigma_z) [m]
%   v_xgp: x group velocity [m/s]
%   v_zgp: z group velocity [m/s]

% Ex:
% k0 = 2*pi/(1540/13.6e6); [rad/m]
% sigma = [58.9110, 73.9967].*1e-6; [m]
% tau = (0:0.1:100).*1e-3; [s]
% test = g1vUS3D(0.7, 0.5, 0.001, 0.020, sigma, k0, tau);

function [g1] = g1vUS2D(F, v_xgp, v_zgp, p, sigma, k0, tau)
% function [g1] = g1vUS3D(v_xgp, v_ygp, v_zgp, p, sigma, k0, tau)
    % F = 1; % No need for F if we start at tau = 0 --> |g1(0)| = 1
    g1 = ( F.*exp(-(v_xgp .* tau).^2 ./ (4 * sigma(1)^2) ...
        - (v_zgp .* tau).^2 ./ (4 * sigma(2)^2)) ...
        .* exp(-(p .* v_zgp .* k0 .* tau).^2) ...
        .* exp(2.*1i.*k0.*tau.*v_zgp) );
end