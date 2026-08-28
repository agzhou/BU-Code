%% Description:
%   Calculate the g1 model as in the vUS paper (Eq. 15 in Jianbo Tang et al., 2020)

% Ex:
% k0 = 2*pi/(1540/13.6e6); [rad/m]
% sigma = [379, 379, 111].*1e-6; [m]
% tau = (0:0.1:100).*1e-3; [s]
% test = g1vUS3D(0.7, 0.5, 0.001, 0.001, 0.020, sigma, k0, tau);

function [g1] = g1vUS1D(F, v, sigma_v, sigma, k0, tau)
    % F = 1; % No need for F if we start at tau = 0 --> |g1(0)| = 1
    g1 = ( F.*exp(-(v .* tau).^2 ./ (4 * sigma^2)) ...
        .* exp(-(sigma_v .* k0 .* tau).^2) ...
        .* exp(2.*1i.*k0.*tau.*v) );
end