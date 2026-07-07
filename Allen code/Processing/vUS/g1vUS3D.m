%% Description:
%   Calculate the g1 model as in the vUS paper (Eq. 15 in Jianbo Tang et al., 2020)

% Ex:
% k0 = 2*pi/(1540/13.6e6); [rad/m]
% sigma = [379, 379, 111].*1e-6; [m]
% tau = (0:0.1:100).*1e-3; [s]
% test = g1vUS3D(0.7, 0.5, 0.001, 0.001, 0.020, sigma, k0, tau);

function [g1] = g1vUS3D(F, p, v_xgp, v_ygp, v_zgp, sigma, k0, tau)

    g1 = ( F.*exp(-(v_xgp .* tau).^2 ./ (4 * sigma(1)^2) ...
        - (v_ygp .* tau).^2 ./ (4 * sigma(2)^2) ...
        - (v_zgp .* tau).^2 ./ (4 * sigma(3)^2)) ...
        .* exp(-(p .* v_zgp .* k0 .* tau).^2) ...
        .* exp(2.*1i.*k0.*tau.*v_zgp) );
end