%% Description:
%   Calculate the g1 model as in the vUS paper (Eq. 15 in Jianbo Tang et al., 2020)
%   This is just a version where parameters are a vector input, which is
%   what fitting functions require

% Inputs:
%   x: v_xgp, v_ygp, v_zgp, p, F (all in SI units)
%   tau: time lag vector [any units, probably seconds]
%   sigma: 1/e values of the PSF shape in x, y, z [m]
%   k0: Wavenumber [rad/m]
%   useF: boolean value --> consider the F parameter or not

% Outputs:
%   g1_split: g1(tau) but split into its real and imaginary components: has shape [nTau, 2]

function [g1_split] = g1vUS3D_vec_split(x, tau, sigma, k0, useF)
% function [g1_split] = g1vUS3D_vec_split(x, tau, sigma, k0, useDC, useF)

    v_xgp = x(1);
    v_ygp = x(2);
    v_zgp = x(3);
    p = x(4);

    % if useDC
    %     DC = x(5); % DC component, can set this to 0 if you don't want to consider that
    % else
        DC = 0;
    % end

    if useF
        % if useDC
        %     F = x(6);
        % else
            F = x(5);
        % end
    else
        F = 1; % No need for F if we start at tau = 0 --> |g1(0)| = 1
    end

    % g1_real = DC + F.*exp(-(v_xgp .* tau).^2 ./ (4 * sigma(1)^2) ...
    %     - (v_ygp .* tau).^2 ./ (4 * sigma(2)^2) ...
    %     - (v_zgp .* tau).^2 ./ (4 * sigma(3)^2)) ...
    %     .* exp(-(p .* v_zgp .* k0 .* tau).^2) ...
    %     .* cos(2.*k0.*tau.*v_zgp);
    % g1_imag = F.*exp(-(v_xgp .* tau).^2 ./ (4 * sigma(1)^2) ...
    %     - (v_ygp .* tau).^2 ./ (4 * sigma(2)^2) ...
    %     - (v_zgp .* tau).^2 ./ (4 * sigma(3)^2)) ...
    %     .* exp(-(p .* v_zgp .* k0 .* tau).^2) ...
    %     .* sin(2.*k0.*tau.*v_zgp);

    % Magnitude only
    g1_real = DC + F.*exp(-(v_xgp .* tau).^2 ./ (4 * sigma(1)^2) ...
        - (v_ygp .* tau).^2 ./ (4 * sigma(2)^2) ...
        - (v_zgp .* tau).^2 ./ (4 * sigma(3)^2)) ...
        .* exp(-(p .* v_zgp .* k0 .* tau).^2);
    g1_imag = zeros(size(g1_real));

    g1_split = [g1_real(:), g1_imag(:)];
end