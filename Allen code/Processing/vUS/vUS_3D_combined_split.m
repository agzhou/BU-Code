function [g1_split] = vUS_3D_combined_split(x, tau, k0)
%% Description:
%   Combined parameter (C) model, split into real and imaginary components
%
% Inputs:
%   x: [C, v_zgp, F, DC] (DC real-valued)
%       C: combined parameter
%       v_zgp: group velocity components [m/s]
%       F: dynamic fraction
%       DC: static/offset component (real-valued)
%   tau: vector of time lags [s]
%   k0: wavenumber [rad/m]
%
% Outputs:
%   g1_split: [numel(tau), 2] = [real(g1), imag(g1)]

    C = x(1); v_zgp = x(2); F = x(3); DC = x(4);
    tau = tau(:);


    g1 = DC + F .* exp(-C.*tau.^2) .* exp(1i .* 2 .* k0 .* v_zgp .* tau);

    g1 = g1(:);
    g1_split = [real(g1), imag(g1)];

end
