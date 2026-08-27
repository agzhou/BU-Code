%% Description:
%   Calculate the simplified g1 model with only tau_c
%   This is just a version where parameters are a vector input, which is
%   what fitting functions require

% Inputs:
%   x: tau_c, optional: [F, DC component] (all in SI units)
%   tau: time lag vector [any units, probably seconds]
%   k0: Wavenumber [rad/m]
%   useF: boolean value --> consider the F parameter (dynamic component) or not
%   useDC: boolean value --> consider the DC parameter (static component) or not

% Outputs:
%   |g1(tau)| model

function [g1_mag] = g1ModelSimplified2D_vec(x, tau, k0, useF, useDC)

    tau_c = x(1);
    if useF
        F = x(2);
    else
        F = 1; % No need for F if we start at tau = 0 --> |g1(0)| = 1
    end

    if useDC
        DC = x(3);
    else
        DC = 0;
    end

    g1_mag = ( DC + F.*exp(-(tau).^2 ./ tau_c^2) );
        % .* exp(2.*1i.*k0.*tau.*v_zgp) );
end