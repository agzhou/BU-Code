% %% Add path to the erfz code -- error function with complex inputs
% codeDir = cd;
% codeDir_split = split(string(codeDir), filesep);
% % AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
% ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
% addpath(genpath(ErrorFunctionCodePath))

% Description: calculate the integrand corresponding to the new g1 model as derived from the Tangelder et
% al. 1986 flow velocity profile, but the general form

% Inputs:
%   x: vector of model parameters
%       v_xgp: x group velocity [m/s]
%       v_zgp: z group velocity [m/s]
%       F: dynamic fraction
%       DC: static fraction
%       k: characterize the "bluntness" of the flow profile
%       a: [0, 1] characterize the wall flow velocity
%   tau: vector of time lags
%   k0: wavenumber [rad/m]
%   sigma: vector of sigma values (sigma_x, sigma_z) [m]

% Outputs:
%   Ifh: integrand as a function handle (wrt f, see derivation)

%%
function [Ifh] = vUS_2D_num_vec(x, tau, k0, sigma)
    v_xgp = x(1); v_zgp = x(2); F = x(3); DC = x(4); k = x(5); a = x(6);
    M = v_xgp.^2./sigma(1)^2 + v_zgp.^2./sigma(2)^2;


    Ifh = @(f) 2./(a.^2 .* k) .* (1 - 2./(k+2).*a.^k) .* exp(-M.*tau.^2 ./4 .* f.^2 + 2.*1i.*k0.*v_zgp.*tau .*f) .* (1 - f.* (1 - 2./(k+2).*a.^k)).^ (2./k - 1);

end

