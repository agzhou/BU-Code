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
%       v_ygp: y group velocity [m/s]
%       v_zgp: z group velocity [m/s]
%       F: dynamic fraction
%       DC: static fraction
%       k: characterize the "bluntness" of the flow profile
%       a: [0, 1] characterize the wall flow velocity
%   tau: vector of time lags
%   k0: wavenumber [rad/m]
%   sigma: vector of sigma values (sigma_narrow, sigma_wide, sigma_axial) [m]

% Outputs:
%   Ifh: integrand as a function handle (wrt f, see derivation)

%%
function [g1] = vUS_3D_num_RCA_vec(x, tau, k0, sigma)
    
    % Parse the input
    v_xgp = x(1); v_ygp = x(2); v_zgp = x(3); F = x(4); DC = x(5); k = x(6); a = x(7);
    % M = v_tgp.^2./sigma(1)^2 + v_zgp.^2./sigma(3)^2;
    sigma_n = sigma(1); sigma_w = sigma(2); sigma_a = sigma(3); % parse sigma vector
    s = sigma_w^2 + sigma_n^2; % Make a new constant for convenience

    % Get limits of integration
    [fmin, fmax] = calc_f_integration_limits(k, a);

    % Output the g1 model integrand as an anonymous function [handle] of f
    Ifh = @(f) exp(-(f.*v_zgp.*tau).^2 ./ (4.*sigma_a.^2) + 2.*1i.*k0.*v_zgp.*tau .*f) .* (1 - f.* (1 - 2./(k+2).*a.^k)).^ (2./k - 1) .* ( exp(-(f.*v_xgp.*tau).^2 ./ (4.*sigma_n.^2) - (f.*v_ygp.*tau).^2 ./ (4.*sigma_w.^2)) + 4.*sigma_n.*sigma_w./s.* exp(-((f.*v_xgp.*tau).^2 + (f.*v_ygp.*tau).^2) ./ (2.*s)) + exp(-(f.*v_xgp.*tau).^2 ./ (4.*sigma_w.^2) - (f.*v_ygp.*tau).^2 ./ (4.*sigma_n.^2)) );

    g1 = DC + F.* 2./(a.^2 .* k) .* (1 - 2./(k+2).*a.^k) ./ (2.*(1 + 2.* sigma_n.*sigma_w./s)) .* integral(Ifh, fmin, fmax, 'ArrayValued', true);
end

