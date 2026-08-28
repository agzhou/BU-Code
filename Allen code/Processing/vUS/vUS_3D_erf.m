% %% Add path to the erfz code -- error function with complex inputs
% codeDir = cd;
% codeDir_split = split(string(codeDir), filesep);
% % AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
% ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
% addpath(genpath(ErrorFunctionCodePath))

% Description: calculate the new g1 model as derived from the Poiseuille
% flow model, e.g., uniform velocity probability distribution

% Inputs:
%   tau: vector of time lags
%   k0: wavenumber [rad/m]
%   sigma: vector of sigma values (sigma_x, sigma_y, sigma_z) [m]
%   v_xgp: x group velocity [m/s]
%   v_ygp: y group velocity [m/s]
%   v_zgp: z group velocity [m/s]
   
%%
function [g1] = vUS_3D_erf(tau, k0, sigma, v_xgp, v_ygp, v_zgp)
    M = v_xgp.^2./sigma(1)^2 + v_ygp.^2./sigma(2)^2 + v_zgp.^2./sigma(3)^2;

    g1 = -sqrt(pi./M)./tau .* exp(-4 .* k0^2 .* v_zgp^2 ./ M) .* ...
         ( erfz(sqrt(M).*tau - 2.*1i.*k0.*v_zgp ./ sqrt(M)) - ...
         erfz(2.*1i.*k0.*v_zgp ./ sqrt(M)));

end

