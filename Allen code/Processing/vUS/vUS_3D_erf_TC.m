% %% Add path to the erfz code -- error function with complex inputs
% codeDir = cd;
% codeDir_split = split(string(codeDir), filesep);
% % AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
% ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
% addpath(genpath(ErrorFunctionCodePath))

% Description: calculate the new g1 model as derived from the Poiseuille
% flow model, e.g., uniform velocity probability distribution.

% This Transverse velocity Combined version assumes sigma_x = sigma_y,
% and defines v_xgp^2 + v_ygp^2 = v_transverse,gp^2

% Inputs:
%   tau: vector of time lags
%   k0: wavenumber [rad/m]
%   sigma: vector of sigma values (sigma_x, sigma_y, sigma_z) [m]
%   v_tgp: transverse group velocity [m/s]
%   v_zgp: z group velocity [m/s]
   
%%
function [g1] = vUS_3D_erf_TC(tau, k0, sigma, v_tgp, v_zgp)
    % Check if the inputs are correct
    if sigma(1) ~= sigma(2)
        error('For this version of the g1 model, sigma_x and sigma_y must be equal.')
    end
    M = v_tgp.^2./sigma(1)^2 + v_zgp.^2./sigma(3)^2;

    % g1 = -1/2 .* sqrt(pi./M)./tau .* exp(-4 .* k0^2 .* v_zgp^2 ./ M) .* ...
    %      ( erfz(sqrt(M).*tau - 2.*1i.*k0.*v_zgp ./ sqrt(M)) + ...
    %      erfz(2.*1i.*k0.*v_zgp ./ sqrt(M)));
    g1 = 1/2 .* sqrt(pi./M)./tau .* exp(-4 .* k0^2 .* v_zgp.^2 ./ M) .* ...
         ( erfz(sqrt(M).*tau - 2.*1i.*k0.*v_zgp ./ sqrt(M)) - ...
         erfz(-2.*1i.*k0.*v_zgp ./ sqrt(M)));

end

