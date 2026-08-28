% %% Add path to the erfz code -- error function with complex inputs
% codeDir = cd;
% codeDir_split = split(string(codeDir), filesep);
% % AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
% ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
% addpath(genpath(ErrorFunctionCodePath))

%%
function [g1] = vUS_1D_erf(tau, k0, sigma, v_max)
    g1 = sqrt(pi).*sigma./(v_max.*tau) .* exp(-4*k0^2*sigma^2) .* ...
         ( erfz(tau.* (v_max/(2*sigma)) + 1i*2*k0*sigma) - erfz(1i*2*k0*sigma));

end

