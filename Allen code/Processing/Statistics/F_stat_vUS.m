%% Calculate F-statistic of two models: one a simplified (nested) version of the other
% This version is specific to how I calculate the fits

% Inputs:
%   modelS: a struct for the simplified model
%   modelF: a struct for the full model
%   pS: # of free parameters in the simplified model
%   pF: # of free parameters in the full model
%   n: (# voxels, 1) vector of how many sample points were used in that voxel's fit
function F = F_stat_vUS(modelS, modelF, pS, pF, n)

    F = (modelS.cost - modelF.cost)./modelF.cost .* (n - pF) ./ (pF - pS);

end
