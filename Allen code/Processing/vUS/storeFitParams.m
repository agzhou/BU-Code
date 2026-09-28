% Description: Store fitted parameter in the struct (vUS fitting)

% Inputs:
%   fit_struct: (initialized) fit parameter struct; see initFitParamStruct.m
%   fpn: fit parameter names (cell array of strings, which become the struct's field names)
%   x: fitted parameter matrix (likely of size # voxels x # fit parameters)
%   fitted_inds: indices corresponding to the fits actually performed (e.g., which voxels are actually fit)
%   PP: processing parameter struct; used to determine how to unstack the fitted parameter values
function [fit_struct] = storeFitParams(fit_struct, fpn, x, fitted_inds, PP)
    for fpi = 1:numel(fpn)
        fit_struct.(fpn{fpi} + "_stacked")(fitted_inds, :) = x(:, fpi);
        fit_struct.(fpn{fpi}) = unstackData(fit_struct.(fpn{fpi} + "_stacked"), PP);
    end

    for fpi = 1:numel(fpn)
        
    end

end