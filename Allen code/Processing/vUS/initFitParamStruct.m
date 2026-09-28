% Description: initialize fit parameter structs for vUS fitting

% Inputs:
%   fpn: fit parameter names (cell array of strings, which become the struct's field names)
%   initialValues: Initial value matrix (probably just zeros)


function [fit_struct] = initFitParamStruct(fpn, initialValues)
    for fpi = 1:numel(fpn)
        fit_struct.(fpn{fpi} + "_stacked") = initialValues;
    end

    for fpi = 1:numel(fpn)
        fit_struct.(fpn{fpi}) = initialValues;
    end

end