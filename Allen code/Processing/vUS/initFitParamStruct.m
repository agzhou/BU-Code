% Description: initialize fit parameter structs for vUS fitting

% Inputs:
%   fpn: fit parameter names (cell array of strings, which become the struct's field names)
%   initialValues: Initial value matrix (probably just zeros)


function [fit_struct] = initFitParamStruct(fpn, initialValues, numCopies)
    if nargin > 2
        for j = 1:numCopies % Replicate the struct some number of times
            for fpi = 1:numel(fpn)
                fit_struct(j).(fpn{fpi} + "_stacked") = initialValues;
            end
        
            for fpi = 1:numel(fpn)
                fit_struct(j).(fpn{fpi}) = initialValues;
            end
        end
    else
        for fpi = 1:numel(fpn)
            fit_struct.(fpn{fpi} + "_stacked") = initialValues;
        end
    
        for fpi = 1:numel(fpn)
            fit_struct.(fpn{fpi}) = initialValues;
        end
    end

end