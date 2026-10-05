%% Calculate F-statistic of two models: one a simplified (nested) version of the other

% Inputs:
%   modelS: the simplified model's predictions
%   modelF: the full model's predictions
%   data: the actual data
%   dim: dimension corresponding to sample points
%   pS: # of free parameters in the simplified model
%   pF: # of free parameters in the full model

function F = F_stat(modelS, modelF, data, dim, pS, pF)
    n = size(data, dim); % Number of sample points
    if ~isequal(size(modelS, dim), size(modelF, dim), size(data, dim))
        error('Error in F_stat.m: input data dimensions are not consistent in the sample point dimension')
    end

    RSS_S = calcRSS(modelS, data, dim);
    RSS_F = calcRSS(modelF, data, dim);

    F = (RSS_S - RSS_F)./RSS_F .* (n - pF) ./ (pF - pS);

end

%% Helper functions
function RSS = calcRSS(model, data, dim)
    RSS = sum(abs( (model - data) ).^2, dim); % Sum of squares of the residuals (wrt the model)
end