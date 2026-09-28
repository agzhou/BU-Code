%% Description: calculate coefficient of determination for some model data, actual data, and the dimension corresponding to samples

function [R2_adj] = calcR2_adj(model, data, dim, p)
    SS_res = sum(abs( (model - data) ).^2, dim); % Sum of squares of the residuals (wrt the model)
    SS_tot = sum( abs(data - mean(data, dim)).^2, dim); % Sum of squares: total
    n = size(data, dim);
    
    R2_adj = 1 - (n-1)/(n-p).*SS_res ./ SS_tot;

end