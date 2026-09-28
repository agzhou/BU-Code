%% Description: calculate coefficient of determination for some model data, actual data, and the dimension corresponding to samples

function [R2] = calcR2(model, data, dim)
    SS_res = sum(abs( (model - data) ).^2, dim); % Sum of squares of the residuals (wrt the model)
    SS_tot = sum( abs(data - mean(data, dim)).^2, dim); % Sum of squares: total
    R2 = 1 - SS_res ./ SS_tot;

end