% Description: Wrapper function to get the objective function for
% lsqnonlin, in the 2D vUS fitting with my new g1 model, using an
% analytical Jacobian

% Optional input: OF_weight: [nTau, 1]

function [obj_fun, J] = vUS_2D_OF_nonsplit(x, tau, k0, sigma, g1_exp, OF_weight)
    if nargin > 6
        if all(size(OF_weight) ~= size(tau)) % Check the size of tau vs. OF_weight
            error("Weight vector must be the same size as tau")
        end
    end

    [g1, J] = vUS_2D_erf_complex_Jac(x, tau, k0, sigma);

    if nargin > 5 % Use the objective function weighting
        obj_fun = (g1 - g1_exp).*(OF_weight); % Weighted normal residuals for the objective function
    else % Don't use the objective function weighting
        obj_fun = g1 - g1_exp; % Normal residuals for the objective function
    end

end