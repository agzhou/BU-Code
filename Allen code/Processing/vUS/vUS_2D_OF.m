% Description: Wrapper function to get the objective function for
% lsqnonlin, in the 2D vUS fitting with my new g1 model, using an
% analytical Jacobian

% Optional input: OF_weight: [nTau, 1]

function [obj_fun, J] = vUS_2D_OF(x, tau, k0, sigma, g1_exp_split, OF_weight)
    if nargin > 6
        if all(size(OF_weight) ~= size(tau)) % Check the size of tau vs. OF_weight
            error("Weight vector must be the same size as tau")
        end
    end

    [g1_split, J] = vUS_2D_erf_vec_split_Jac(x, tau, k0, sigma);

    if nargin > 5 % Use the objective function weighting
        obj_fun = (g1_split - g1_exp_split).*(OF_weight); % Weighted normal residuals for the objective function
    else % Don't use the objective function weighting
        obj_fun = g1_split - g1_exp_split; % Normal residuals for the objective function
    end

end