% Description: Wrapper function to get the objective function for
% lsqnonlin, in the 3D vUS fitting with my new g1 model, with the numerically
% integrated general form (k, a free), and the Transverse velocities Combined

% Optional input: OF_weight: [nTau, 1]

function [obj_fun] = vUS_3D_num_OF(x, tau, k0, sigma, g1_exp_split, OF_weight)
    if nargin > 5
        if all(size(OF_weight) ~= size(tau)) % Check the size of tau vs. OF_weight
            error("Weight vector must be the same size as tau")
        end
    end

    if any(isnan(g1_exp_split), "all")
        error('There are NaN values in the inputted experimental g1 data')
    end

    splitOrNot = true;
    [g1_split] = vUS_3D_num_wrapper(x, tau, k0, sigma, splitOrNot);

    if nargin > 5 % Use the objective function weighting
        obj_fun = (g1_split - g1_exp_split).*(OF_weight); % Weighted normal residuals for the objective function
    else % Don't use the objective function weighting
        obj_fun = g1_split - g1_exp_split; % Normal residuals for the objective function
    end

end