% Simple helper function to calculate the limits of integration, for my
% general version of the g1 model with numerical integration

function [fmin, fmax] = calc_f_integration_limits(k, a)

    fmin = (1-a.^k) .* (k+2)./(k+2-2.*a.^k);
    fmax = (k+2)./(k+2-2.*a.^k);

end