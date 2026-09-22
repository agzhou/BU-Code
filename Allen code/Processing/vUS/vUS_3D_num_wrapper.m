% splitOrNot: boolean to split (or not) the g1 into separate real and
% imaginary components
function [g1] = vUS_3D_num_wrapper(x, tau, k0, sigma, splitOrNot)
    % Parse the input
    v_tgp = x(1); v_zgp = x(2); F = x(3); DC = x(4); k = x(5); a = x(6);
    if nargin < 5
        splitOrNot = false;
    end
    if size(tau, 2) > 1
        error('tau must be a column vector')
    end
    
    [fmin, fmax] = calc_f_integration_limits(k, a);
    g1 = DC + F .* integral(vUS_3D_num_vec(x, tau, k0, sigma), fmin, fmax, 'ArrayValued', true);

    if splitOrNot
        g1 = g1(:);
        g1 = [real(g1), imag(g1)];
    end
end