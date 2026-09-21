function [g1] = vUS_2D_num_wrapper(x, tau, k0, sigma)
    % Parse the input
    v_xgp = x(1); v_zgp = x(2); F = x(3); DC = x(4); k = x(5); a = x(6);
    
    [fmin, fmax] = calc_f_integration_limits(k, a);
    g1 = DC + F .* integral(vUS_2D_num_vec(x, tau, k0, sigma), fmin, fmax, 'ArrayValued', true);
end