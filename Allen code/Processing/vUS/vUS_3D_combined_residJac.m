function [r, J] = vUS_3D_combined_residJac(x, tau, k0, ydataReal, ydataImag, xscale)
%% Description:
%   Residual + analytic Jacobian for the combined-parameter model, for direct use with lsqnonlin.
%   Mirrors vUS_3D_quad_residJac.m / vUS_3D_quad_RCA_xy_residJac.m's convention:
%       xscale = [1, 1e-2, 1, 1];   % C in units of 1/s^2, v_zgp in units of 10 mm/s
%       fun  = @(xs) vUS_3D_combined_residJac(xs, tau, k0, yReal, yImag, xscale);
%       opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true);
%       xs = lsqnonlin(fun, x0./xscale, lb./xscale, ub./xscale, opts);
%       x  = xs .* xscale;   % x = [C, v_zgp, F, DC]
%
% Inputs:
%   x: [C, v_zgp, F, DC] divided by xscale (if xscale is given)
%   tau, k0: see vUS_3D_combined_vec.m
%   ydataReal, ydataImag: real/imag parts of the observed g1(tau), same length as tau
%   xscale: (optional) [1, 4] parameter scale; the model is evaluated at x .* xscale and J is returned
%       w.r.t. the scaled x. Default ones(1, 4).
%
% Outputs:
%   r: [2*numel(tau), 1] stacked residual, [real part; imaginary part]
%   J: [2*numel(tau), 4] analytic Jacobian of r w.r.t. x

    tau = tau(:); ydataReal = ydataReal(:); ydataImag = ydataImag(:);
    if nargin < 6 || isempty(xscale)
        xscale = ones(1, 4);
    end
    xscale = xscale(:).';
    xTrue = x(:).' .* xscale;

    if nargout > 1
        [g1, Jc] = vUS_3D_combined_complex_Jac(xTrue, tau, k0);
        Jc = Jc .* xscale;                       % chain rule for the scaling
        J = [real(Jc); imag(Jc)];
    else
        g1 = vUS_3D_combined_vec(xTrue, tau, k0);
    end

    r = [real(g1) - ydataReal; imag(g1) - ydataImag];
end
