function [r, J] = vUS_3D_quad_RCA_xy_residJac(x, tau, k0, sigma, s, w, ydataReal, ydataImag, xscale)
%% Description:
%   Residual + analytic Jacobian for the two-step RCA-PSF g1 model (fixed-quadrature form), for direct
%   use with lsqnonlin. Direct RCA-PSF analog of vUS_3D_quad_residJac.m (7 parameters instead of 6):
%       [s, w] = gaussLegendre01(48);
%       xscale = [1e-2, 1e-2, 1e-2, 1, 1, 1, 1];   % v_xgp, v_ygp, v_zgp in units of 10 mm/s
%       fun  = @(xs) vUS_3D_quad_RCA_xy_residJac(xs, tau, k0, sigma, s, w, yReal, yImag, xscale);
%       opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true);
%       xs = lsqnonlin(fun, x0./xscale, lb./xscale, ub./xscale, opts);
%       x  = xs .* xscale;   % x = [v_xgp, v_ygp, v_zgp, F, DC, k, a]
%
%   Mirrors vUS_3D_quad_residJac.m's convention. Since the residual is model - data and data doesn't
%   depend on x, the residual's Jacobian equals the model's (vUS_3D_quad_RCA_xy_complex_Jac.m).
%
% Inputs:
%   x: [v_xgp, v_ygp, v_zgp, F, DC, k, a] divided by xscale (if xscale is given)
%   tau, k0, sigma, s, w: see vUS_3D_quad_RCA_xy_vec.m
%   ydataReal, ydataImag: real/imag parts of the observed g1(tau), same length as tau
%   xscale: (optional) [1, 7] parameter scale; the model is evaluated at x .* xscale and J is returned
%       w.r.t. the scaled x. Default ones(1, 7).
%
% Outputs:
%   r: [2*numel(tau), 1] stacked residual, [real part; imaginary part]
%   J: [2*numel(tau), 7] analytic Jacobian of r w.r.t. x

    tau = tau(:); ydataReal = ydataReal(:); ydataImag = ydataImag(:);
    if nargin < 9 || isempty(xscale)
        xscale = ones(1, 7);
    end
    xscale = xscale(:).';
    xTrue = x(:).' .* xscale;

    if nargout > 1
        [g1, Jc] = vUS_3D_quad_RCA_xy_complex_Jac(xTrue, tau, k0, sigma, s, w);
        Jc = Jc .* xscale;                       % chain rule for the scaling
        J = [real(Jc); imag(Jc)];
    else
        g1 = vUS_3D_quad_RCA_xy_vec(xTrue, tau, k0, sigma, s, w);
    end

    r = [real(g1) - ydataReal; imag(g1) - ydataImag];
end
