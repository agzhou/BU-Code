function r = vUS_3D_quad_RCA_vessel_residual(x, R, tau, k0, sigma, s, w, nth, ydataReal, ydataImag, xscale)
%% Description:
%   Residual for the general (any-vessel-orientation) RCA g1 model
%   (vUS_3D_quad_RCA_vessel_vec.m), for direct use with lsqnonlin:
%       [s, w] = gaussLegendre01(48);
%       xscale = [1e-2, 1e-2, 1, 1, 1, 1, 1];   % v_tgp, v_zgp in units of 10 mm/s
%       fun  = @(xs) vUS_3D_quad_RCA_vessel_residual(xs, R, tau, k0, sigma, s, w, 32, yReal, yImag, xscale);
%       opts = optimoptions('lsqnonlin', 'Display', 'off');   % finite-difference Jacobian --
%           no analytic Jacobian is provided (see the header note below)
%       xs = lsqnonlin(fun, x0./xscale, lb./xscale, ub./xscale, opts);
%       x  = xs .* xscale;   % x = [v_tgp, v_zgp, F, DC, k, a, phi]
%
%   R (vessel radius) is a FIXED input, not fit: unlike v_tgp/v_zgp/phi, R is not something g1(tau)
%   alone identifies well (it trades off against k, a, and the dynamic fraction F), and in practice
%   it would come from a vessel segmentation / assumed value, not from fitting a single g1 curve.
%   Pass the same R used to simulate/collect the data; sweeping R and refitting is the way to check
%   sensitivity to that assumption.
%
%   No analytic Jacobian: vUS_3D_quad_RCA_vec.m's Jacobian is generated symbolically (see
%   generate_vUS_3D_quad_Jac.m) because the quadrature nodes there are parameter-free, so d/dx
%   commutes with the s-quadrature. Here the theta quadrature is ALSO parameter-free (fixed
%   trapezoid nodes), so the same commuting argument still holds in principle -- but the model
%   itself is algebraically heavier (a 2x2 sum of J_jk terms per node, each its own nested
%   exp/sqrt), so a symbolic derivation was not attempted; lsqnonlin's finite-difference Jacobian
%   is used instead. This is slower per iteration and per-voxel (not batched), suitable for single
%   fits and small-scale comparisons, not yet for whole-volume fitting.
%
% Inputs:
%   x: [v_tgp, v_zgp, F, DC, k, a, phi] divided by xscale (if xscale is given) -- same convention
%       and order as vUS_3D_quad_RCA_vec.m
%   R: vessel radius [m], fixed (not part of x)
%   tau, k0, sigma, s, w, nth: see vUS_3D_quad_RCA_vessel_vec.m
%   ydataReal, ydataImag: real/imag parts of the observed g1(tau), same length as tau
%   xscale: (optional) [1, 7] parameter scale; the model is evaluated at x .* xscale. Default ones(1, 7).
%
% Outputs:
%   r: [2*numel(tau), 1] stacked residual, [real part; imaginary part]

    tau = tau(:); ydataReal = ydataReal(:); ydataImag = ydataImag(:);
    if nargin < 11 || isempty(xscale)
        xscale = ones(1, 7);
    end
    xTrue = x(:).' .* xscale(:).';

    g1 = vUS_3D_quad_RCA_vessel_vec([xTrue, R], tau, k0, sigma, s, w, nth);

    r = [real(g1) - ydataReal; imag(g1) - ydataImag];
end
