function [g1, J] = vUS_3D_quad_complex_Jac(x, tau, k0, sigma, s, w)
%% Description:
%   [g1, J] for the corrected Tangelder-profile g1 model (fixed-quadrature
%   form), with the analytic Jacobian w.r.t. all six parameters
%   [v_tgp, v_zgp, F, DC, k, a]. Hand-written assembly around the
%   symbolically generated vUS_3D_quad_complex_Jac_raw.m (see
%   generate_vUS_3D_quad_Jac.m for the derivation).
%
%       g1  = DC + F * G,              G = E*w = int_0^1 exp(phi) ds
%       dg1/dF  = G
%       dg1/dDC = 1
%       dg1/dp  = F * (dE/dp)*w        for p in {v_tgp, v_zgp, k, a}
%
% Inputs: same as vUS_3D_quad_vec.m
%   x: [v_tgp, v_zgp, F, DC, k, a];  tau: [nTau, 1] [s];  k0: [rad/m]
%   sigma: [sigma_x, sigma_y, sigma_z] [m];  s: [1, N], w: [N, 1] from gaussLegendre01.m
%
% Outputs:
%   g1: [nTau, 1] complex
%   J:  [nTau, 6] complex Jacobian of g1 w.r.t. [v_tgp, v_zgp, F, DC, k, a]
%
% Note: at a = 0 exactly the k- and a-columns are degenerate (the profile is
% plug flow, f == 1, and k has no effect). a is clamped to realmin here so
% a^k*log(a) evaluates to 0 instead of NaN; keep a bounded away from 0 in
% the fit if k needs to be identifiable.

    if sigma(1) ~= sigma(2)
        error('vUS_3D_quad_complex_Jac:sigmaMismatch', ...
            'For the transverse-combined model, sigma_x and sigma_y must be equal.')
    end

    if nargout < 2
        g1 = vUS_3D_quad_vec(x, tau, k0, sigma, s, w);
        return
    end

    tau = tau(:);
    s = s(:).';
    w = w(:);
    x = x(:);
    sigma = sigma(:);   % the generated code indexes sigma(3,:), which needs a column
    x(6) = max(x(6), realmin);

    F = x(3); DC = x(4);
    [E, dE_dvt, dE_dvz, dE_dk, dE_da] = vUS_3D_quad_complex_Jac_raw(x, tau, k0, sigma, s);

    % Every generated output carries the exp(phi) factor, so all five are
    % full-size [nTau, N] arrays and the quadrature is a matrix-vector product.
    G     = E      * w;
    dG_vt = dE_dvt * w;
    dG_vz = dE_dvz * w;
    dG_k  = dE_dk  * w;
    dG_a  = dE_da  * w;

    g1 = DC + F * G;
    J  = [F*dG_vt, F*dG_vz, G, ones(size(G)), F*dG_k, F*dG_a];
end
