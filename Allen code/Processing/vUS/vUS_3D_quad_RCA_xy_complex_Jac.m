function [g1, J] = vUS_3D_quad_RCA_xy_complex_Jac(x, tau, k0, sigma, s, w)
%% Description:
%   [g1, J] for the two-step RCA-PSF g1 model (fixed-quadrature form), with the analytic Jacobian
%   w.r.t. all seven parameters [v_xgp, v_ygp, v_zgp, F, DC, k, a]. Hand-written assembly around the
%   symbolically generated vUS_3D_quad_RCA_xy_complex_Jac_raw.m (see generate_vUS_3D_quad_RCA_xy_Jac.m
%   for the derivation). Direct RCA-PSF analog of vUS_3D_quad_complex_Jac.m.
%
%       g1  = DC + F * G,              G = I*w = int_0^1 Integrand ds
%       dg1/dF  = G
%       dg1/dDC = 1
%       dg1/dp  = F * (dI/dp)*w        for p in {v_xgp, v_ygp, v_zgp, k, a}
%
% Inputs: same as vUS_3D_quad_RCA_xy_vec.m
%   x: [v_xgp, v_ygp, v_zgp, F, DC, k, a];  tau: [nTau, 1] [s];  k0: [rad/m]
%   sigma: [sigma_narrow, sigma_wide, sigma_axial] [m];  s: [1, N], w: [N, 1] from gaussLegendre01.m
%
% Outputs:
%   g1: [nTau, 1] complex
%   J:  [nTau, 7] complex Jacobian of g1 w.r.t. [v_xgp, v_ygp, v_zgp, F, DC, k, a]
%
% Note: at a = 0 exactly the k- and a-columns are degenerate (the profile is plug flow, f == 1, and k
% has no effect). a is clamped to realmin here so a^k*log(a) evaluates to 0 instead of NaN; keep a
% bounded away from 0 in the fit if k needs to be identifiable.

    if nargout < 2
        g1 = vUS_3D_quad_RCA_xy_vec(x, tau, k0, sigma, s, w);
        return
    end

    tau = tau(:);
    s = s(:).';
    w = w(:);
    x = x(:);
    sigma = sigma(:);   % the generated code indexes sigma(1,:) etc., which needs a column
    x(7) = max(x(7), realmin);

    F = x(4); DC = x(5);
    [I, dI_dvx, dI_dvy, dI_dvz, dI_dk, dI_da] = vUS_3D_quad_RCA_xy_complex_Jac_raw(x, tau, k0, sigma, s);

    % Every generated output carries the Integrand factor, so all six are full-size [nTau, N] arrays
    % and the quadrature is a matrix-vector product.
    G      = I       * w;
    dG_vx  = dI_dvx  * w;
    dG_vy  = dI_dvy  * w;
    dG_vz  = dI_dvz  * w;
    dG_k   = dI_dk   * w;
    dG_a   = dI_da   * w;

    g1 = DC + F * G;
    J  = [F*dG_vx, F*dG_vy, F*dG_vz, G, ones(size(G)), F*dG_k, F*dG_a];
end
