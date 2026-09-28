function g1 = vUS_3D_quad_RCA_xy_vec(x, tau, k0, sigma, s, w)
%% Description:
%   Fast, fixed-quadrature counterpart of vUS_3D_num_RCA_vec.m: the two-step RCA-PSF g1 model
%   (single-scatterer g1,s from the RCA PSF -- sum of two "orthogonal" lateral Gaussians x an axial
%   Gaussian -- convolved with the general Tangelder velocity profile), evaluated with fixed
%   Gauss-Legendre quadrature over s = (r/R)^2 instead of MATLAB's adaptive integral() over f. Same
%   change of variables as vUS_3D_quad_vec.m relative to vUS_3D_num_vec.m: it removes the (1 - f*c)^(2/k-1)
%   power-law weight from the integrand (singular near f = fmax when k > 2), which is what makes
%   adaptive quadrature slow and occasionally fussy.
%
%   Model: with S = sigma_narrow^2 + sigma_wide^2, c = 2*sigma_narrow*sigma_wide/S,
%
%       g1,s(tau, f) = exp(i 2 k0 v_zgp tau f) * exp(-(v_zgp tau f)^2 / (4 sigma_axial^2)) *
%                      [ (E1 + E4)/2 + c*E23 ] / (1 + c)
%       E1  = exp(-(f v_xgp tau)^2/(4 sigma_narrow^2) - (f v_ygp tau)^2/(4 sigma_wide^2))
%       E4  = exp(-(f v_xgp tau)^2/(4 sigma_wide^2)   - (f v_ygp tau)^2/(4 sigma_narrow^2))
%       E23 = exp(-((f v_xgp tau)^2 + (f v_ygp tau)^2) / (2 S))
%       g1(tau) = DC + F * int_0^1 g1,s(tau, f(s)) ds,   f(s) = fmax*(1 - a^k s^(k/2)),  fmax = (k+2)/(k+2-2a^k)
%
%   Unlike the single-Gaussian model, this depends on v_xgp and v_ygp separately, not only on
%   v_xgp^2 + v_ygp^2 -- see RCA_g1s_derivation.md ("v_transverse cannot be cleanly compacted"). It IS
%   invariant to the sign of either component and to swapping them (v_xgp, v_ygp) <-> (v_ygp, v_xgp),
%   since it only ever depends on v_xgp^2 and v_ygp^2. For sigma_narrow == sigma_wide this reduces
%   exactly to vUS_3D_quad_vec.m with v_tgp^2 = v_xgp^2 + v_ygp^2 (checked in
%   generate_vUS_3D_quad_RCA_xy_Jac.m and test_vUS_3D_quad_RCA_xy_batch.m).
%
% Inputs:
%   x: [v_xgp, v_ygp, v_zgp, F, DC, k, a] (DC real-valued) -- same order as vUS_3D_num_RCA_vec.m
%       v_xgp, v_ygp, v_zgp: group velocity components [m/s]
%       F: dynamic fraction;  DC: static fraction
%       k: profile bluntness (typically 2 to 3);  a in [0, 1]: wall-speed coefficient
%   tau: time lags [s]
%   k0: wavenumber [rad/m]
%   sigma: [sigma_narrow, sigma_wide, sigma_axial] [m], FIELD (amplitude) sigmas of the RCA PSF model,
%       from fit_RCA_PSF_v2.m
%   s, w: Gauss-Legendre nodes [1, N] and weights [N, 1] on [0, 1], from gaussLegendre01.m
%
% Outputs:
%   g1: [numel(tau), 1] complex

    tau = tau(:);
    s = s(:).';
    w = w(:);

    v_xgp = x(1); v_ygp = x(2); v_zgp = x(3); F = x(4); DC = x(5); k = x(6); a = x(7);
    sn = sigma(1); sw = sigma(2); sa = sigma(3);

    S  = sn^2 + sw^2;
    c  = 2*sn*sw/S;
    w1 = 1/(2*(1 + c));  w23 = c/(1 + c);   % w1 = w4, and 2*w1 + w23 = 1

    Mt1  = v_xgp^2/sn^2 + v_ygp^2/sw^2;
    Mt4  = v_xgp^2/sw^2 + v_ygp^2/sn^2;
    Mt23 = 2*(v_xgp^2 + v_ygp^2)/S;
    Mz   = v_zgp^2/sa^2;

    fmax = (k + 2) / (k + 2 - 2*a^k);
    f = fmax * (1 - a^k * s.^(k/2));                                   % [1, N]
    T2f2 = tau.^2 .* f.^2;                                             % [nTau, N]

    EZ = exp(2i*k0*v_zgp .* tau .* f - (Mz/4).*T2f2);
    D1 = exp(-(Mt1/4).*T2f2); D4 = exp(-(Mt4/4).*T2f2); D23 = exp(-(Mt23/4).*T2f2);

    g1 = DC + F * ( EZ .* ( w1*(D1 + D4) + w23*D23 ) * w );
end
