function g1 = vUS_3D_quad_RCA_vec(x, tau, k0, sigma, s, w)
%% Description:
%   g1 model for the row-column array (RCA) PSF, evaluated with fixed
%   Gauss-Legendre quadrature over s = (r/R)^2 like vUS_3D_quad_vec.m. The
%   difference is the single-scatterer g1 (g1,s), which comes from the RCA PSF
%   model (sum of two "orthogonal" lateral Gaussians x an axial Gaussian):
%
%       PSF = A [ G(x;sn) G(y;sw) + G(x;sw) G(y;sn) ] G(z;sa),  G(u;s) = exp(-u^2/(2 s^2))
%
%   For a scatterer group velocity v = f*v_gp (see vUS_3D_quad_vec.m):
%
%       g1,s(tau, f) = exp(i 2 k0 v_z tau f) * exp(-v_z^2 tau^2 f^2 / (4 sa^2))
%                      * [ (E1 + E4)/2 + c*E23 ] / (1 + c)
%       E1  = exp(-tau^2 f^2 ( v_x^2/(4 sn^2) + v_y^2/(4 sw^2) ))   (narrow-in-x term, both times)
%       E4  = exp(-tau^2 f^2 ( v_x^2/(4 sw^2) + v_y^2/(4 sn^2) ))   (wide-in-x term, both times)
%       E23 = exp(-tau^2 f^2 ( v_x^2 + v_y^2 ) / (2 (sn^2 + sw^2)))  (narrow/wide cross terms)
%       c   = 2 sn sw / (sn^2 + sw^2)
%
%   and g1(tau) = DC + F * int_0^1 g1,s(tau, f(s)) ds, with the same f(s) as
%   vUS_3D_quad_vec.m. Each of the three terms is a Gaussian in f, i.e., the
%   old integrand with its own M:
%       M_1  = v_x^2/sn^2 + v_y^2/sw^2 + v_z^2/sa^2
%       M_4  = v_x^2/sw^2 + v_y^2/sn^2 + v_z^2/sa^2
%       M_23 = 2 (v_x^2 + v_y^2)/(sn^2 + sw^2) + v_z^2/sa^2
%
%   Unlike the single-Gaussian model, g1 depends on the direction of the
%   transverse flow relative to the array axes (through phi below), so the
%   transverse-combined form with only v_tgp is not exact. g1 is unchanged by
%   phi -> -phi and phi -> pi/2 - phi (x <-> y symmetry of the RCA PSF and
%   sign symmetry), so phi is only identifiable on [0, pi/4]. The short-lag
%   decay is independent of phi; phi only changes the shape of the tail.
%
%   For sn == sw this reduces exactly to vUS_3D_quad_vec.m (any phi).
%
% Inputs:
%   x: [v_tgp, v_zgp, F, DC, k, a, phi] (DC real-valued)
%       v_tgp: transverse group velocity [m/s], v_tgp^2 = v_xgp^2 + v_ygp^2
%       v_zgp: axial group velocity [m/s]
%       F: dynamic fraction;  DC: static fraction
%       k: profile bluntness (typically 2 to 3);  a in [0, 1]: wall-speed coefficient
%       phi: direction of the transverse flow in the x-y plane, measured from the
%           x axis [rad]: v_xgp = v_tgp cos(phi), v_ygp = v_tgp sin(phi)
%   tau: time lags [s]
%   k0: wavenumber [rad/m]
%   sigma: [sigma_narrow, sigma_wide, sigma_axial] [m], FIELD (amplitude) sigmas of the
%       RCA PSF model, from fit_RCA_PSF_v2.m
%   s, w: Gauss-Legendre nodes [1, N] and weights [N, 1] on [0, 1], from
%       gaussLegendre01.m
%
% Outputs:
%   g1: [numel(tau), 1] complex

    tau = tau(:);
    s = s(:).';
    w = w(:);

    v_tgp = x(1); v_zgp = x(2); F = x(3); DC = x(4); k = x(5); a = x(6); phi = x(7);
    sn = sigma(1); sw = sigma(2); sa = sigma(3);

    v_xgp2 = (v_tgp*cos(phi))^2;
    v_ygp2 = (v_tgp*sin(phi))^2;

    S  = sn^2 + sw^2;
    c  = 2*sn*sw/S;
    Mz = v_zgp^2/sa^2;
    M1  = v_xgp2/sn^2 + v_ygp2/sw^2 + Mz;
    M4  = v_xgp2/sw^2 + v_ygp2/sn^2 + Mz;
    M23 = 2*(v_xgp2 + v_ygp2)/S + Mz;

    fmax = (k + 2) / (k + 2 - 2*a^k);
    f = fmax * (1 - a^k * s.^(k/2));                                   % [1, N]
    T2f2 = tau.^2 .* f.^2;                                             % [nTau, N]

    E = exp(2i*k0*v_zgp .* tau .* f) .* ...
        ( (exp(-M1/4.*T2f2) + exp(-M4/4.*T2f2))/2 + c*exp(-M23/4.*T2f2) ) / (1 + c);

    g1 = DC + F * (E * w);
end
