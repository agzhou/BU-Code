function g1 = vUS_3D_quad_RCA_vessel_vec(x, tau, k0, sigma, s, w, nth)
%% Description:
%   General (any vessel orientation) g1 model for the RCA PSF, from the "General vessel
%   orientation" derivation in RCA_g1s_derivation.md. Generalizes vUS_3D_quad_RCA_vec.m (which
%   assumes the PSF is much wider than the vessel, so the vessel's flow profile can be averaged
%   independently of the PSF) by explicitly weighting each point of the vessel cross-section by
%   the PSF intensity there. Reduces to vUS_3D_quad_RCA_vec.m exactly as the vessel radius R grows
%   past the PSF width, and further to vUS_3D_quad_vec.m when sigma_narrow = sigma_wide (any R).
%
%   The vessel axis direction is n_hat = v_gp/|v_gp| (the local velocity is always parallel to the
%   bulk group velocity; only its magnitude scales with the flow profile f(rho)), and the vessel
%   axis is assumed to pass through the voxel center (no lateral offset). Points on the vessel
%   cross-section are written r = rho*u_hat(theta) + zeta*n_hat, for any fixed orthonormal
%   (e_xi, e_eta) spanning the plane perpendicular to n_hat; the along-axis coordinate zeta is
%   integrated out in closed form (a scatterer's motion is purely along zeta, so this is an
%   ordinary Gaussian integral, exactly Step 1/2's completing-the-square trick), leaving a 2D
%   quadrature over the vessel cross-section (rho via Gauss-Legendre in s = (rho/R)^2, theta via
%   the trapezoid rule, spectrally accurate for a smooth periodic integrand).
%
% Inputs:
%   x: [v_tgp, v_zgp, F, DC, k, a, phi, R] (DC real-valued) -- the same 7 parameters, in the same
%       order, as vUS_3D_quad_RCA_vec.m, with R appended as an 8th
%       v_tgp, v_zgp, phi: transverse group speed [m/s], axial group velocity [m/s], and the
%           transverse flow direction [rad] (v_xgp = v_tgp*cos(phi), v_ygp = v_tgp*sin(phi)) --
%           same convention as vUS_3D_quad_RCA_vec.m. n_hat = v_gp/|v_gp| is computed from these.
%       F, DC: dynamic fraction, static fraction
%       k, a: Tangelder profile bluntness and wall-speed coefficient
%       R: vessel radius [m] -- NEW relative to vUS_3D_quad_RCA_vec.m; the vessel-vs-PSF size that
%           model implicitly assumes is infinite
%   tau: time lags [s]
%   k0: wavenumber [rad/m]
%   sigma: [sigma_narrow, sigma_wide, sigma_axial] [m], field sigmas of the RCA PSF model
%   s, w: Gauss-Legendre nodes [1, N] and weights [N, 1] on [0, 1] (radial quadrature, same
%       convention/nodes as vUS_3D_quad_vec.m / vUS_3D_quad_RCA_vec.m)
%   nth: number of trapezoid-rule nodes for the theta (azimuthal) quadrature (default 32)
%
% Outputs:
%   g1: [numel(tau), 1] complex

    if nargin < 7 || isempty(nth)
        nth = 32;
    end

    tau = tau(:);
    v_tgp = x(1); v_zgp = x(2); F = x(3); DC = x(4); k = x(5); a = x(6); phi = x(7); R = x(8);
    sn = sigma(1); sw = sigma(2); sa = sigma(3);

    vgp_vec = [v_tgp*cos(phi); v_tgp*sin(phi); v_zgp];
    vgp_mag = norm(vgp_vec);
    if vgp_mag < eps % No flow at all: every point is stationary, g1,s = 1 everywhere
        g1 = (DC + F) * ones(size(tau));
        return
    end
    n_hat = vgp_vec / vgp_mag;

    % Orthonormal basis for the plane perpendicular to n_hat (arbitrary choice -- the theta
    % quadrature integrates over the full circle, so the final result doesn't depend on it)
    ref = [0; 0; 1];
    if abs(dot(ref, n_hat)) > 0.9
        ref = [1; 0; 0];
    end
    e_xi = ref - dot(ref, n_hat)*n_hat; e_xi = e_xi / norm(e_xi);
    e_eta = cross(n_hat, e_xi);

    Lam1 = diag([1/sn^2, 1/sw^2, 1/sa^2]);
    Lam2 = diag([1/sw^2, 1/sn^2, 1/sa^2]);

    % Radial (rho) and azimuthal (theta) quadrature grids, flattened to a single point list
    rho = R * sqrt(s(:)); % [N, 1]
    area_w = (R^2/2) * w(:); % ds weight for rho drho = (R^2/2) ds (per radial node)
    thg = (0:nth-1)/nth * 2*pi; % Uniform trapezoid nodes on the periodic interval [0, 2*pi)
    theta_w = (2*pi/nth); % Uniform weight (trapezoid rule on a periodic function)

    [RHO, TH] = ndgrid(rho, thg); % [N, nth]
    WPT = area_w .* ones(1, nth) * theta_w; % [N, nth] combined 2D quadrature weight
    RHO = RHO(:).'; TH = TH(:).'; WPT = WPT(:).'; % Flatten to [1, N*nth]

    u_hat = e_xi*cos(TH) + e_eta*sin(TH); % [3, N*nth]

    % K^{zeta zeta} (scalar per k), L_k(rho,theta), Q_k(rho,theta) -- see the derivation
    Kzz1 = n_hat.'*Lam1*n_hat; Kzz2 = n_hat.'*Lam2*n_hat;
    L1 = RHO .* (n_hat.'*Lam1*u_hat); L2 = RHO .* (n_hat.'*Lam2*u_hat); % [1, N*nth]
    Q1 = RHO.^2 .* sum(u_hat.*(Lam1*u_hat), 1); Q2 = RHO.^2 .* sum(u_hat.*(Lam2*u_hat), 1);

    Kzz = [Kzz1, Kzz2]; Lc = {L1, L2}; Qc = {Q1, Q2};

    % Denominator: J_jk at shift delta = 0, summed over j,k in {1,2}, quadrature-weighted
    Jden = zeros(1, size(RHO, 2));
    for j = 1:2
        for kk = 1:2
            Ajk = (Kzz(j)+Kzz(kk))/2;
            Bjk = Lc{j} + Lc{kk};
            Cjk = (Qc{j}+Qc{kk})/2;
            Jden = Jden + sqrt(pi./Ajk) .* exp(Bjk.^2./(4*Ajk) - Cjk);
        end
    end
    Denominator = sum(Jden .* WPT);

    % Numerator: J_jk at shift delta = f(rho)*|v_gp|*tau, for every (tau, point) pair
    fmax = (k + 2) / (k + 2 - 2*a^k);
    f_rho = fmax * (1 - a^k * (RHO/R).^k); % [1, N*nth]
    delta = f_rho .* vgp_mag .* tau; % [nTau, N*nth] via implicit expansion
    phase = exp(1i*2*k0*v_zgp .* f_rho .* tau); % [nTau, N*nth]

    Jnum = zeros(size(delta));
    for j = 1:2
        for kk = 1:2
            Ajk = (Kzz(j)+Kzz(kk))/2;
            Bjk = Lc{j} + Lc{kk} + Kzz(kk).*delta;
            Cjk = (Qc{j}+Qc{kk})/2 + Lc{kk}.*delta + Kzz(kk).*delta.^2/2;
            Jnum = Jnum + sqrt(pi./Ajk) .* exp(Bjk.^2./(4*Ajk) - Cjk);
        end
    end
    Numerator = sum(phase .* Jnum .* WPT, 2); % [nTau, 1]

    g1 = DC + F * (Numerator / Denominator);
end
