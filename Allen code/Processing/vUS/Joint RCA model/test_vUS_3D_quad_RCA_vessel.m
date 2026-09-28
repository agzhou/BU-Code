%% Description:
%   Validation of vUS_3D_quad_RCA_vessel_vec.m (general, any-vessel-orientation g1 model). Checks:
%     1. Reduces to vUS_3D_quad_RCA_vec.m as R shrinks well below the PSF, for AXIAL flow
%        (v_tgp = 0). This is the case worked out in closed form in the derivation.
%     2. Reduces to vUS_3D_quad_vec.m as R shrinks, for a MIXED flow direction, when
%        sigma_narrow = sigma_wide (isotropic PSF). Isotropy alone is enough, no axis alignment
%        needed: the (sigma_narrow - sigma_wide) factor that keeps the two-step model exact for any
%        orientation vanishes here.
%     2b. Negative result, documented rather than asserted: for a MIXED flow direction (both
%        transverse and axial) with the ACTUAL anisotropic sigma, the joint model does NOT reduce
%        to vUS_3D_quad_RCA_vec.m even as R -> 0. Physically, R -> 0 collapses the vessel to a
%        1-D filament along the (tilted) flow axis, not to "the whole PSF moving uniformly" (what
%        vUS_3D_quad_RCA_vec.m's own g1,s assumes) -- these only coincide when the flow axis is
%        aligned with an image axis or the PSF is laterally isotropic. This is expected, not a bug
%        (see the brute-force checks in #3, which hold for this exact case at the actual R).
%     3. Matches independent brute-force 3D numerical integration (integral3, no reference to the
%        model's own closed-form quadrature) for several vessel orientations, including purely
%        transverse flow, AT THE ACTUAL (non-shrunk) R -- the check that matters for real use
%     4. g1(0) = DC + F

clearvars
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
addpath(genpath(fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing\")))

sn = 57.6e-6; sw = 286.5e-6; sa = 50.1e-6;
sigma = [sn, sw, sa];
k0 = 2*pi / (1540/13.6e6);
[s, w] = gaussLegendre01(60);
tol = 1e-6;
allPass = true;
report = @(name, err, tolv) fprintf('%-70s %10.2e  %s\n', name, err, ternary(err < tolv, 'pass', 'FAIL'));

%% 1. Axial flow (v_tgp = 0): reduces to vUS_3D_quad_RCA_vec.m as R shrinks below the PSF
tau = (0:40)'/2500;
x7 = [0, 5e-3, 0.8, 0.10, 2.5, 0.7, 0]; % [v_tgp, v_zgp, F, DC, k, a, phi], v_tgp = 0
g_twostep = vUS_3D_quad_RCA_vec(x7, tau, k0, sigma, s, w);
for R = [20e-6, 5e-6, 1e-6]
    g_vessel = vUS_3D_quad_RCA_vessel_vec([x7, R], tau, k0, sigma, s, w, 32);
    e = max(abs(g_twostep - g_vessel));
    allPass = allPass && e < 1e-2 * (R/20e-6); % Tighten the tolerance as R shrinks
    report(sprintf('1. axial flow, R = %.1f um vs two-step model', R*1e6), e, 1e-2 * (R/20e-6));
end

%% 2. Mixed direction, isotropic sigma: reduces to vUS_3D_quad_vec.m as R shrinks, any orientation
sigma_iso = [130.1474, 130.1474, 44.7070] * 1e-6;
x7b = [8e-3, 5e-3, 0.8, 0.10, 2.5, 0.7, deg2rad(35)];
g_iso6 = vUS_3D_quad_vec(x7b(1:6), tau, k0, sigma_iso, s, w);
for R = [20e-6, 5e-6, 1e-6]
    g_vessel = vUS_3D_quad_RCA_vessel_vec([x7b, R], tau, k0, sigma_iso, s, w, 32);
    e = max(abs(g_iso6 - g_vessel));
    allPass = allPass && e < 1e-2 * (R/20e-6);
    report(sprintf('2. mixed direction, isotropic sigma, R = %.1f um vs single-Gaussian', R*1e6), e, 1e-2 * (R/20e-6));
end

%% 2b. Mixed direction, anisotropic sigma: does NOT reduce to the two-step model as R -> 0
% (documented, not asserted as a pass/fail -- see the header)
x7c = [8e-3, 5e-3, 0.8, 0.10, 2.5, 0.7, deg2rad(20)];
g_twostep_c = vUS_3D_quad_RCA_vec(x7c, tau, k0, sigma, s, w);
g_vessel_tiny = vUS_3D_quad_RCA_vessel_vec([x7c, 5e-8], tau, k0, sigma, s, w, 32); % R -> ~0
e2b = max(abs(g_twostep_c - g_vessel_tiny));
fprintf('%-70s %10.2e  %s\n', '2b. mixed direction, anisotropic sigma, R -> 0 vs two-step model', e2b, '(expected to stay finite, not a failure)');

%% 3. Independent brute-force 3D integration, several orientations
Lam1 = diag([1/sn^2, 1/sw^2, 1/sa^2]);
Lam2 = diag([1/sw^2, 1/sn^2, 1/sa^2]);
R3 = 200e-6;
tau3 = [1, 3, 6]/2500;
k3 = 2.5; a3 = 0.7;
fmax3 = (k3+2)/(k3+2-2*a3^k3);

cases = { ...
    'axis-aligned (v along z)',      0,  20e-3; ...
    'tilted, transverse+axial mix',  40, 20e-3; ...
    'purely transverse (v_z = 0)',   90, 0; ...
};
for ci = 1:size(cases, 1)
    name = cases{ci, 1}; tiltdeg = cases{ci, 2};
    switch tiltdeg
        case 0, vtgp_case = 0; vzgp_case = 20e-3;
        case 90, vtgp_case = 20e-3; vzgp_case = 0;
        otherwise, vtgp_case = 20e-3*sind(tiltdeg); vzgp_case = 20e-3*cosd(tiltdeg);
    end
    phi_case = deg2rad(15);
    x8 = [vtgp_case, vzgp_case, 1, 0, k3, a3, phi_case, R3];
    g_model = vUS_3D_quad_RCA_vessel_vec(x8, tau3, k0, sigma, s, w, 40);

    vgp_vec = [vtgp_case*cos(phi_case); vtgp_case*sin(phi_case); vzgp_case];
    vgp_mag = norm(vgp_vec); n_hat = vgp_vec/max(vgp_mag, eps);
    ref = [0;0;1]; if abs(dot(ref,n_hat)) > 0.9, ref = [1;0;0]; end
    e_xi = ref - dot(ref,n_hat)*n_hat; e_xi = e_xi/norm(e_xi); e_eta = cross(n_hat, e_xi);
    rvec = @(xi,eta,zeta) e_xi*xi + e_eta*eta + n_hat*zeta;
    Amp = @(r, Lam) exp(-0.5*(r(1,:).^2*Lam(1,1) + r(2,:).^2*Lam(2,2) + r(3,:).^2*Lam(3,3)));
    f_of_rho = @(rho) fmax3*(1 - a3^k3*(rho/R3).^k3);
    Zmax = 12*max([sn sw sa]);

    denom_polar = @(rho,th,z) arrayfun(@(rho_,th_,z_) ...
        (Amp(rvec(rho_*cos(th_),rho_*sin(th_),z_),Lam1)+Amp(rvec(rho_*cos(th_),rho_*sin(th_),z_),Lam2))^2 * rho_, rho,th,z);
    Denom_bf = integral3(@(rho,th,z) reshape(denom_polar(rho(:),th(:),z(:)), size(rho)), ...
        0, R3, 0, 2*pi, -Zmax, Zmax, 'AbsTol', 1e-30, 'RelTol', 1e-7, 'Method', 'iterated');

    g_bf = zeros(size(tau3));
    for ti = 1:numel(tau3)
        tv = tau3(ti);
        num_polar = @(rho,th,z) arrayfun(@(rho_,th_,z_) local_num(rho_,th_,z_,rvec,Amp,Lam1,Lam2,f_of_rho,vgp_vec,tv,k0), rho,th,z) .* rho;
        Num_bf = integral3(@(rho,th,z) reshape(num_polar(rho(:),th(:),z(:)), size(rho)), ...
            0, R3, 0, 2*pi, -Zmax, Zmax, 'AbsTol', 1e-30, 'RelTol', 1e-7, 'Method', 'iterated');
        g_bf(ti) = Num_bf/Denom_bf;
    end
    e = max(abs(g_bf(:) - g_model(:)));
    allPass = allPass && e < 1e-6;
    report(sprintf('3. %s vs brute-force integral3', name), e, 1e-6);
end

%% 4. g1(0) = DC + F
x8 = [x7, 200e-6];
g0 = vUS_3D_quad_RCA_vessel_vec(x8, 0, k0, sigma, s, w, 32);
e = abs(g0 - (x8(3)+x8(4)));
allPass = allPass && e < tol;
report('4. g1(0) = DC + F', e, tol);

fprintf('\n');
if allPass
    fprintf('All checks passed.\n');
else
    warning('One or more checks failed -- inspect before use.');
end

function val = local_num(rho, th, z, rvec, Amp, Lam1, Lam2, f_of_rho, vgp_vec, tau, k0)
    r = rvec(rho*cos(th), rho*sin(th), z);
    d = f_of_rho(rho) * norm(vgp_vec) * (vgp_vec/norm(vgp_vec)) * tau;
    A1 = Amp(r, Lam1); A2 = Amp(r, Lam2);
    A1p = Amp(r+d, Lam1); A2p = Amp(r+d, Lam2);
    val = (A1+A2)*(A1p+A2p)*exp(1i*2*k0*dot(d, [0;0;1]));
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
