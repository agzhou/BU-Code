%% Description:
%   Validation of vUS_3D_quad_RCA_vec.m (g1 model for the RCA PSF: sum of two
%   "orthogonal" lateral Gaussians x an axial Gaussian). Checks:
%     1. sigma_narrow == sigma_wide reduces exactly to vUS_3D_quad_vec.m (any phi)
%     2. g1(0) = DC + F
%     3. Plug flow (a = 0): the model equals g1,s, checked against brute-force
%        numerical integration of the sIQ autocorrelation (i = j terms, i.e.,
%        the point where the derivation assumes the phase-scrambled cross
%        terms vanish), for several velocities and directions
%     4. General k, a: quadrature over s vs. a direct integral over the vessel
%        cross-section (rho = r/R)
%     5. k = 2, a = 1: equals the weighted sum of three vUS_3D_erf.m
%        evaluations (the old closed-form model), one per Gaussian term
%     6. Symmetries: phi -> -phi, phi -> pi/2 - phi
%   Also prints how much g1 depends on the transverse flow direction.

clearvars

%% Add paths
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
AllenCodeDir = join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\";
addpath(genpath(fullfile(AllenCodeDir + "Processing\")))
addpath(genpath(fullfile(AllenCodeDir + "ErrorFunction\")))

%% Test parameters (RC15gV, 5 x 2 angles from -6 to 6 deg; sigmas from fit_RCA_PSF_v2.m)
sigma   = [57.6, 286.5, 50.1] * 1e-6;       % [sigma_narrow, sigma_wide, sigma_axial], field sigmas [m]
k0      = 2*pi / (1540/13.6e6);             % [rad/m]
tau     = (0:40)' / 2500;                   % [s], 2.5 kHz frame rate
[s, w]  = gaussLegendre01(48);
tol     = 1e-9;
allPass = true;
report  = @(name, err, tol) fprintf('%-62s %10.2e  %s\n', name, err, ternary(err < tol, 'pass', 'FAIL'));

x0 = [8e-3, 5e-3, 0.8, 0.10, 2.5, 0.7, deg2rad(20)];   % [v_tgp, v_zgp, F, DC, k, a, phi]

%% 1. sigma_narrow == sigma_wide reduces to the single-Gaussian model
sigma_iso = [130.1474, 130.1474, 44.7070] * 1e-6;
for phi = deg2rad([0, 20, 45, 70])
    x = x0; x(7) = phi;
    g_new = vUS_3D_quad_RCA_vec(x, tau, k0, sigma_iso, s, w);
    g_old = vUS_3D_quad_vec(x(1:6), tau, k0, sigma_iso, s, w);
    e = max(abs(g_new - g_old));
    allPass = allPass && e < tol;
    report(sprintf('1. sn == sw vs vUS_3D_quad_vec, phi = %2.0f deg', rad2deg(phi)), e, tol);
end

%% 2. g1(0) = DC + F
g0 = vUS_3D_quad_RCA_vec(x0, 0, k0, sigma, s, w);
e = abs(g0 - (x0(4) + x0(3)));
allPass = allPass && e < tol;
report('2. g1(0) = DC + F', e, tol);

%% 3. Plug flow (a = 0) vs brute-force numerical autocorrelation of the sIQ model
G  = @(u, sg) exp(-u.^2./(2*sg.^2));
PSFlat = @(X, Y) G(X, sigma(1)).*G(Y, sigma(2)) + G(X, sigma(2)).*G(Y, sigma(1));
Lx = 9*sigma(2); nx = 1201;
xg = linspace(-Lx, Lx, nx); dxg = xg(2) - xg(1);
[Xg, Yg] = ndgrid(xg, xg);
P0 = PSFlat(Xg, Yg);
den_lat = sum(P0(:).^2) * dxg^2;
zg = linspace(-8*sigma(3), 8*sigma(3) + 4e-4, 3001); dzg = zg(2) - zg(1);
den_ax = sum(G(zg, sigma(3)).^2) * dzg;

tau3 = (0:4:40)' / 2500;
vcases = [ 8e-3, 5e-3, deg2rad(0); ...
           8e-3, 5e-3, deg2rad(20); ...
           8e-3, -3e-3, deg2rad(45); ...
           15e-3, 0, deg2rad(70) ];          % [v_tgp, v_zgp, phi]
for i = 1:size(vcases, 1)
    vt = vcases(i,1); vz = vcases(i,2); phi = vcases(i,3);
    d = [vt*cos(phi), vt*sin(phi), vz] .* 1;
    g_bf = zeros(size(tau3));
    for m = 1:numel(tau3)
        dm = d * tau3(m);
        Pt = PSFlat(Xg + dm(1), Yg + dm(2));
        num_lat = sum(P0(:).*Pt(:)) * dxg^2;
        num_ax  = sum(G(zg, sigma(3)).*G(zg + dm(3), sigma(3))) * dzg;
        g_bf(m) = num_lat*num_ax/(den_lat*den_ax) * exp(1i*2*k0*dm(3));
    end
    g_mod = vUS_3D_quad_RCA_vec([vt, vz, 1, 0, 2.5, 0, phi], tau3, k0, sigma, s, w);
    e = max(abs(g_mod - g_bf));
    allPass = allPass && e < 1e-6;
    report(sprintf('3. a = 0 vs brute force, v_t=%.0f mm/s, v_z=%.0f, phi=%2.0f deg', vt*1e3, vz*1e3, rad2deg(phi)), e, 1e-6);
end

%% 4. General k, a: quadrature over s vs. direct integral over rho = r/R
for cs = [2.5 0.7; 2 0.5; 3 0.9; 2.5 1; 2.5 0.05]'
    kk = cs(1); aa = cs(2);
    x = x0; x(5) = kk; x(6) = aa;
    fr = @(rho) (1 - (aa*rho).^kk) ./ (1 - 2*aa^kk/(kk+2));
    gd = integral(@(rho) 2*rho .* g1s_formula(tau, x(1), x(2), x(7), fr(rho), k0, sigma), 0, 1, ...
        'ArrayValued', true, 'RelTol', 1e-12, 'AbsTol', 1e-14);
    gq = vUS_3D_quad_RCA_vec(x, tau, k0, sigma, s, w);
    e = max(abs(x(4) + x(3)*gd - gq));
    allPass = allPass && e < 1e-8;
    report(sprintf('4. quadrature vs disc integral, k = %.2f, a = %.2f', kk, aa), e, 1e-8);
end

%% 5. k = 2, a = 1: weighted sum of three vUS_3D_erf.m (old closed form) evaluations
tau5 = tau(2:end);                                   % vUS_3D_erf divides by tau
S = sigma(1)^2 + sigma(2)^2; c = 2*sigma(1)*sigma(2)/S; sigma_e = sqrt(S/2);
for phi = deg2rad([0, 25, 45])
    x = [x0(1:4), 2, 1, phi]; vt = x(1); vz = x(2);
    vx = vt*cos(phi); vy = vt*sin(phi);
    G1  = vUS_3D_erf(tau5, k0, [sigma(1), sigma(2), sigma(3)], vx, vy, vz);
    G4  = vUS_3D_erf(tau5, k0, [sigma(2), sigma(1), sigma(3)], vx, vy, vz);
    G23 = vUS_3D_erf(tau5, k0, [sigma_e, sigma_e, sigma(3)], vx, vy, vz);
    g_erf = x(4) + x(3)*((G1 + G4)/2 + c*G23)/(1 + c);
    g_q   = vUS_3D_quad_RCA_vec(x, tau5, k0, sigma, s, w);
    e = max(abs(g_erf(:) - g_q(:)));
    allPass = allPass && e < 1e-8;
    report(sprintf('5. k=2, a=1 vs weighted sum of 3 erf models, phi = %2.0f deg', rad2deg(phi)), e, 1e-8);
end

%% 6. Symmetries in phi
xa = x0; xa(7) = deg2rad(20);
gA = vUS_3D_quad_RCA_vec(xa, tau, k0, sigma, s, w);
xb = xa; xb(7) = -xa(7);                gB = vUS_3D_quad_RCA_vec(xb, tau, k0, sigma, s, w);
xc = xa; xc(7) = pi/2 - xa(7);          gC = vUS_3D_quad_RCA_vec(xc, tau, k0, sigma, s, w);
e = max([max(abs(gA - gB)), max(abs(gA - gC))]);
allPass = allPass && e < tol;
report('6. g1(phi) = g1(-phi) = g1(pi/2 - phi)', e, tol);

fprintf('\n');
if allPass
    fprintf('All checks passed.\n');
else
    warning('One or more checks failed -- inspect before use.');
end

%% Info: dependence of g1 on the transverse flow direction (plug flow, no axial flow)
fprintf('\n|g1| for plug flow (a = 0), v_z = 0, v_t = 8 mm/s, vs. transverse flow direction\n');
fprintf('%-10s %10s %10s %10s %10s   |  single-Gaussian (sigma = %.0f um)\n', 'tau [ms]', 'phi = 0', 'phi = 22.5', 'phi = 45', 'max diff', 130.1474);
tauI = [1 2 4 6 8 12] * 1e-3;
for t = tauI
    gg = arrayfun(@(p) abs(vUS_3D_quad_RCA_vec([8e-3, 0, 1, 0, 2.5, 0, deg2rad(p)], t, k0, sigma, s, w)), [0 22.5 45]);
    g_sg = abs(vUS_3D_quad_vec([8e-3, 0, 1, 0, 2.5, 0], t, k0, sigma_iso, s, w));
    fprintf('%-10.1f %10.4f %10.4f %10.4f %10.4f   |  %10.4f\n', t*1e3, gg, max(gg) - min(gg), g_sg);
end
% Expanding g1,s to first order in tau^2 gives an isotropic short-lag decay, 1 - v_t^2 tau^2/(4 sigma_eq^2), with
%   1/sigma_eq^2 = (1/sn^2 + 1/sw^2)/(2(1 + c)) + 2c/((1 + c) S)
sigma_eq = 1/sqrt((1/sigma(1)^2 + 1/sigma(2)^2)/(2*(1 + c)) + 2*c/((1 + c)*S));
fprintf('\nShort-lag (tau -> 0) equivalent single-Gaussian lateral sigma: %.1f um (independent of phi)\n', sigma_eq*1e6);

%% Local functions
function g = g1s_formula(tau, vt, vz, phi, f, k0, sigma)
    % g1,s for group velocity f*[vt cos(phi), vt sin(phi), vz], written directly from the derivation
    % (independent of vUS_3D_quad_RCA_vec.m's M_j form). tau: [nTau, 1]; f: scalar or [1, N]
    sn = sigma(1); sw = sigma(2); sa = sigma(3);
    vx = f.*vt.*cos(phi); vy = f.*vt.*sin(phi); vzf = f.*vz;
    S = sn^2 + sw^2; c = 2*sn*sw/S;
    E1  = exp(-tau.^2.*(vx.^2/(4*sn^2) + vy.^2/(4*sw^2)));
    E4  = exp(-tau.^2.*(vx.^2/(4*sw^2) + vy.^2/(4*sn^2)));
    E23 = exp(-tau.^2.*(vx.^2 + vy.^2)/(2*S));
    g = ((E1 + E4)/2 + c*E23)/(1 + c) .* exp(-(vzf.*tau).^2/(4*sa^2) + 1i*2*k0*vzf.*tau);
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
