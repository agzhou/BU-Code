%% Description:
%   Symbolically derives the analytic Jacobian of the two-step RCA-PSF g1 model (fixed-quadrature
%   form, vUS_3D_quad_RCA_xy_vec.m) w.r.t. ALL SEVEN free parameters
%   [v_xgp, v_ygp, v_zgp, F, DC, k, a] (DC a real additive offset), and auto-generates a fast numeric
%   MATLAB function via matlabFunction. Direct RCA-PSF analog of generate_vUS_3D_quad_Jac.m: same
%   change-of-variables (s = (r/R)^2), same F/DC-by-hand shortcut, same one-output-per-quantity
%   generation style (see that file's header for the full rationale).
%
%   Model (see RCA_g1s_derivation.md and vUS_3D_quad_RCA_xy_vec.m):
%       g1(tau) = DC + F * int_0^1 Integrand(tau, f(s)) ds
%       Integrand = EZ * (w1*D1 + w1*D4 + w23*D23),   EZ = exp(i 2 k0 v_zgp tau f - (Mz/4) tau^2 f^2)
%       D1 = exp(-(Mt1/4) tau^2 f^2), D4 = exp(-(Mt4/4) tau^2 f^2), D23 = exp(-(Mt23/4) tau^2 f^2)
%       Mt1 = v_xgp^2/sn^2 + v_ygp^2/sw^2,  Mt4 = v_xgp^2/sw^2 + v_ygp^2/sn^2,  Mt23 = 2(v_xgp^2+v_ygp^2)/S
%       Mz = v_zgp^2/sa^2,  S = sn^2+sw^2,  c = 2 sn sw/S,  w1 = 1/(2(1+c)),  w23 = c/(1+c)
%       f(s) = fmax*(1 - a^k s^(k/2)),  fmax = (k+2)/(k+2-2 a^k)   (same f(s) as the single-Gaussian model)
%
%   F and DC enter linearly and are handled by hand, exactly as in generate_vUS_3D_quad_Jac.m:
%       dg1/dF = int Integrand ds,   dg1/dDC = 1
%   so only Integrand and its 5 nonlinear-parameter derivatives dIntegrand/d[v_xgp,v_ygp,v_zgp,k,a] need
%   to be generated.
%
%   Run this script once (and again any time the model formula changes) to (re)generate:
%       vUS_3D_quad_RCA_xy_complex_Jac_raw.m -- AUTO-GENERATED internal helper: Integrand and its 5
%           derivatives for every (tau, node), do not hand-edit
%   The hand-written vUS_3D_quad_RCA_xy_complex_Jac.m assembles [g1, J] from it, and
%   vUS_3D_quad_RCA_xy_residJac.m plugs it into lsqnonlin.
%
%   The self-checks at the bottom compare (1) the model (vUS_3D_quad_RCA_xy_vec.m) against
%   vUS_3D_num_RCA_vec.m -- the ORIGINAL, independently-written adaptive-integral model this whole file
%   accelerates -- and (2) the generated Jacobian against central finite differences of
%   vUS_3D_quad_RCA_xy_vec.m.
%
% Requires: Symbolic Math Toolbox

%% 1. Build the model symbolically
x     = sym('x', [7, 1], 'real');       % x = [v_xgp; v_ygp; v_zgp; F; DC; k; a]
sigma = sym('sigma', [3, 1], 'real');   % sigma = [sigma_narrow; sigma_wide; sigma_axial]
syms tau k0 real
syms s positive                         % quadrature node in (0, 1)
assumeAlso(x(6) > 0)                    % k > 0
assumeAlso(x(7) > 0)                    % a > 0 (a = 0 is handled numerically in the wrapper)

v_xgp = x(1); v_ygp = x(2); v_zgp = x(3); Fp = x(4); DC = x(5); k = x(6); a = x(7);
sn = sigma(1); sw = sigma(2); sa = sigma(3);

S  = sn^2 + sw^2;
c  = 2*sn*sw/S;
w1 = 1/(2*(1 + c));
w23 = c/(1 + c);

Mt1  = v_xgp^2/sn^2 + v_ygp^2/sw^2;
Mt4  = v_xgp^2/sw^2 + v_ygp^2/sn^2;
Mt23 = 2*(v_xgp^2 + v_ygp^2)/S;
Mz   = v_zgp^2/sa^2;

ak   = a^k;
Dd   = k + 2 - 2*ak;
fmax = (k + 2)/Dd;
f    = fmax*(1 - ak*s^(k/2));

EZ = exp(1i*2*k0*v_zgp*tau*f - (Mz/4)*tau^2*f^2);
D1 = exp(-(Mt1/4)*tau^2*f^2);
D4 = exp(-(Mt4/4)*tau^2*f^2);
D23 = exp(-(Mt23/4)*tau^2*f^2);
Integrand = EZ*(w1*D1 + w1*D4 + w23*D23);   % one quadrature node's contribution to int_0^1 (...) ds

%% 2. Symbolic Jacobian of the integrand w.r.t. the nonlinear parameters
nl = [v_xgp, v_ygp, v_zgp, k, a];
dI = jacobian(Integrand, nl);   % 1x5, complex-valued

%% 3. Generate the raw numeric code -- Integrand + 5 SEPARATE derivative outputs.
% Vars = {x, tau, k0, sigma, s}: call it with tau [nTau, 1] and s [1, N] and every output comes back
% [nTau, N] via implicit expansion (all generated operations are elementwise).
vUSDir = fileparts(mfilename('fullpath'));
rawFile = fullfile(vUSDir, 'vUS_3D_quad_RCA_xy_complex_Jac_raw.m');
matlabFunction(Integrand, dI(1), dI(2), dI(3), dI(4), dI(5), ...
    'File', rawFile, ...
    'Vars', {x, tau, k0, sigma, s}, ...
    'Outputs', {'I', 'dI_dvx', 'dI_dvy', 'dI_dvz', 'dI_dk', 'dI_da'}, ...
    'Optimize', true);

txt = fileread(rawFile);
banner = sprintf(['%%%% AUTO-GENERATED (internal helper) by generate_vUS_3D_quad_RCA_xy_Jac.m on %s.\n', ...
    '%%%% Do not call directly or hand-edit -- use vUS_3D_quad_RCA_xy_complex_Jac.m,\n', ...
    '%%%% and re-run the generator script to update.\n\n'], char(datetime('now')));
fid = fopen(rawFile, 'w');
fwrite(fid, [banner, txt]);
fclose(fid);
rehash path
fprintf('Generated %s\n', rawFile);

%% 4. Self-check (1): the fast quadrature model vs. the original adaptive-integral model
% (vUS_3D_num_RCA_vec.m) it accelerates -- both were written independently.
sigma_test = [57.6, 286.5, 50.1] * 1e-6;               % field sigmas, RC15gV 5 x 2 angles -6 to 6 deg (fit_RCA_PSF_v2.m)
k0_test    = 2*pi / (1540/13.6e6);                     % [rad/m]
tau_test   = (1:40)' / 2500;                           % [s], 2.5 kHz frame rate
[s_test, w_test] = gaussLegendre01(48);

fprintf('\n%-40s %14s\n', 'model vs. vUS_3D_num_RCA_vec.m', 'max abs. err');
for cs = [2.5 0.7; 2 0.5; 3 0.9; 2.5 1; 2.5 0.05]'
    kk = cs(1); aa = cs(2);
    xx = [8e-3; -5e-3; 8e-3; 1; 0; kk; aa];
    gnum = vUS_3D_num_RCA_vec(xx.', tau_test, k0_test, sigma_test);           % x row for vUS_3D_num_RCA_vec.m
    gq   = vUS_3D_quad_RCA_xy_vec(xx, tau_test, k0_test, sigma_test, s_test, w_test);
    fprintf('k = %-4.2f  a = %-5.2f %30.3e\n', kk, aa, max(abs(gnum(:) - gq(:))));
end

%% 5. Self-check (2): analytic Jacobian vs. central finite differences of
% the independently-written vUS_3D_quad_RCA_xy_vec.m. Points cover typical flow, the k and a bounds,
% near-stagnant flow (small M), small a, and vx/vy swapped or signed differently (checking the
% sign/swap invariance doesn't hide a Jacobian bug).
testPoints = { ...
    [ 8e-3; -5e-3; -8e-3;  0.8; 0.10; 2.5; 0.70], 'typical flow'; ...
    [ 5e-3;  8e-3;  5e-3;  0.9; 0.00; 2.0; 1.00], 'k = 2, a = 1 (k-only limit)'; ...
    [12e-3; -3e-3; -3e-3;  0.6; 0.20; 3.0; 0.40], 'k = 3, moderate a'; ...
    [ 1e-4;  1e-4;  1e-4;  0.5; 0.05; 2.5; 0.90], 'near-stagnant (small M)'; ...
    [ 8e-3; -5e-3; -8e-3;  0.8; 0.10; 2.5; 0.05], 'small a'; ...
    [-8e-3;  5e-3;  15e-3; 0.9; 0.00; 3.0; 0.90], 'fast axial, negative vxgp'};

hList = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
fprintf('\n%-32s %14s %14s\n', 'test point', 'max rel. err', 'max abs. err');
allPass = true;
for i = 1:size(testPoints, 1)
    x0 = testPoints{i, 1};
    [~, J_analytic] = vUS_3D_quad_RCA_xy_complex_Jac(x0, tau_test, k0_test, sigma_test, s_test, w_test);

    bestRelErr = Inf(1, 7);
    bestAbsErr = Inf(1, 7);
    for h = hList
        for p = 1:7
            dx = zeros(7, 1);
            step = h * max(abs(x0(p)), 1e-6);
            dx(p) = step;

            g1_plus  = vUS_3D_quad_RCA_xy_vec(x0 + dx, tau_test, k0_test, sigma_test, s_test, w_test);
            g1_minus = vUS_3D_quad_RCA_xy_vec(x0 - dx, tau_test, k0_test, sigma_test, s_test, w_test);
            fd_col = (g1_plus - g1_minus) / (2*step);
            fd_col = [real(fd_col); imag(fd_col)];

            Jc = [real(J_analytic(:,p)); imag(J_analytic(:,p))];

            colScale = max(abs(fd_col));
            sig = abs(fd_col) > 1e-2 * colScale;
            relErr = abs(Jc(sig) - fd_col(sig)) ./ abs(fd_col(sig));
            bestRelErr(p) = min(bestRelErr(p), max(relErr));
            bestAbsErr(p) = min(bestAbsErr(p), max(abs(Jc - fd_col)) / max(colScale, eps));
        end
    end

    mRel = max(bestRelErr); mAbs = max(bestAbsErr);
    fprintf('%-32s %14.3e %14.3e\n', testPoints{i, 2}, mRel, mAbs);
    allPass = allPass && (mRel < 1e-4);
end

%% 6. Self-check (3): isotropic reduction (sigma_narrow == sigma_wide) matches the ALREADY-VALIDATED
% single-Gaussian Jacobian (vUS_3D_quad_complex_Jac.m) exactly -- ties the new 7-parameter Jacobian to
% the existing, independently-checked 6-parameter one.
sigma_iso = [130.1474, 130.1474, 44.7070] * 1e-6;
x6 = [8e-3, -8e-3, 0.8, 0.10, 2.5, 0.7];                          % [v_tgp, v_zgp, F, DC, k, a]
x7 = [x6(1), 0, x6(2), x6(3), x6(4), x6(5), x6(6)];               % same physical point: v_xgp = v_tgp, v_ygp = 0
[g6, J6] = vUS_3D_quad_complex_Jac(x6, tau_test, k0_test, sigma_iso, s_test, w_test);
[g7, J7] = vUS_3D_quad_RCA_xy_complex_Jac(x7, tau_test, k0_test, sigma_iso, s_test, w_test);
dG = max(abs(g6 - g7));
dJ = max(abs([J6(:,1), J6(:,2:6)] - [J7(:,1), J7(:,3:7)]), [], 'all');   % v_xgp<->v_tgp, v_zgp,F,DC,k,a line up; v_ygp has no counterpart
check6 = dG < 1e-12 && dJ < 1e-9;
if check6, pf6 = 'PASS'; else, pf6 = 'FAIL'; end
fprintf('\n6. isotropic reduction vs vUS_3D_quad_complex_Jac.m: max|dg1| = %.3e, max|dJ| = %.3e -> %s\n', dG, dJ, pf6);
allPass = allPass && check6;

if allPass
    fprintf('\nAll test points passed (analytic vs. finite-difference of vUS_3D_quad_RCA_xy_vec.m agree, and the isotropic reduction matches vUS_3D_quad_complex_Jac.m).\n');
else
    warning('One or more test points show large analytic-vs-finite-difference disagreement -- inspect before use.');
end
