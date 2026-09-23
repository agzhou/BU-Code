%% Description:
%   Symbolically derives the analytic Jacobian of the CORRECTED
%   Tangelder-profile g1 model (fixed-quadrature form, vUS_3D_quad_vec.m)
%   w.r.t. ALL SIX free parameters [v_tgp, v_zgp, F, DC, k, a] (DC a real
%   additive offset), and auto-generates a fast numeric MATLAB function via
%   matlabFunction.
%
%   Model (see the derivation write-up):
%       g1(tau) = DC + F * int_0^1 exp(phi) ds,
%       phi     = -M tau^2 f^2/4 + i 2 k0 v_zgp tau f
%       f(s)    = fmax * (1 - a^k s^(k/2)),   fmax = (k+2)/(k+2-2 a^k)
%       M       = v_tgp^2/sigma_t^2 + v_zgp^2/sigma_z^2
%
%   The integration variable s = (r/R)^2 is parameter-free (that is what
%   the change of variables in the derivation buys), so the Gauss-Legendre
%   nodes/weights do not depend on any fit parameter. d/dx therefore
%   commutes with the quadrature: the Jacobian of the discretised model is
%   exactly the quadrature of the Jacobian of the integrand, and only
%   the integrand E = exp(phi) and its 4 nonlinear-parameter derivatives
%   dE/d[v_tgp, v_zgp, k, a] need to be generated. F and DC enter linearly
%   and are handled in the hand-written wrapper:
%       dg1/dF = int E ds,   dg1/dDC = 1
%
%   Run this script once (and again any time the model formula changes) to
%   (re)generate:
%       vUS_3D_quad_complex_Jac_raw.m -- AUTO-GENERATED internal helper: E and
%           dE/d{v_tgp, v_zgp, k, a} for every (tau, node), do not hand-edit
%   The hand-written vUS_3D_quad_complex_Jac.m assembles [g1, J] from it, and
%   vUS_3D_quad_residJac.m plugs it into lsqnonlin. Same design as
%   generate_vUS_3D_erf_Jac.m: one separate matlabFunction output per
%   quantity (no internal concatenation).
%
%   The self-checks at the bottom compare (1) the model against a direct
%   integral over the vessel cross-section (checks the formula itself,
%   independent of the quadrature) and (2) the generated Jacobian against
%   central finite differences of vUS_3D_quad_vec.m.
%
% Requires: Symbolic Math Toolbox

%% 1. Build the model symbolically
x     = sym('x', [6, 1], 'real');       % x = [v_tgp; v_zgp; F; DC; k; a]
sigma = sym('sigma', [3, 1], 'real');   % sigma = [sigma_x; sigma_y; sigma_z]
syms tau k0 real
syms s positive                         % quadrature node in (0, 1)
assumeAlso(x(5) > 0)                    % k > 0
assumeAlso(x(6) > 0)                    % a > 0 (a = 0 is handled numerically in the wrapper)

v_tgp = x(1); v_zgp = x(2); Fp = x(3); DC = x(4); k = x(5); a = x(6);
sx = sigma(1); sz = sigma(3);           % sigma_x = sigma_y = sigma_t

M    = v_tgp^2/sx^2 + v_zgp^2/sz^2;
ak   = a^k;
D    = k + 2 - 2*ak;
fmax = (k + 2)/D;
f    = fmax*(1 - ak*s^(k/2));
phi  = -M*tau^2/4*f^2 + 2*1i*k0*v_zgp*tau*f;
E    = exp(phi);                        % one quadrature node's contribution to int_0^1 (...) ds

%% 2. Symbolic Jacobian of the integrand w.r.t. the nonlinear parameters
nl = [v_tgp, v_zgp, k, a];
dE = jacobian(E, nl);                   % 1x4, complex-valued

% Optional: print the building blocks (chain-rule pieces) in closed form
printBlocks = false;
if printBlocks
    disp('dphi/dv_tgp ='); disp(simplify(diff(phi, v_tgp)))
    disp('dphi/dv_zgp ='); disp(simplify(diff(phi, v_zgp)))
    disp('dfmax/dk    ='); disp(simplify(diff(fmax, k)))
    disp('dfmax/da    ='); disp(simplify(diff(fmax, a)))
    disp('df/dk       ='); disp(simplify(diff(f, k)))
    disp('df/da       ='); disp(simplify(diff(f, a)))
end

%% 3. Generate the raw numeric code -- E + 4 SEPARATE derivative outputs.
% Vars = {x, tau, k0, sigma, s}: call it with tau [nTau, 1] and s [1, N] and
% every output comes back [nTau, N] via implicit expansion (all generated
% operations are elementwise).
vUSDir = fileparts(mfilename('fullpath'));
rawFile = fullfile(vUSDir, 'vUS_3D_quad_complex_Jac_raw.m');
matlabFunction(E, dE(1), dE(2), dE(3), dE(4), ...
    'File', rawFile, ...
    'Vars', {x, tau, k0, sigma, s}, ...
    'Outputs', {'E', 'dE_dvt', 'dE_dvz', 'dE_dk', 'dE_da'}, ...
    'Optimize', true);

txt = fileread(rawFile);
banner = sprintf(['%%%% AUTO-GENERATED (internal helper) by generate_vUS_3D_quad_Jac.m on %s.\n', ...
    '%%%% Do not call directly or hand-edit -- use vUS_3D_quad_complex_Jac.m,\n', ...
    '%%%% and re-run the generator script to update.\n\n'], char(datetime('now')));
fid = fopen(rawFile, 'w');
fwrite(fid, [banner, txt]);
fclose(fid);
rehash path
fprintf('Generated %s\n', rawFile);

%% 4. Self-check (1): the model formula vs. a direct integral over the disc
% Ground truth needs no pdf algebra and no quadrature in s:
%   v(rho) = vmax (1 - (a rho)^k),  vgp = vmax (1 - 2 a^k/(k+2)),  f = v/vgp,
%   g1 = int_0^1 2 rho exp(phi(f(rho))) d rho     (rho = r/R)
sigma_test = [368.1124; 368.1124; 126.4505] * 1e-6;   % RC15gV PSF (5 x 2 angles, -6 to 6 deg)
k0_test    = 2*pi / (1540/13.6e6);                    % [rad/m]
tau_test   = (1:40)' / 2500;                          % [s], 2.5 kHz frame rate, tau1..tau40
[s_test, w_test] = gaussLegendre01(48);

fprintf('\n%-34s %14s\n', 'model vs. direct disc integral', 'max abs. err');
for cs = [2.5 0.7; 2 0.5; 3 0.9; 2.5 1; 2.5 0.05]'
    kk = cs(1); aa = cs(2);
    xx = [8e-3; -8e-3; 1; 0; kk; aa];
    Mt = xx(1)^2/sigma_test(1)^2 + xx(2)^2/sigma_test(3)^2;
    fr = @(rho) (1 - (aa*rho).^kk) ./ (1 - 2*aa^kk/(kk+2));
    gd = integral(@(rho) 2*rho .* exp(-Mt/4*tau_test.^2.*fr(rho).^2 + 2i*k0_test*xx(2)*tau_test.*fr(rho)), ...
        0, 1, 'ArrayValued', true, 'RelTol', 1e-12, 'AbsTol', 1e-14);
    gq = vUS_3D_quad_vec(xx, tau_test, k0_test, sigma_test, s_test, w_test);
    fprintf('k = %-4.2f  a = %-5.2f %24.3e\n', kk, aa, max(abs(gd - gq)));
end

%% 5. Self-check (2): analytic Jacobian vs. central finite differences of
% the independently-written vUS_3D_quad_vec.m. Points cover typical flow, the
% k and a bounds, near-stagnant flow (small M) and small a.
testPoints = { ...
    [ 8e-3; -8e-3;  0.8; 0.10; 2.5; 0.70], 'typical flow'; ...
    [ 3e-3;  5e-3;  0.9; 0.00; 2.0; 1.00], 'k = 2, a = 1 (k-only limit)'; ...
    [12e-3; -3e-3;  0.6; 0.20; 3.0; 0.40], 'k = 3, moderate a'; ...
    [ 1e-4;  1e-4;  0.5; 0.05; 2.5; 0.90], 'near-stagnant (small M)'; ...
    [ 8e-3; -8e-3;  0.8; 0.10; 2.5; 0.05], 'small a'; ...
    [ 5e-3;  15e-3; 0.9; 0.00; 3.0; 0.90], 'fast axial, k = 3'};

hList = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
fprintf('\n%-30s %14s %14s\n', 'test point', 'max rel. err', 'max abs. err');
allPass = true;
for i = 1:size(testPoints, 1)
    x0 = testPoints{i, 1};
    [~, J_analytic] = vUS_3D_quad_complex_Jac(x0, tau_test, k0_test, sigma_test, s_test, w_test);

    bestRelErr = Inf(1, 6);
    bestAbsErr = Inf(1, 6);
    for h = hList
        for p = 1:6
            dx = zeros(6, 1);
            step = h * max(abs(x0(p)), 1e-6);
            dx(p) = step;
            g1_plus  = vUS_3D_quad_vec(x0 + dx, tau_test, k0_test, sigma_test, s_test, w_test);
            g1_minus = vUS_3D_quad_vec(x0 - dx, tau_test, k0_test, sigma_test, s_test, w_test);
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
    fprintf('%-30s %14.3e %14.3e\n', testPoints{i, 2}, mRel, mAbs);
    allPass = allPass && (mRel < 1e-4);
end

if allPass
    fprintf('\nAll test points passed (analytic vs. finite-difference of vUS_3D_quad_vec.m agree).\n');
else
    warning('One or more test points show large analytic-vs-finite-difference disagreement -- inspect before use.');
end
