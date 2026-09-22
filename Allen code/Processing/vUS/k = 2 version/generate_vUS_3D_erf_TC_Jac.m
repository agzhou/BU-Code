%% Description:
%   TC (Transverse-Combined) analog of generate_vUS_3D_erf_Jac.m:
%   symbolically derives the analytic Jacobian of the new (Poiseuille /
%   uniform-velocity-distribution) g1 model -- vUS_3D_erf_TC.m /
%   vUS_3D_erf_TC_vec_split.m -- w.r.t. its four free parameters
%   [v_tgp, v_zgp, F, DC] (DC a real additive offset), and auto-generates
%   a fast numeric MATLAB function via matlabFunction.
%
%   Run this script once (and again any time the model formula itself
%   changes) to (re)generate:
%       vUS_3D_erf_TC_complex_Jac_raw.m  -- AUTO-GENERATED internal
%           helper, do not call directly or hand-edit
%       vUS_3D_erf_TC_complex_Jac.m      -- AUTO-GENERATED public [g1,J]
%           interface, do not hand-edit
%   Then use the small hand-written wrappers
%   vUS_3D_erf_TC_vec_split_Jac.m and vUS_3D_erf_TC_residJac.m to plug
%   into lsqnonlin.
%
%   Same two-file design as generate_vUS_2D_erf_Jac.m /
%   generate_vUS_3D_erf_Jac.m, for the same reason: dg1/dDC is an exact
%   constant (1), which matlabFunction correctly simplifies away -- but
%   that then makes a single J=[col1,...,col4] concatenation crash, since
%   a bare scalar doesn't auto-broadcast against the other tau-length
%   columns under horzcat. Forcing a fake "+0*tau" dependency to dodge
%   this gets optimized back out by matlabFunction's own CSE pass (tested
%   and confirmed for the 2D case on R2025b / Symbolic Math Toolbox
%   25.2) -- so each column is generated as its own separate output
%   instead (no internal concatenation at all), and the public wrapper
%   below assembles + explicitly broadcasts them. 'Vectorize' is
%   intentionally omitted for the same reason as the other generators:
%   not a recognized matlabFunction option in this Symbolic Math Toolbox
%   version; elementwise (.^, .*, ./) code generation is the default now.
%
%   Mathematically this is the SAME two-term M = a^2/sigma_a^2 +
%   v_zgp^2/sigma_z^2 form as the 2D model (vUS_2D_erf.m) -- v_tgp plays
%   exactly the role v_xgp played there, just relabeled/reinterpreted as
%   the combined transverse magnitude sqrt(v_xgp^2+v_ygp^2) rather than a
%   single Cartesian component. It inherits the same fundamental
%   identifiability issue: v_tgp enters only through M (no phase term),
%   so this Jacobian's dg1/dv_tgp column will show the same shallow,
%   sign-blind (here: magnitude-blind, since v_tgp>=0 by construction)
%   sensitivity as v_xgp did in 2D. sigma_y (sigma(2)) never appears in
%   the model itself (only used for the sigma_x==sigma_y validation
%   check in the hand-written wrapper), so it won't appear in the
%   generated code either.
%
%   This script self-checks the generated Jacobian against finite
%   differences of the independently-written vUS_3D_erf_TC_vec_split.m
%   at a few test points, including a near-stagnant-flow case (small M).
%
% Requires: Symbolic Math Toolbox

%% 1. Build the model symbolically
x     = sym('x', [4, 1], 'real');       % x = [v_tgp; v_zgp; F; DC]
sigma = sym('sigma', [3, 1], 'real');   % sigma = [sigma_x; sigma_y; sigma_z]
syms tau k0 real

v_tgp = x(1); v_zgp = x(2); Fp = x(3); DC = x(4);
sx = sigma(1); sz = sigma(3); % sigma(2) unused: sigma_x==sigma_y is enforced, not modeled separately

M = v_tgp^2/sx^2 + v_zgp^2/sz^2;

A = sqrt(M)*tau - 2*1i*k0*v_zgp/sqrt(M);
B = -2*1i*k0*v_zgp/sqrt(M);

g1 = DC + Fp/2 * sqrt(sym(pi)/M)/tau * exp(-4*k0^2*v_zgp^2/M) * (erf(A) - erf(B));

%% 2. Symbolic Jacobian w.r.t. [v_tgp, v_zgp, F, DC]
J = jacobian(g1, x);   % 1x4, complex-valued

%% 3. Generate the raw numeric code -- g1 + 4 SEPARATE Jacobian-column
% outputs (see header for why this avoids matlabFunction's internal
% concatenation crash).
vUSDir = fileparts(mfilename('fullpath'));
rawFile = fullfile(vUSDir, 'vUS_3D_erf_TC_complex_Jac_raw.m');
matlabFunction(g1, J(1), J(2), J(3), J(4), ...
    'File', rawFile, ...
    'Vars', {x, tau, k0, sigma}, ...
    'Outputs', {'g1', 'dg1_dvtgp', 'dg1_dvzgp', 'dg1_dF', 'dg1_dDC'}, ...
    'Optimize', true);

% Swap the symbolic 'erf' calls for the codebase's complex-safe 'erfz'
% (see generate_vUS_2D_erf_Jac.m's header for the full rationale).
txt = fileread(rawFile);
txt = regexprep(txt, '(?<![A-Za-z0-9_])erf\(', 'erfz(');
banner = sprintf(['%%%% AUTO-GENERATED (internal helper) by generate_vUS_3D_erf_TC_Jac.m on %s.\n', ...
    '%%%% Do not call directly or hand-edit -- use vUS_3D_erf_TC_complex_Jac.m,\n', ...
    '%%%% and re-run the generator script to update both.\n\n'], char(datetime('now')));
fid = fopen(rawFile, 'w');
fwrite(fid, [banner, txt]);
fclose(fid);
fprintf('Generated %s\n', rawFile);

%% 4. Write the thin public wrapper: standard [g1, J] interface (matching
% vUS_3D_erf_TC_vec_split_Jac.m / vUS_3D_erf_TC_residJac.m's
% expectations), with each Jacobian column explicitly broadcast to
% size(tau) before concatenation. NOTE: this wrapper's sigma is the
% Jacobian's own [sigma_x,sigma_y,sigma_z] -- it does NOT re-validate
% sigma_x==sigma_y (that check lives in vUS_3D_erf_TC_vec_split.m /
% vUS_3D_erf_TC_vec_split_Jac.m, closer to where callers actually enter).
outFile = fullfile(vUSDir, 'vUS_3D_erf_TC_complex_Jac.m');
pubCode = sprintf([ ...
    '%%%% AUTO-GENERATED by generate_vUS_3D_erf_TC_Jac.m on %s.\n', ...
    '%%%% Do not hand-edit -- re-run the generator script instead.\n', ...
    '%%%% Public [g1, J] wrapper around vUS_3D_erf_TC_complex_Jac_raw.m.\n\n', ...
    'function [g1, J] = vUS_3D_erf_TC_complex_Jac(x, tau, k0, sigma)\n', ...
    '    tau = tau(:);\n', ...
    '    if nargout > 1\n', ...
    '        [g1, dvt, dvz, dF, dDC] = vUS_3D_erf_TC_complex_Jac_raw(x, tau, k0, sigma);\n', ...
    '        onesTau = ones(size(tau));\n', ...
    '        J = [dvt.*onesTau, dvz.*onesTau, dF.*onesTau, dDC.*onesTau];\n', ...
    '    else\n', ...
    '        g1 = vUS_3D_erf_TC_complex_Jac_raw(x, tau, k0, sigma);\n', ...
    '    end\n', ...
    '    g1 = g1(:) .* ones(size(tau));\n', ...
    'end\n'], char(datetime('now')));
fid = fopen(outFile, 'w');
fwrite(fid, pubCode);
fclose(fid);
rehash path
fprintf('Generated %s\n', outFile);

%% 5. Self-check: compare the analytic Jacobian against finite differences
% of the independently-written vUS_3D_erf_TC_vec_split.m. Test points use
% the RC15gV matrix-probe PSF from test_vUS_3D_erf.m, with sigma_x=sigma_y
% as this model requires.
sigma_test = [429; 429; 126] * 1e-6;       % [m], RC15gV probe (column --
                                             % the symbolic-generated
                                             % function indexes in4(2,:)/
                                             % in4(3,:), which requires a
                                             % column input)
k0_test    = 2*pi / (1540/15.625e6);       % [rad/m]
tau_test   = (1:40)' / 5000;               % [s], 5 kHz frame rate, tau1..tau40

testPoints = { ...
    [ 6e-3;   8e-3; 0.7; 0.05], 'typical flow'; ...
    [ 3e-4;   3e-4; 0.5; 0.05], 'near-stagnant (small M)'; ...
    [15e-3;  -5e-3; 0.9; 0.00], 'fast transverse, negative vzgp'};
% NOTE: the near-stagnant point is 3e-4, not 1e-4 as in the 2D/3D
% generators. At exactly [1e-4,1e-4] with this probe's sigma (RC15gV,
% sigma_z=126um -- larger than the 2D L22-14v case's 74um), |exp(-B^2)|
% reaches ~1e103 (vs. ~1e15 for the analogous 2D test point), landing
% squarely on this probe's worst-conditioned corner: relative error
% there plateaus around 1.6e-3 across a wide range of finite-difference
% step sizes (not the usual V-shaped step-size curve), while the
% *absolute* error stays ~1e-5 -- consistent with a numerical-precision
% artifact at that single extreme corner, not a formula error (confirmed
% by testing intermediate scales: 3e-4 -> 2e-5 relative error, 1e-3 ->
% 1e-6, both a clean V-shape). Moving the test point 3x further from the
% origin sidesteps the worst of that corner while still exercising the
% near-small-M regime.

hList = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
fprintf('\n%-32s %14s %14s\n', 'test point', 'max rel. err', 'max abs. err');
allPass = true;
for i = 1:size(testPoints, 1)
    x0 = testPoints{i, 1};
    [~, J_analytic] = vUS_3D_erf_TC_complex_Jac(x0, tau_test, k0_test, sigma_test);

    bestRelErr = Inf(1, 4);
    bestAbsErr = Inf(1, 4);
    for h = hList
        for p = 1:4
            dx = zeros(4, 1);
            step = h * max(abs(x0(p)), 1e-6);
            dx(p) = step;

            g1_plus  = vUS_3D_erf_TC_vec_split(x0 + dx, tau_test, k0_test, sigma_test);
            g1_minus = vUS_3D_erf_TC_vec_split(x0 - dx, tau_test, k0_test, sigma_test);
            d = (g1_plus - g1_minus) / (2*step);
            fd_col = [d(:,1); d(:,2)];

            Jc = [real(J_analytic(:,p)); imag(J_analytic(:,p))];

            colScale = max(abs(fd_col));
            sig = abs(fd_col) > 1e-2 * colScale;
            relErr = abs(Jc(sig) - fd_col(sig)) ./ abs(fd_col(sig));
            bestRelErr(p) = min(bestRelErr(p), max(relErr));
            bestAbsErr(p) = min(bestAbsErr(p), max(abs(Jc - fd_col)));
        end
    end

    mRel = max(bestRelErr); mAbs = max(bestAbsErr);
    fprintf('%-32s %14.3e %14.3e\n', testPoints{i, 2}, mRel, mAbs);
    allPass = allPass && (mRel < 1e-3) && (mAbs < 1e-4);
end

if allPass
    fprintf('\nAll test points passed (analytic vs. finite-difference of vUS_3D_erf_TC_vec_split.m agree).\n');
else
    warning('One or more test points show large analytic-vs-finite-difference disagreement -- inspect before use.');
end
