%% Description:
%   Symbolically derives the analytic Jacobian of the new (Poiseuille /
%   uniform-velocity-distribution) g1 model -- vUS_2D_erf.m /
%   vUS_2D_erf_vec_split.m -- w.r.t. its four free parameters
%   [v_xgp, v_zgp, F, DC] (DC a real additive offset), and auto-generates
%   a fast numeric MATLAB function from the result via matlabFunction.
%
%   Run this script once (and again any time the model formula itself
%   changes) to (re)generate:
%       vUS_2D_erf_complex_Jac_raw.m  -- AUTO-GENERATED internal helper,
%           do not call directly or hand-edit
%       vUS_2D_erf_complex_Jac.m      -- AUTO-GENERATED public [g1,J]
%           interface, do not hand-edit
%   Then use the small hand-written wrappers vUS_2D_erf_vec_split_Jac.m
%   and vUS_2D_erf_residJac.m to plug into lsqnonlin (see their headers
%   for usage).
%
%   Rationale: hand-differentiating this model has several chain-rule
%   layers -- M(v_xgp, v_zgp), the exp(-4*k0^2*v_zgp^2/M) prefactor, and
%   two erf arguments that both depend on v_zgp AND M -- and is easy to
%   get subtly wrong (see g1vUS2D_Jac.m's header for how much bookkeeping
%   the OLD model's hand-derived version needed, and that model has one
%   fewer parameter). Symbolic differentiation removes that risk.
%
%   Two-file design (raw helper + thin public wrapper), instead of a
%   single matlabFunction(g1, J, ...) call: dg1/dDC is an exact constant
%   (1), which matlabFunction correctly simplifies away -- but that then
%   makes its internal J=[col1,col2,col3,col4] concatenation crash,
%   since a bare scalar doesn't auto-broadcast against the other
%   tau-length columns under horzcat (unlike elementwise operators).
%   Forcing a fake "+0*tau" dependency to dodge this gets optimized back
%   out by matlabFunction's own CSE pass (tested, confirmed on R2025b /
%   Symbolic Math Toolbox 25.2) -- so instead each column is generated as
%   its own separate output (no internal concatenation at all), and the
%   public wrapper below assembles + explicitly broadcasts them.
%
%   This script also self-checks the generated Jacobian against finite
%   differences at a few test points, including a near-stagnant-flow
%   case (small M) -- the catastrophic-cancellation regime flagged in
%   review, where analytic differentiation matters most -- before you
%   rely on it.
%
% Requires: Symbolic Math Toolbox

%% 1. Build the model symbolically
x     = sym('x', [4, 1], 'real');       % x = [v_xgp; v_zgp; F; DC]
sigma = sym('sigma', [2, 1], 'real');   % sigma = [sigma_x; sigma_z]
syms tau k0 real

v_xgp = x(1); v_zgp = x(2); Fp = x(3); DC = x(4);
sx = sigma(1); sz = sigma(2);

M = v_xgp^2/sx^2 + v_zgp^2/sz^2;

A = sqrt(M)*tau - 2*1i*k0*v_zgp/sqrt(M);
B = -2*1i*k0*v_zgp/sqrt(M);

g1 = DC + Fp/2 * sqrt(sym(pi)/M)/tau * exp(-4*k0^2*v_zgp^2/M) * (erf(A) - erf(B));

%% 2. Symbolic Jacobian w.r.t. [v_xgp, v_zgp, F, DC]
% Differentiation w.r.t. a REAL parameter commutes with real()/imag(), so
% this single complex-valued Jacobian is all that's needed -- the
% real/imag split for lsqnonlin happens downstream in
% vUS_2D_erf_vec_split_Jac.m without any further differentiation.
J = jacobian(g1, x);   % 1x4, complex-valued

%% 3. Generate the raw numeric code -- g1 + 4 SEPARATE Jacobian-column
% outputs (see header for why this avoids matlabFunction's internal
% concatenation crash). Note: 'Vectorize' is intentionally omitted --
% it's not a recognized matlabFunction option as of Symbolic Math Toolbox
% 25.2 (R2025b); elementwise (.^, .*, ./) code generation is the default
% now, confirmed by inspecting the generated output.
vUSDir = fileparts(mfilename('fullpath'));
rawFile = fullfile(vUSDir, 'vUS_2D_erf_complex_Jac_raw.m');
matlabFunction(g1, J(1), J(2), J(3), J(4), ...
    'File', rawFile, ...
    'Vars', {x, tau, k0, sigma}, ...
    'Outputs', {'g1', 'dg1_dvxgp', 'dg1_dvzgp', 'dg1_dF', 'dg1_dDC'}, ...
    'Optimize', true);

% Swap the symbolic 'erf' calls for the codebase's complex-safe 'erfz'.
% matlabFunction emits erf(...); MATLAB's built-in erf's complex-argument
% support varies by release, and erfz (Allen Code/ErrorFunction) is what
% the rest of this codebase already uses/trusts for complex arguments --
% keep everything on one implementation for numerical consistency, and so
% this generated file's g1 values match vUS_2D_erf_vec_split.m's exactly.
txt = fileread(rawFile);
txt = regexprep(txt, '(?<![A-Za-z0-9_])erf\(', 'erfz(');
banner = sprintf(['%%%% AUTO-GENERATED (internal helper) by generate_vUS_2D_erf_Jac.m on %s.\n', ...
    '%%%% Do not call directly or hand-edit -- use vUS_2D_erf_complex_Jac.m,\n', ...
    '%%%% and re-run the generator script to update both.\n\n'], char(datetime('now')));
fid = fopen(rawFile, 'w');
fwrite(fid, [banner, txt]);
fclose(fid);
fprintf('Generated %s\n', rawFile);

%% 4. Write the thin public wrapper: standard [g1, J] interface (matching
% vUS_2D_erf_vec_split_Jac.m / vUS_2D_erf_residJac.m's expectations),
% with each Jacobian column explicitly broadcast to size(tau) before
% concatenation.
outFile = fullfile(vUSDir, 'vUS_2D_erf_complex_Jac.m');
pubCode = sprintf([ ...
    '%%%% AUTO-GENERATED by generate_vUS_2D_erf_Jac.m on %s.\n', ...
    '%%%% Do not hand-edit -- re-run the generator script instead.\n', ...
    '%%%% Public [g1, J] wrapper around vUS_2D_erf_complex_Jac_raw.m.\n\n', ...
    'function [g1, J] = vUS_2D_erf_complex_Jac(x, tau, k0, sigma)\n', ...
    '    tau = tau(:);\n', ...
    '    if nargout > 1\n', ...
    '        [g1, dvx, dvz, dF, dDC] = vUS_2D_erf_complex_Jac_raw(x, tau, k0, sigma);\n', ...
    '        onesTau = ones(size(tau));\n', ...
    '        J = [dvx.*onesTau, dvz.*onesTau, dF.*onesTau, dDC.*onesTau];\n', ...
    '    else\n', ...
    '        g1 = vUS_2D_erf_complex_Jac_raw(x, tau, k0, sigma);\n', ...
    '    end\n', ...
    '    g1 = g1(:) .* ones(size(tau));\n', ...
    'end\n'], char(datetime('now')));
fid = fopen(outFile, 'w');
fwrite(fid, pubCode);
fclose(fid);
rehash path
fprintf('Generated %s\n', outFile);

%% 5. Self-check: compare the analytic Jacobian against finite differences
% Test points include a near-stagnant-flow case (small M, the
% catastrophic-cancellation regime) plus typical and fast/negative-vxgp
% cases, using this dataset's actual PSF/transducer parameters.
sigma_test = [58.9110, 73.9967] * 1e-6;    % [m], from vUS_2D_newmodel.m
k0_test    = 2*pi / (1540/15.625e6);       % [rad/m], L22-14v @ 15.625 MHz
tau_test   = (1:40)' / 5000;               % [s], 5 kHz frame rate, tau1..tau40

testPoints = { ...
    [ 3e-3;   8e-3; 0.7; 0.05], 'typical flow'; ...
    [ 1e-4;   1e-4; 0.5; 0.05], 'near-stagnant (small M)'; ...
    [-5e-3;  15e-3; 0.9; 0.00], 'fast axial, negative vxgp'};

% Intermediate terms (e.g. exp(-B^2)) can be huge (~1e15) even when g1
% itself is small -- k0*sigma is an O(1-10) geometric constant of the
% transducer/PSF, independent of flow speed -- so a single fixed
% finite-difference step is not robust across all test points. Sweep
% several step sizes per parameter and keep the best (minimum error)
% found, the standard way to check a gradient against finite differences
% for a function with large dynamic range.
hList = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
fprintf('\n%-32s %14s %14s\n', 'test point', 'max rel. err', 'max abs. err');
allPass = true;
for i = 1:size(testPoints, 1)
    x0 = testPoints{i, 1};
    [~, J_analytic] = vUS_2D_erf_complex_Jac(x0, tau_test, k0_test, sigma_test);

    bestRelErr = Inf(1, 4);
    bestAbsErr = Inf(1, 4);
    for h = hList
        for p = 1:4
            dx = zeros(4, 1);
            step = h * max(abs(x0(p)), 1e-6);
            dx(p) = step;
            g1_plus  = vUS_2D_erf_complex_Jac(x0 + dx, tau_test, k0_test, sigma_test);
            g1_minus = vUS_2D_erf_complex_Jac(x0 - dx, tau_test, k0_test, sigma_test);
            fd_col = (g1_plus - g1_minus) / (2*step);

            % Score relative error only on entries that are a
            % non-negligible fraction of this column's own peak
            % magnitude -- near-zero entries (e.g. Im(dg1/dv_xgp) as
            % tau->0, since g1(tau=0)=1 regardless of parameters)
            % inflate relative error on noise both estimates already
            % agree is tiny. Absolute error has no denominator, so
            % nothing is masked there.
            colScale = max(abs(fd_col));
            sig = abs(fd_col) > 1e-2 * colScale;
            relErr = abs(J_analytic(sig,p) - fd_col(sig)) ./ abs(fd_col(sig));
            bestRelErr(p) = min(bestRelErr(p), max(relErr));
            bestAbsErr(p) = min(bestAbsErr(p), max(abs(J_analytic(:,p) - fd_col)));
        end
    end

    mRel = max(bestRelErr); mAbs = max(bestAbsErr);
    fprintf('%-32s %14.3e %14.3e\n', testPoints{i, 2}, mRel, mAbs);
    allPass = allPass && (mRel < 1e-3) && (mAbs < 1e-4);
end

if allPass
    fprintf('\nAll test points passed (analytic vs. finite-difference agree to <1e-3 relative).\n');
else
    warning('One or more test points show large analytic-vs-finite-difference disagreement -- inspect before use.');
end
