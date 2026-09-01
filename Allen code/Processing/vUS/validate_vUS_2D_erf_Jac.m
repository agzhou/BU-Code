%% Description:
%   Validates vUS_2D_erf_complex_Jac.m's analytic Jacobian against
%   central finite differences of the EXISTING, independently-written
%   vUS_2D_erf_vec_split.m (not against itself, so this is a genuine
%   cross-check). Uses only base MATLAB -- no Symbolic Math or
%   Optimization Toolbox required, so it can run anywhere this codebase
%   runs.
%
%   Test points include a near-stagnant-flow case (small M), the
%   catastrophic-cancellation regime flagged as a numerical-fragility
%   concern in review -- this is exactly where an analytic Jacobian
%   should behave best relative to finite differences of the already-
%   fragile function value.
%
%   Because intermediate terms (e.g. exp(-B^2)) can be huge (~1e15) even
%   when g1 itself is small -- a consequence of k0*sigma being an O(1-10)
%   geometric constant of the transducer/PSF, independent of flow speed
%   -- a single fixed finite-difference step size is not robust across
%   all test points: too small and the subtraction drowns in rounding
%   noise, too large and truncation error dominates. This script sweeps
%   several step sizes per parameter and reports the best (minimum
%   error) one found, which is the standard way to check an analytic
%   gradient against finite differences for a function with large
%   dynamic range.

sigma_test = [58.9110, 73.9967] * 1e-6;    % [m], from vUS_2D_newmodel.m
k0_test    = 2*pi / (1540/15.625e6);       % [rad/m], L22-14v @ 15.625 MHz
tau_test   = (1:40)' / 5000;               % [s], 5 kHz frame rate, tau1..tau40

testPoints = { ...
    [ 3e-3;   8e-3; 0.7; 0.05], 'typical flow'; ...
    [ 1e-4;   1e-4; 0.5; 0.05], 'near-stagnant (small M)'; ...
    [-5e-3;  15e-3; 0.9; 0.00], 'fast axial, negative vxgp'; ...
    [ 1e-2;   1e-3; 0.3; 0.10], 'fast lateral, slow axial'};

hList = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]; % relative step scales to sweep
fprintf('%-32s %14s %14s\n', 'test point', 'max rel. err', 'max abs. err');
allPass = true;
for i = 1:size(testPoints, 1)
    x0 = testPoints{i, 1};

    [~, J_analytic] = vUS_2D_erf_vec_split_Jac(x0, tau_test, k0_test, sigma_test);

    bestRelErr = Inf(1, 4); % best (min) error found per parameter, across the h sweep
    bestAbsErr = Inf(1, 4);
    for h = hList
        for p = 1:4
            dx = zeros(4, 1);
            step = h * max(abs(x0(p)), 1e-6);
            dx(p) = step;

            g1_plus  = vUS_2D_erf_vec_split(x0 + dx, tau_test, k0_test, sigma_test);
            g1_minus = vUS_2D_erf_vec_split(x0 - dx, tau_test, k0_test, sigma_test);

            d = (g1_plus - g1_minus) / (2*step); % [nTau, 2]
            fd_col = [d(:,1); d(:,2)];           % stack real over imag, matching J's convention

            % Only score entries that are a non-negligible fraction of
            % this column's own peak magnitude for the relative check --
            % entries near zero (e.g. Im(dg1/dv_xgp) as tau->0, since
            % g1(tau=0)=1 regardless of parameters) inflate relative
            % error on noise that both estimates already agree is tiny
            % and physically unimportant. The absolute-error check has
            % no denominator, so nothing is masked there.
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
    fprintf('\nAll test points passed (analytic vs. finite-difference of vUS_2D_erf_vec_split.m agree).\n');
else
    warning('One or more test points show large analytic-vs-finite-difference disagreement -- inspect before use.');
end
