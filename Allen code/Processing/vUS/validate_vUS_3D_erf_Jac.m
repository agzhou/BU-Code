%% Description:
%   3D analog of validate_vUS_2D_erf_Jac.m. Validates
%   vUS_3D_erf_complex_Jac.m's analytic Jacobian against central finite
%   differences of the independently-written vUS_3D_erf_vec_split.m.
%   Uses only base MATLAB -- no Symbolic Math or Optimization Toolbox
%   required, so it can run anywhere this codebase runs, independent of
%   whether the symbolic generator is available.
%
%   Test points include a near-stagnant-flow case (small M) and use the
%   RC15gV matrix-probe PSF from test_vUS_3D_erf.m -- a genuinely
%   different geometry than the 2D L22-14v case (larger, near-isotropic
%   lateral sigma vs. a smaller axial sigma) -- plus several step sizes
%   per parameter (see generate_vUS_2D_erf_Jac.m's header for why a
%   single fixed step isn't robust: k0*sigma sets an O(1-10) scale
%   independent of flow speed, so intermediate terms can be huge even
%   when g1 itself is small).

sigma_test = [429; 429; 126] * 1e-6;       % [m], RC15gV probe (column --
                                             % vUS_3D_erf_complex_Jac.m
                                             % indexes in4(2,:), which
                                             % requires a column input)
k0_test    = 2*pi / (1540/15.625e6);       % [rad/m]
tau_test   = (1:40)' / 5000;               % [s], 5 kHz frame rate, tau1..tau40

testPoints = { ...
    [ 3e-3;  -4e-3;   8e-3; 0.7; 0.05], 'typical flow'; ...
    [ 1e-4;   1e-4;   1e-4; 0.5; 0.05], 'near-stagnant (small M)'; ...
    [-5e-3;   6e-3;  15e-3; 0.9; 0.00], 'fast axial, negative vxgp'; ...
    [ 1e-2;  -1e-2;   1e-3; 0.3; 0.10], 'fast in-plane, slow axial'};

hList = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
fprintf('%-32s %14s %14s\n', 'test point', 'max rel. err', 'max abs. err');
allPass = true;
for i = 1:size(testPoints, 1)
    x0 = testPoints{i, 1};

    [~, J_analytic] = vUS_3D_erf_vec_split_Jac(x0, tau_test, k0_test, sigma_test);

    bestRelErr = Inf(1, 5);
    bestAbsErr = Inf(1, 5);
    for h = hList
        for p = 1:5
            dx = zeros(5, 1);
            step = h * max(abs(x0(p)), 1e-6);
            dx(p) = step;

            g1_plus  = vUS_3D_erf_vec_split(x0 + dx, tau_test, k0_test, sigma_test);
            g1_minus = vUS_3D_erf_vec_split(x0 - dx, tau_test, k0_test, sigma_test);

            d = (g1_plus - g1_minus) / (2*step); % [nTau, 2]
            fd_col = [d(:,1); d(:,2)];           % stack real over imag, matching J's convention

            % Only score entries that are a non-negligible fraction of
            % this column's own peak magnitude for the relative check --
            % entries near zero (e.g. Im(dg1/dv_xgp) as tau->0, since
            % g1(tau=0)=1 regardless of parameters) inflate relative
            % error on noise both estimates already agree is tiny. The
            % absolute-error check has no denominator, so nothing is
            % masked there.
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
    fprintf('\nAll test points passed (analytic vs. finite-difference of vUS_3D_erf_vec_split.m agree).\n');
else
    warning('One or more test points show large analytic-vs-finite-difference disagreement -- inspect before use.');
end
