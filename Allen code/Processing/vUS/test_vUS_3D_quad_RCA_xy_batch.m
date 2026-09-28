%% Description:
%   Self-contained check of the accelerated + batched two-step RCA-PSF g1 model (no data files needed),
%   the RCA-PSF analog of test_vUS_3D_quad_batch.m:
%     1. vUS_3D_quad_RCA_xy_vec.m (fast, fixed quadrature) vs. vUS_3D_num_RCA_vec.m (the original,
%        independently-written adaptive-integral model this whole file accelerates)
%     2. vUS_3D_quad_RCA_xy_batchModel.m's analytic Jacobian vs. central finite differences of
%        vUS_3D_quad_RCA_xy_vec.m -- an INDEPENDENT check (vec.m contains no Jacobian code at all)
%     3. Isotropic reduction (sigma_narrow == sigma_wide): vUS_3D_quad_RCA_xy_batchModel.m must match
%        the already-validated single-Gaussian vUS_3D_quad_batchModel.m exactly
%     4. (best-effort) batched model/Jacobian vs. the symbolically-derived vUS_3D_quad_RCA_xy_complex_Jac.m
%        -- SKIPPED with a note if vUS_3D_quad_RCA_xy_complex_Jac_raw.m has not been generated yet (needs
%        the Symbolic Math Toolbox: run generate_vUS_3D_quad_RCA_xy_Jac.m once)
%     5a. fitBatchedLM.m with the analytic Jacobian (vUS_3D_quad_RCA_xy_batchModel.m) vs. fitBatchedLM.m
%         with its OWN forward-difference Jacobian on the same model: isolates the analytic Jacobian's
%         effect from any difference in optimizer algorithm
%     5b. fitBatchedLM.m (analytic) vs. a per-voxel lsqnonlin loop (numerical Jacobian, a DIFFERENT
%         algorithm): final cost, wall time, and recovery of v_zgp and the transverse speed
%         v_t = sqrt(v_xgp^2 + v_ygp^2), reported for context (looser pass bar than 5a -- see its header)
%   Prints PASS/FAIL/SKIP per check. Requires the Optimization Toolbox (for lsqnonlin); check 4 also
%   needs the Symbolic Math Toolbox (only to have been run once, not at test time).

rng(1);
sigma = [57.6, 286.5, 50.1]*1e-6;              % field sigmas, RC15gV 5 x 2 angles -6 to 6 deg (fit_RCA_PSF_v2.m)
k0 = 2*pi/(1540/13.6e6);                       % 13.6 MHz, c = 1540 m/s
frameRate = 2500; nTau = 60; tau = (0:nTau-1)'/frameRate; t1i = 2;
[s, w] = gaussLegendre01(48); xscale = [1e-2 1e-2 1e-2 1 1 1 1];

%% 1. Fast quadrature model vs. the original adaptive-integral model (vUS_3D_num_RCA_vec.m)
fprintf('1. vUS_3D_quad_RCA_xy_vec.m vs. vUS_3D_num_RCA_vec.m:\n');
maxd1 = 0;
for trial = 1:8
    xr = [40e-3*(2*rand-1), 40e-3*(2*rand-1), 40e-3*(2*rand-1), 0.4+0.6*rand, 0.3*rand, 2+rand, 0.05+0.9*rand];
    gnum = vUS_3D_num_RCA_vec(xr, tau(2:end), k0, sigma);              % vUS_3D_num_RCA_vec.m is undefined/slow at tau = 0
    gq   = vUS_3D_quad_RCA_xy_vec(xr.', tau(2:end), k0, sigma, s, w);
    maxd1 = max(maxd1, max(abs(gnum(:) - gq(:))));
end
% Threshold is integral()'s own default RelTol (1e-6), not vUS_3D_quad_RCA_xy_vec.m's accuracy: doubling
% the quadrature nodes (48 -> 192) changes its output by ~1e-10, so the ~1e-8 gap seen here is
% vUS_3D_num_RCA_vec.m's adaptive-integral tolerance floor, not disagreement between the two models.
ok1 = maxd1 < 1e-6;
fprintf('   max|dg1| over 8 random points = %.2e -> %s\n', maxd1, passfail(ok1));

%% 2. Batched model's analytic Jacobian vs. central finite differences of vUS_3D_quad_RCA_xy_vec.m
% (vec.m shares no code with the Jacobian, so this is a genuinely independent check.)
fprintf('\n2. vUS_3D_quad_RCA_xy_batchModel.m Jacobian vs. finite differences of vUS_3D_quad_RCA_xy_vec.m:\n');
tau_fd = (1:40)'/2500;
testPoints = { ...
    [ 8e-3; -5e-3; -8e-3;  0.8; 0.10; 2.5; 0.70], 'typical flow'; ...
    [ 5e-3;  8e-3;  5e-3;  0.9; 0.00; 2.0; 1.00], 'k = 2, a = 1'; ...
    [12e-3; -3e-3; -3e-3;  0.6; 0.20; 3.0; 0.40], 'k = 3, moderate a'; ...
    [ 1e-4;  1e-4;  1e-4;  0.5; 0.05; 2.5; 0.90], 'near-stagnant (small M)'; ...
    [ 8e-3; -5e-3; -8e-3;  0.8; 0.10; 2.5; 0.05], 'small a'; ...
    [-8e-3;  5e-3;  15e-3; 0.9; 0.00; 3.0; 0.90], 'fast axial, negative vxgp'};
hList = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
ok2 = true;
fprintf('   %-30s %14s\n', 'test point', 'max rel. err');
for i = 1:size(testPoints, 1)
    x0 = testPoints{i, 1};
    [~, Ja] = vUS_3D_quad_RCA_xy_batchModel(x0, tau_fd, k0, sigma, s, w);   % Bc = 1
    Ja = reshape(Ja, numel(tau_fd), 7);
    bestRel = Inf(1, 7);
    for h = hList
        for p = 1:7
            dx = zeros(7, 1); dx(p) = h * max(abs(x0(p)), 1e-6);
            gp = vUS_3D_quad_RCA_xy_vec(x0 + dx, tau_fd, k0, sigma, s, w);
            gm = vUS_3D_quad_RCA_xy_vec(x0 - dx, tau_fd, k0, sigma, s, w);
            fd = (gp - gm) / (2*dx(p));
            fd_col = [real(fd); imag(fd)];  Jc = [real(Ja(:,p)); imag(Ja(:,p))];
            colScale = max(abs(fd_col));
            sig = abs(fd_col) > 1e-2 * colScale;
            relErr = max(abs(Jc(sig) - fd_col(sig)) ./ abs(fd_col(sig)));
            bestRel(p) = min(bestRel(p), relErr);
        end
    end
    fprintf('   %-30s %14.3e\n', testPoints{i, 2}, max(bestRel));
    ok2 = ok2 && max(bestRel) < 1e-4;
end
fprintf('   -> %s\n', passfail(ok2));

%% 3. Isotropic reduction: vUS_3D_quad_RCA_xy_batchModel.m vs. vUS_3D_quad_batchModel.m
fprintf('\n3. isotropic reduction (sigma_narrow == sigma_wide) vs. the existing vUS_3D_quad_batchModel.m:\n');
sigma_iso = [130.1474, 130.1474, 44.7070]*1e-6;
nChk = 20;
X6 = [ (2+18*rand(1,nChk))*1e-3; (-25+50*rand(1,nChk))*1e-3; 0.4+0.6*rand(1,nChk); 0.3*rand(1,nChk); 2+rand(1,nChk); 0.05+0.9*rand(1,nChk) ];
X7 = [X6(1,:); zeros(1,nChk); X6(2:6,:)];                                 % v_xgp = v_tgp, v_ygp = 0
[Y6, J6] = vUS_3D_quad_batchModel(X6, tau_fd, k0, sigma_iso, s, w);
[Y7, J7] = vUS_3D_quad_RCA_xy_batchModel(X7, tau_fd, k0, sigma_iso, s, w);
dY = max(abs(Y6 - Y7), [], 'all');
dJ3 = max(abs(cat(2, J6(:,1,:), J6(:,2:6,:)) - cat(2, J7(:,1,:), J7(:,3:7,:))), [], 'all');   % v_xgp<->v_tgp column; v_ygp has no counterpart
ok3 = dY < 1e-12 && dJ3 < 1e-9;
fprintf('   max|dY| = %.2e, max|dJ| = %.2e -> %s\n', dY, dJ3, passfail(ok3));

%% 4. (best-effort) batched model vs. the symbolically-derived single-voxel Jacobian
fprintf('\n4. batched model vs. vUS_3D_quad_RCA_xy_complex_Jac.m (needs the generated _raw.m file):\n');
if exist('vUS_3D_quad_RCA_xy_complex_Jac_raw', 'file') == 2
    maxdr4 = 0; maxdJ4 = 0;
    for i = 1:size(testPoints, 1)
        x0 = testPoints{i, 1};
        [Y4, J4] = vUS_3D_quad_RCA_xy_batchModel(x0, tau_fd, k0, sigma, s, w);
        [g4, Jc4] = vUS_3D_quad_RCA_xy_complex_Jac(x0.', tau_fd, k0, sigma, s, w);
        maxdr4 = max(maxdr4, max(abs(Y4(:) - g4(:))));
        maxdJ4 = max(maxdJ4, max(abs(reshape(J4, numel(tau_fd), 7) - Jc4), [], 'all'));
    end
    ok4 = maxdr4 < 1e-9 && maxdJ4 < 1e-9;
    fprintf('   max|dg1| = %.2e, max|dJ| = %.2e -> %s\n', maxdr4, maxdJ4, passfail(ok4));
else
    ok4 = true;   % not a failure -- just not generated in this environment
    fprintf('   SKIPPED: vUS_3D_quad_RCA_xy_complex_Jac_raw.m not found. Run generate_vUS_3D_quad_RCA_xy_Jac.m once\n');
    fprintf('   (needs the Symbolic Math Toolbox) to generate it and enable this cross-check.\n');
end

%% Synthetic voxels for the fit comparison: model + complex noise, random fit windows, random flow direction
nV = 200; win = randi([8 40], nV, 1); tdi = t1i + win - 1;
vt_truth = (2 + 18*rand(nV,1))*1e-3; phi_truth = (pi/2)*rand(nV,1);           % v_t in [2,20] mm/s, phi in [0, 90 deg]
truth = [vt_truth.*cos(phi_truth), vt_truth.*sin(phi_truth), (-25 + 50*rand(nV,1))*1e-3, ...
         0.5 + 0.35*rand(nV,1), 0.15*rand(nV,1), 2 + rand(nV,1), 0.3 + 0.7*rand(nV,1)];
g1 = zeros(nV, nTau);
for i = 1:nV
    g1(i,:) = (vUS_3D_quad_RCA_xy_vec(truth(i,:).', tau, k0, sigma, s, w) + 0.03*(randn(nTau,1) + 1i*randn(nTau,1))/sqrt(2)).';
end
Vz0 = truth(:,3) + 3e-3*randn(nV,1);

%% 5a. fitBatchedLM with the analytic Jacobian (vUS_3D_quad_RCA_xy_batchModel.m) vs. fitBatchedLM with
% its OWN forward-finite-difference Jacobian: same solver algorithm, same model, so any gap is purely
% about the Jacobian, not a difference in optimizer. This is the cleanest correctness/quality check for
% the analytic Jacobian, and needs no generated file.
fprintf('\n5a. fitBatchedLM: analytic Jacobian (vUS_3D_quad_RCA_xy_batchModel.m) vs. its own forward-difference Jacobian:\n');
lb = [-50e-3, -50e-3, -50e-3, 0, 0, 2, 0];  ub = [50e-3, 50e-3, 50e-3, 1, 1, 3, 1];
x0 = [5e-3*ones(nV,1), 5e-3*ones(nV,1), Vz0, ones(nV,1), zeros(nV,1), 2.5*ones(nV,1), ones(nV,1)];
model = @(X, cols) vUS_3D_quad_RCA_xy_batchModel(X, tau(cols), k0, sigma, s, w);
optsFB = struct('window', tdi, 't1i', t1i, 'jacobian', 'analytic', 'xscale', xscale);
tic
[x, cost, iters, conv] = fitBatchedLM(model, x0, g1, lb, ub, optsFB);
tBatch = toc;

modelNoJ = @(X, cols) vUS_3D_quad_RCA_xy_batchModel(X, tau(cols), k0, sigma, s, w);   % same model, called with nargout = 1 below
optsFD = struct('window', tdi, 't1i', t1i, 'jacobian', 'forward', 'xscale', xscale);
[xFD, costFD, ~, convFD] = fitBatchedLM(modelNoJ, x0, g1, lb, ub, optsFD);
ratio5a = cost./costFD;
ok5a = mean(conv) >= 0.99 && mean(convFD) >= 0.99 && abs(median(ratio5a) - 1) < 1e-3;
fprintf('   converged: analytic %.1f%%, forward-diff %.1f%% | cost ratio: median %.6f, analytic worse by >1%% in %d, better in %d (of %d) -> %s\n', ...
    100*mean(conv), 100*mean(convFD), median(ratio5a), sum(ratio5a > 1.01), sum(ratio5a < 0.99), nV, passfail(ok5a));

%% 5b. fitBatchedLM (analytic Jacobian) vs. a per-voxel lsqnonlin loop with a NUMERICAL Jacobian (lsqnonlin
% has no way to use vUS_3D_quad_RCA_xy_batchModel.m's Jacobian directly, and the symbolically-generated
% single-voxel one needs check 4's file). This compares two different ALGORITHMS, not just the Jacobian,
% so some fraction of voxels going the other way is expected (confirmed by 5a: the batched solver itself
% performs identically with or without an analytic Jacobian) -- reported for wall-time context, with a
% looser pass bar than 5a.
fprintf('\n5b. fitBatchedLM (analytic) vs. per-voxel lsqnonlin (numerical Jacobian, different algorithm):\n');
o = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', false, 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
cL = zeros(nV,1); xL = zeros(nV,7);
tic
for q = 1:nV
    inds = t1i:tdi(q); gv = g1(q,:).'; gv_split = [real(gv(inds)); imag(gv(inds))];
    x0q = [5e-3, 5e-3, Vz0(q), 1, 0, 2.5, 1]; lbq = [lb(1:2), Vz0(q)-0.05, lb(4:end)]; ubq = [ub(1:2), Vz0(q)+0.05, ub(4:end)];
    resfun = @(xs) reshape([real(vUS_3D_quad_RCA_xy_vec(xs(:).*xscale(:), tau(inds), k0, sigma, s, w)); ...
                             imag(vUS_3D_quad_RCA_xy_vec(xs(:).*xscale(:), tau(inds), k0, sigma, s, w))], [], 1) - gv_split;
    [xs, cL(q)] = lsqnonlin(resfun, x0q./xscale, lbq./xscale, ubq./xscale, o);
    xL(q,:) = xs.*xscale;
end
tSerial = toc;

ratio = cost./cL;
vt_batch = sqrt(x(:,1).^2 + x(:,2).^2); vt_lsq = sqrt(xL(:,1).^2 + xL(:,2).^2);
fprintf('   converged %.1f%% | cost/lsqnonlin cost: median %.6f, batched worse by >1%% in %d of %d voxels, better by >1%% in %d\n', ...
    100*mean(conv), median(ratio), sum(ratio > 1.01), nV, sum(ratio < 0.99));
fprintf('   median |v_zgp - truth|: batched %.2f mm/s, lsqnonlin %.2f mm/s | median |v_t - truth|: batched %.2f mm/s, lsqnonlin %.2f mm/s\n', ...
    1e3*median(abs(x(:,3) - truth(:,3))), 1e3*median(abs(xL(:,3) - truth(:,3))), ...
    1e3*median(abs(vt_batch - vt_truth)), 1e3*median(abs(vt_lsq - vt_truth)));
fprintf('   wall time: batched %.2f s, serial lsqnonlin (numerical Jacobian) %.2f s (%.1fx)\n', tBatch, tSerial, tSerial/tBatch);
ok5b = mean(conv) >= 0.99 && abs(median(ratio) - 1) < 0.05;
fprintf('   -> %s (needs: >=99%% converged, median cost ratio within 5%% of 1; a fair number of voxels going\n', passfail(ok5b));
fprintf('      either way is expected here -- see 5a for the apples-to-apples Jacobian check)\n');

fprintf('\n');
if ok1 && ok2 && ok3 && ok4 && ok5a && ok5b, fprintf('All checks passed (or skipped where the Symbolic Math Toolbox artifact was unavailable).\n'); else, warning('One or more checks failed -- inspect before use.'); end

function s = passfail(tf)
    if tf, s = 'PASS'; else, s = 'FAIL'; end
end
