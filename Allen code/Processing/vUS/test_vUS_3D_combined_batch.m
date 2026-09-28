%% Description:
%   Self-contained check of the accelerated + batched combined-parameter model (no data files needed),
%   the analog of test_vUS_3D_quad_batch.m / test_vUS_3D_quad_RCA_xy_batch.m for the much simpler
%   C-only model (vUS_3D_combined_split.m):
%     1. vUS_3D_combined_vec.m (plain complex) vs. vUS_3D_combined_split.m (the existing, original
%        [real, imag] implementation) -- must agree to round-off (same formula)
%     2. vUS_3D_combined_batchModel.m's analytic Jacobian vs. central finite differences of
%        vUS_3D_combined_vec.m -- an INDEPENDENT check (vec.m contains no Jacobian code at all)
%     3a. fitBatchedLM.m with the analytic Jacobian vs. fitBatchedLM.m with its own forward-difference
%         Jacobian on the same model: isolates the analytic Jacobian's effect from the optimizer
%     3b. fitBatchedLM.m (analytic) vs. a per-voxel lsqnonlin loop using vUS_3D_combined_split.m directly
%         (numerical Jacobian, a DIFFERENT algorithm): final cost, wall time, and recovery of v_zgp and C
%   Prints PASS/FAIL per check. Requires the Optimization Toolbox (for lsqnonlin).

rng(1);
k0 = 2*pi/(1540/13.6e6);                       % 13.6 MHz, c = 1540 m/s
frameRate = 2500; nTau = 60; tau = (0:nTau-1)'/frameRate; t1i = 2;
xscale = [500, 1e-2, 1, 1];                    % C in units of 500 (its typical magnitude), v_zgp x10mm/s, F, DC

%% 1. Plain complex model vs. the existing [real, imag]-split implementation
fprintf('1. vUS_3D_combined_vec.m vs. vUS_3D_combined_split.m:\n');
maxd1 = 0;
for trial = 1:8
    xr = [200 + 2000*rand, 40e-3*(2*rand-1), 0.4+0.6*rand, 0.3*rand];
    gv = vUS_3D_combined_vec(xr, tau, k0);
    gs = vUS_3D_combined_split(xr, tau, k0);
    maxd1 = max(maxd1, max(abs(gv - complex(gs(:,1), gs(:,2)))));
end
ok1 = maxd1 < 1e-12;
fprintf('   max|dg1| over 8 random points = %.2e -> %s\n', maxd1, passfail(ok1));

%% 2. Batched model's analytic Jacobian vs. central finite differences of vUS_3D_combined_vec.m
fprintf('\n2. vUS_3D_combined_batchModel.m Jacobian vs. finite differences of vUS_3D_combined_vec.m:\n');
tau_fd = (1:40)'/2500;
testPoints = { ...
    [ 500;  8e-3; 0.8; 0.10], 'typical flow'; ...
    [  50; -20e-3; 0.9; 0.00], 'slow decay, fast axial'; ...
    [3000;  1e-4; 0.5; 0.20], 'fast decay, near-stagnant axial'; ...
    [ 800; -5e-3; 0.6; 0.05], 'moderate'};
hList = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5];
ok2 = true;
fprintf('   %-34s %14s\n', 'test point', 'max rel. err');
for i = 1:size(testPoints, 1)
    x0 = testPoints{i, 1};
    [~, Ja] = vUS_3D_combined_batchModel(x0, tau_fd, k0);              % Bc = 1
    Ja = reshape(Ja, numel(tau_fd), 4);
    bestRel = Inf(1, 4);
    for h = hList
        for p = 1:4
            dx = zeros(4, 1); dx(p) = h * max(abs(x0(p)), 1e-6);
            gp = vUS_3D_combined_vec(x0 + dx, tau_fd, k0);
            gm = vUS_3D_combined_vec(x0 - dx, tau_fd, k0);
            fd = (gp - gm) / (2*dx(p));
            fd_col = [real(fd); imag(fd)];  Jc = [real(Ja(:,p)); imag(Ja(:,p))];
            colScale = max(abs(fd_col));
            sig = abs(fd_col) > 1e-2 * colScale;
            relErr = max(abs(Jc(sig) - fd_col(sig)) ./ abs(fd_col(sig)));
            bestRel(p) = min(bestRel(p), relErr);
        end
    end
    fprintf('   %-34s %14.3e\n', testPoints{i, 2}, max(bestRel));
    ok2 = ok2 && max(bestRel) < 1e-4;
end
fprintf('   -> %s\n', passfail(ok2));

%% Synthetic voxels for the fit comparison: model + complex noise, random fit windows
nV = 200; win = randi([8 40], nV, 1); tdi = t1i + win - 1;
truth = [50 + 2950*rand(nV,1), (-25 + 50*rand(nV,1))*1e-3, 0.5 + 0.35*rand(nV,1), 0.15*rand(nV,1)];
g1 = zeros(nV, nTau);
for i = 1:nV
    g1(i,:) = (vUS_3D_combined_vec(truth(i,:), tau, k0) + 0.03*(randn(nTau,1) + 1i*randn(nTau,1))/sqrt(2)).';
end
Vz0 = truth(:,2) + 3e-3*randn(nV,1);

%% 3a. fitBatchedLM: analytic Jacobian vs. its own forward-difference Jacobian (same solver, same model)
fprintf('\n3a. fitBatchedLM: analytic Jacobian (vUS_3D_combined_batchModel.m) vs. its own forward-difference Jacobian:\n');
lb = [0, -50e-3, 0, 0];  ub = [5000, 50e-3, 1, 1];
x0 = [500*ones(nV,1), Vz0, ones(nV,1), zeros(nV,1)];
model = @(X, cols) vUS_3D_combined_batchModel(X, tau(cols), k0);
optsFB = struct('window', tdi, 't1i', t1i, 'jacobian', 'analytic', 'xscale', xscale);
tic
[x, cost, ~, conv] = fitBatchedLM(model, x0, g1, lb, ub, optsFB);
tBatch = toc;

optsFD = struct('window', tdi, 't1i', t1i, 'jacobian', 'forward', 'xscale', xscale);
[~, costFD, ~, convFD] = fitBatchedLM(model, x0, g1, lb, ub, optsFD);
ratio3a = cost./costFD;
ok3a = mean(conv) >= 0.99 && mean(convFD) >= 0.99 && abs(median(ratio3a) - 1) < 1e-3;
fprintf('   converged: analytic %.1f%%, forward-diff %.1f%% | cost ratio: median %.6f, analytic worse by >1%% in %d, better in %d (of %d) -> %s\n', ...
    100*mean(conv), 100*mean(convFD), median(ratio3a), sum(ratio3a > 1.01), sum(ratio3a < 0.99), nV, passfail(ok3a));

%% 3b. fitBatchedLM (analytic) vs. per-voxel lsqnonlin using vUS_3D_combined_split.m (numerical Jacobian)
fprintf('\n3b. fitBatchedLM (analytic) vs. per-voxel lsqnonlin with vUS_3D_combined_split.m (numerical Jacobian):\n');
o = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', false, 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
cL = zeros(nV,1); xL = zeros(nV,4);
tic
for q = 1:nV
    inds = t1i:tdi(q); gv = g1(q,:).'; gv_split = [real(gv(inds)), imag(gv(inds))];
    x0q = [500, Vz0(q), 1, 0]; lbq = [lb(1), Vz0(q)-0.05, lb(3:4)]; ubq = [ub(1), Vz0(q)+0.05, ub(3:4)];
    resfun = @(xs) reshape(vUS_3D_combined_split(xs(:).*xscale(:), tau(inds), k0) - gv_split, [], 1);
    [xs, cL(q)] = lsqnonlin(resfun, x0q./xscale, lbq./xscale, ubq./xscale, o);
    xL(q,:) = xs.*xscale;
end
tSerial = toc;

ratio = cost./cL;
fprintf('   converged %.1f%% | cost/lsqnonlin cost: median %.6f, batched worse by >1%% in %d of %d voxels, better by >1%% in %d\n', ...
    100*mean(conv), median(ratio), sum(ratio > 1.01), nV, sum(ratio < 0.99));
fprintf('   median |v_zgp - truth|: batched %.2f mm/s, lsqnonlin %.2f mm/s | median |C - truth|: batched %.2f, lsqnonlin %.2f\n', ...
    1e3*median(abs(x(:,2) - truth(:,2))), 1e3*median(abs(xL(:,2) - truth(:,2))), ...
    median(abs(x(:,1) - truth(:,1))), median(abs(xL(:,1) - truth(:,1))));
fprintf('   wall time: batched %.3f s, serial lsqnonlin (numerical Jacobian) %.2f s (%.1fx)\n', tBatch, tSerial, tSerial/tBatch);
ok3b = mean(conv) >= 0.99 && abs(median(ratio) - 1) < 0.05;
fprintf('   -> %s (needs: >=99%% converged, median cost ratio within 5%% of 1)\n', passfail(ok3b));

fprintf('\n');
if ok1 && ok2 && ok3a && ok3b, fprintf('All checks passed.\n'); else, warning('One or more checks failed -- inspect before use.'); end

function s = passfail(tf)
    if tf, s = 'PASS'; else, s = 'FAIL'; end
end
