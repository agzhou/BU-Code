%% Description:
%   Self-contained check of vUS_2D_quad_batchModel.m, the 2D (v_xgp, v_zgp) batched quadrature model,
%   with the input structure of vUS_2D_newmodel.m (sigma = [sigma_x, sigma_z], start point
%   [Vx0 Vz0 FR0 DCR0] + [k a], lb/ub = [0 Vz0-0.01 0 0]/[30e-3 Vz0+0.01 1 1] + bounds on k and a):
%     1. value and analytic Jacobian vs. the per-voxel vUS_3D_quad_vec.m / vUS_3D_quad_complex_Jac.m
%        (with sigma_y = sigma_x); sigma ordering and Jacobian column order/units vs. finite differences
%     2. value vs. the numerically integrated vUS_2D_num_wrapper.m (an independent implementation), g1(0) = DC + F
%     3. k = 2, a = 1 reproduces vUS_2D_erf_vec.m (the 4-parameter model vUS_2D_newmodel.m fits); needs erfz
%     4. sigma must have 2 elements (a 3D sigma is rejected)
%     5. fitBatchedLM.m fits in the vUS_2D_newmodel.m structure: all six parameters free, and k = 2, a = 1 fixed
%        (lb = ub) vs. a per-voxel lsqnonlin loop on the erf model that the script uses today
%   Prints PASS/FAIL per check. Requires the Optimization Toolbox (for lsqnonlin).

rng(4);
here = fileparts(mfilename('fullpath'));
addPaths = {fullfile(here, '..', '..', 'ErrorFunction'), fullfile(here, 'k = 2 version')};   % erfz and vUS_2D_erf_vec.m, as in vUS_2D_newmodel.m
for i = 1:numel(addPaths)
    if isfolder(addPaths{i}), addpath(addPaths{i}); end
end
try, erfz(1); hasErf = true; catch, hasErf = false; end

sigma = [60.5514 73.9511]*1e-6;                % [sigma_x, sigma_z], L22-14v (vUS_2D_newmodel.m); sigma_x ~= sigma_z so an ordering slip shows
k0 = 2*pi/(1540/15.625e6);                     % 15.625 MHz, c = 1540 m/s
frameRate = 5000; nTau = ceil(20e-3*frameRate); tau = (0:nTau-1)'/frameRate; t1i = 2;
[s, w] = gaussLegendre01(48); xscale = [1e-2 1e-2 1 1 1 1];
sigma3 = [sigma(1) sigma(1) sigma(2)];         % the same model in the 3D functions
allOk = true;

%% 1. value + Jacobian vs. the per-voxel 3D reference code, and vs. finite differences
nV = 100;
Xt = [(2 + 18*rand(nV,1))*1e-3, (-25 + 50*rand(nV,1))*1e-3, 0.5 + 0.4*rand(nV,1), 0.15*rand(nV,1), 2 + rand(nV,1), 0.3 + 0.7*rand(nV,1)].';
[Y, J] = vUS_2D_quad_batchModel(Xt, tau, k0, sigma, s, w);  Y1 = vUS_2D_quad_batchModel(Xt, tau, k0, sigma, s, w);
dY = 0; dJ = 0; Jmax = 0;
for q = 1:nV
    [gq, Jq] = vUS_3D_quad_complex_Jac(Xt(:,q).', tau, k0, sigma3, s, w);
    dY = max(dY, max(abs(Y(:,q) - gq)));  dJ = max(dJ, max(abs(J(:,:,q) - Jq), [], 'all'));  Jmax = max(Jmax, max(abs(Jq), [], 'all'));
end
ok = dY < 1e-12 && dJ/Jmax < 1e-12 && max(abs(Y1 - Y), [], 'all') < 1e-12;
allOk = report('1a. batch model/Jacobian vs vUS_3D_quad_vec / complex_Jac with sigma = [sx sx sz]: max|dY| = %.1e, max|dJ|/max|J| = %.1e (single-output call agrees: %d)', ok, dY, dJ/Jmax, max(abs(Y1 - Y), [], 'all') < 1e-12) && allOk;

% central differences of the model itself, per parameter and in physical units: checks the column order (col 1 = d/dv_xgp) and units
dq = 20; relErr = zeros(1, 6);
for j = 1:6
    h = 1e-5*xscale(j);  Xp = Xt(:,1:dq);  Xm = Xt(:,1:dq);  Xp(j,:) = Xp(j,:) + h;  Xm(j,:) = Xm(j,:) - h;
    Jfd = (vUS_2D_quad_batchModel(Xp, tau, k0, sigma, s, w) - vUS_2D_quad_batchModel(Xm, tau, k0, sigma, s, w))/(2*h);
    Jj = reshape(J(:,j,1:dq), numel(tau), dq);
    relErr(j) = max(abs(Jfd - Jj), [], 'all')/max(abs(Jj), [], 'all');
end
ok = all(relErr < 1e-6);
allOk = report('1b. Jacobian column j vs central differences of the model, [v_xgp v_zgp F DC k a]: max relative error %s', ok, mat2str(relErr, 2)) && allOk;
% swapping sigma_x and sigma_z must change the model (guards the ordering check above)
Yswap = vUS_2D_quad_batchModel(Xt, tau, k0, sigma([2 1]), s, w);
allOk = report('1c. swapping sigma_x and sigma_z changes the model (max |dY| = %.1e), so the ordering in 1a is tested', max(abs(Yswap - Y), [], 'all') > 1e-3, max(abs(Yswap - Y), [], 'all')) && allOk;

%% 2. independent reference: numerically integrated vUS_2D_num_wrapper (integral(), the 2D sigma convention)
nN = 12; dN = 0; g0 = 0;
for q = 1:nN
    x = Xt(:,q).';
    gn = vUS_2D_num_wrapper(x, tau, k0, sigma);
    dN = max(dN, max(abs(Y(:,q) - gn)));  g0 = max(g0, abs(Y(1,q) - (x(3) + x(4))));
end
ok = dN < 1e-6 && g0 < 1e-12;
allOk = report('2. vs vUS_2D_num_wrapper (integral): max|dY| = %.1e (integral() tolerance), |g1(0) - (DC + F)| = %.1e', ok, dN, g0) && allOk;

%% 3. k = 2, a = 1 is the erf model of vUS_2D_erf_vec.m
if hasErf
    dE = 0;
    for q = 1:20
        x4 = Xt(1:4, q).';
        ge = vUS_2D_erf_vec(x4, tau(2:end), k0, sigma);               % tau = 0 is 0/0 in the closed form
        gq = vUS_2D_quad_batchModel([x4 2 1].', tau, k0, sigma, s, w);
        dE = max(dE, max(abs(ge(:) - gq(2:end))));
    end
    allOk = report('3. quadrature model at k = 2, a = 1 vs vUS_2D_erf_vec: max|dg1| = %.1e (48-node quadrature)', dE < 1e-6, dE) && allOk;
else
    fprintf('3. SKIPPED: erfz not found (add Allen code/ErrorFunction to the path)\n');
end

%% 4. sigma input structure
ok = throws(@() vUS_2D_quad_batchModel(Xt, tau, k0, sigma3, s, w), 'sigmaSize') && throws(@() vUS_2D_quad_batchModel(Xt, tau, k0, sigma(1), s, w), 'sigmaSize');
Yc = vUS_2D_quad_batchModel(Xt, tau, k0, sigma(:), s, w);            % a column sigma is fine
okCol = max(abs(Yc - Y1), [], 'all') < 1e-14;
allOk = report('4. a 3-element (3D) or 1-element sigma is rejected with a clear error: %d; a column [sigma_x; sigma_z] is accepted: %d', ok && okCol, ok, okCol) && allOk;

%% 5. fitBatchedLM in the structure of vUS_2D_newmodel.m: a stacked plane [nz*nx, nTau], per-voxel Vz0 and fit window
nz = 20; nx = 15; nP = nz*nx; win = randi([8 40], nP, 1); tdi = t1i + win - 1;
truth = [(2 + 18*rand(nP,1))*1e-3, (-25 + 50*rand(nP,1))*1e-3, 0.5 + 0.35*rand(nP,1), 0.15*rand(nP,1), 2 + rand(nP,1), 0.3 + 0.7*rand(nP,1)];
noise = @() 0.03*(randn(nP, nTau) + 1i*randn(nP, nTau))/sqrt(2);
Vz0 = truth(:,2) + 3e-3*randn(nP,1);  Vx0 = 5e-3*ones(nP,1);         % as in the script: findVzPhaseDiff output and a uniform v_xgp guess
model = @(X, cols) vUS_2D_quad_batchModel(X, tau(cols), k0, sigma, s, w);
lsqo = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
nq = 40;

% 5a. all six parameters free
allCols = (1:nTau).';
g1 = (model(truth.', allCols)).' + noise();
x0 = [Vx0, Vz0, ones(nP,1), zeros(nP,1), 2.5*ones(nP,1), ones(nP,1)];
lb = [zeros(nP,1), Vz0 - 0.01, zeros(nP,2), 2*ones(nP,1), zeros(nP,1)];  ub = [30e-3*ones(nP,1), Vz0 + 0.01, ones(nP,2), 3*ones(nP,1), ones(nP,1)];
o = struct('window', tdi, 't1i', t1i, 'xscale', xscale);
o.jacobian = 'analytic';  [xA, cA, ~, cvA] = fitBatchedLM(model, x0, g1, lb, ub, o);
o.jacobian = 'forward';   [xF, cF, ~, cvF] = fitBatchedLM(model, x0, g1, lb, ub, o);
% the wiring check: the ORIGINAL specialised solver on the same data through the 3D functions (sigma_y = sigma_x) must give the same fit
[xO, cO] = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma3, [0 NaN 0 0 2 0], [30e-3 NaN 1 1 3 1], struct('x0', [5e-3 NaN 1 0 2.5 1], 'xscale', xscale));
cL = zeros(nq,1);
for q = 1:nq
    inds = t1i:tdi(q); gv = g1(q,:).';
    [~, cL(q)] = lsqnonlin(@(xs) vUS_3D_quad_residJac(xs, tau(inds), k0, sigma3, s, w, real(gv(inds)), imag(gv(inds)), xscale), x0(q,:)./xscale, lb(q,:)./xscale, ub(q,:)./xscale, ...
        optimoptions(lsqo, 'SpecifyObjectiveGradient', true));
end
rA = cA(1:nq)./cL;  rAF = cF./cA;
ok = max(abs(xA(:) - xO(:))) < 1e-12 && mean(cvA) >= 0.99 && median(rA) < 1 + 1e-6 && median(abs(rAF - 1)) < 1e-6 && all(xA(:,3) >= 0 & xA(:,3) <= 1) && all(xA(:,1) <= 30e-3 + 1e-12);
allOk = report('5a. six free parameters, %d voxels: identical to the original vUS_3D_quad_fitBatched (max|dx| = %.1e); converged %.1f%%; forward-difference vs analytic cost: median |ratio - 1| = %.1e; median |v_zgp - truth| = %.2f mm/s', ...
    ok, nP, max(abs(xA(:) - xO(:))), 100*mean(cvA), median(abs(rAF - 1)), 1e3*median(abs(xA(:,2) - truth(:,2)))) && allOk;
fprintf('    INFO (not a pass/fail criterion): cost/lsqnonlin on %d voxels: median %.8f, LM ends >1%% worse in %d (%.0f%%) -- with the uniform Vx0 = 5 mm/s start most of these end exactly on lb(v_xgp) = 0; see the header of vUS_2D_quad_batchModel.m\n', nq, median(rA), sum(rA > 1.01), 100*mean(rA > 1.01));

% 5b. the script's current 4-parameter fit: k = 2, a = 1 fixed (lb = ub), vs lsqnonlin on the erf model the script uses
tr4 = truth;  tr4(:,5) = 2;  tr4(:,6) = 1;
g14 = (model(tr4.', allCols)).' + noise();
x04 = x0;  x04(:,5:6) = repmat([2 1], nP, 1);
lb4 = lb;  lb4(:,5:6) = repmat([2 1], nP, 1);  ub4 = ub;  ub4(:,5:6) = repmat([2 1], nP, 1);   % the script's lb/ub for the first four, k and a fixed
o4 = o;  o4.jacobian = 'forward';
[x4f, c4, ~, cv4] = fitBatchedLM(model, x04, g14, lb4, ub4, o4);
cL4 = zeros(nq,1);
for q = 1:nq
    inds = t1i:tdi(q); gv = g14(q,:).';
    if hasErf, fun = @(p) vUS_2D_erf_vec(p, tau(inds), k0, sigma) - gv(inds);
    else,      fun = @(p) vUS_2D_quad_batchModel([p(:); 2; 1], tau(inds), k0, sigma, s, w) - gv(inds);
    end
    [~, cL4(q)] = lsqnonlin(@(p) [real(fun(p)); imag(fun(p))], x04(q,1:4), lb4(q,1:4), ub4(q,1:4), lsqo);
end
r4 = c4(1:nq)./cL4;
ok = all(x4f(:,5) == 2) && all(x4f(:,6) == 1) && mean(cv4) >= 0.99 && median(r4) < 1 + 1e-6 && sum(r4 > 1.01) <= 0.05*nq;
src = 'vUS_2D_quad_batchModel'; if hasErf, src = 'vUS_2D_erf_vec (the script''s model)'; end
allOk = report('5b. k = 2, a = 1 fixed: k, a returned exactly = %d; converged %.1f%%; cost/lsqnonlin on %s (%d voxels): median %.8f, worse by >1%% in %d', ...
    ok, all(x4f(:,5) == 2) && all(x4f(:,6) == 1), 100*mean(cv4), src, nq, median(r4), sum(r4 > 1.01)) && allOk;

if allOk, fprintf('\nAll checks passed.\n'); else, warning('One or more checks failed -- inspect before use.'); end
for i = 1:numel(addPaths)
    if isfolder(addPaths{i}), rmpath(addPaths{i}); end
end

function tf = report(fmt, ok, varargin)
% print one check line with PASS/FAIL; ok is the verdict, varargin the values for the format
    if ok, s = 'PASS'; else, s = 'FAIL'; end
    fprintf([fmt ' -> %s\n'], varargin{:}, s);
    tf = logical(ok);
end

function tf = throws(fh, idSuffix)
% true if fh() errors with an identifier ending in idSuffix
    try
        fh(); tf = false;
    catch e
        tf = endsWith(e.identifier, idSuffix);
    end
end
