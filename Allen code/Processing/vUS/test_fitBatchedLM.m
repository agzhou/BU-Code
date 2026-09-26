%% Description:
%   Self-contained check of the generic batched solver (no data files needed):
%     1. vUS_3D_quad_batchModel.m (all voxels at once) vs. the per-voxel vUS_3D_quad_vec.m /
%        vUS_3D_quad_complex_Jac.m -- value and analytic Jacobian
%     2. fitBatchedLM.m with the vUS model reproduces vUS_3D_quad_fitBatched.m (analytic Jacobian: exactly;
%        forward/central finite differences: same quality), and is independent of the chunk size
%     3. models the solver has never seen: a real-valued and a complex-valued model, each vs. a per-voxel
%        lsqnonlin loop, with finite-difference and with analytic Jacobians and per-voxel start points
%     4. bounds, fixed parameters, per-voxel start points, fit windows (last index / mask / holes), weights,
%        model handles with fewer inputs
%     5. bad input is rejected with clear errors; non-finite models and empty windows do not crash the batch
%     6. batchModelFromSingle.m (per-voxel model through the batched solver), single precision
%     7. timing (informational only)
%   Prints PASS/FAIL per check. Requires the Optimization Toolbox (for lsqnonlin).

rng(1);
sigma = [368.1124 368.1124 126.4505]*1e-6;     % RC15gV PSF (5 x 2 angles, -6 to 6 deg)
k0 = 2*pi/0.000113235294117647;                % 13.6 MHz, c = 1540 m/s
frameRate = 4900; nTau = 60; tau = (0:nTau-1)'/frameRate; t1i = 2;
[s, w] = gaussLegendre01(48); xscale = [1e-2 1e-2 1 1 1 1];
allOk = true;

%% Synthetic vUS voxels: model + complex noise, random fit windows
nV = 200; win = randi([8 40], nV, 1); tdi = t1i + win - 1;
truth = [(2 + 18*rand(nV,1))*1e-3, (-25 + 50*rand(nV,1))*1e-3, 0.5 + 0.35*rand(nV,1), 0.15*rand(nV,1), 2 + rand(nV,1), 0.3 + 0.7*rand(nV,1)];
g1 = zeros(nV, nTau);
for i = 1:nV
    g1(i,:) = (vUS_3D_quad_vec(truth(i,:), tau, k0, sigma, s, w) + 0.03*(randn(nTau,1) + 1i*randn(nTau,1))/sqrt(2)).';
end
Vz0 = truth(:,2) + 3e-3*randn(nV,1);
x0v = [5e-3*ones(nV,1), Vz0, ones(nV,1), zeros(nV,1), 2.5*ones(nV,1), ones(nV,1)];                      % same start as vUS_3D_quad_fitBatched
lbv = [zeros(nV,1), Vz0 - 0.01, zeros(nV,1), zeros(nV,1), 2*ones(nV,1), zeros(nV,1)];                    % ... and its default bounds
ubv = [sqrt(2)*30e-3*ones(nV,1), Vz0 + 0.01, ones(nV,1), ones(nV,1), 3*ones(nV,1), ones(nV,1)];
vModel = @(X, cols) vUS_3D_quad_batchModel(X, tau(cols), k0, sigma, s, w);
vOpts = struct('window', tdi, 't1i', t1i, 'xscale', xscale);

%% 1. batched vUS model + Jacobian vs. the per-voxel reference code
Xt = [truth(:,1:2), 0.7*ones(nV,1), 0.1*ones(nV,1), 2.4*ones(nV,1), 0.6*ones(nV,1)].';   % generic interior point (physical units)
[Y, J] = vUS_3D_quad_batchModel(Xt, tau, k0, sigma, s, w);  Y1 = vUS_3D_quad_batchModel(Xt, tau, k0, sigma, s, w);
dY = 0; dJ = 0; Jmax = 0;
for q = 1:nV
    [gq, Jq] = vUS_3D_quad_complex_Jac(Xt(:,q).', tau, k0, sigma, s, w);
    dY = max(dY, max(abs(Y(:,q) - gq)));  dJ = max(dJ, max(abs(J(:,:,q) - Jq), [], 'all'));  Jmax = max(Jmax, max(abs(Jq), [], 'all'));
end
ok = dY < 1e-12 && dJ/Jmax < 1e-12 && max(abs(Y1 - Y), [], 'all') < 1e-12;
allOk = report('1. vUS batch model/Jacobian vs vUS_3D_quad_vec / complex_Jac: max|dY| = %.1e, max|dJ|/max|J| = %.1e (single-output call agrees: %d)', ok, dY, dJ/Jmax, max(abs(Y1 - Y), [], 'all') < 1e-12) && allOk;

%% 2. generic solver + vUS model vs. vUS_3D_quad_fitBatched
[xo, co, io] = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma);
o = vOpts;  o.jacobian = 'analytic';
[xa, ca, ia, cva] = fitBatchedLM(vModel, x0v, g1, lbv, ubv, o);
ok = max(abs(xa(:) - xo(:))) < 1e-12 && isequal(ia, io);
allOk = report('2a. analytic Jacobian reproduces vUS_3D_quad_fitBatched: max|dx| = %.1e, identical iteration counts = %d', ok, max(abs(xa(:) - xo(:))), isequal(ia, io)) && allOk;

for mode = {'forward', 'central'}
    o = vOpts;  o.jacobian = mode{1};
    [xf, cf, ifw, cvf] = fitBatchedLM(vModel, x0v, g1, lbv, ubv, o);
    r = cf./co;
    ok = mean(cvf) >= 0.99 && abs(median(r) - 1) < 1e-6 && sum(r > 1.01) <= 0.05*nV;
    allOk = report('2b. %s finite differences (no Jacobian supplied): converged %.1f%%, cost/old cost: median %.8f, worse by >1%% in %d of %d voxels, better in %d; median iterations %g (old %g)', ...
        ok, mode{1}, 100*mean(cvf), median(r), sum(r > 1.01), nV, sum(r < 0.99), median(ifw), median(io)) && allOk;
end

o = vOpts;  o.jacobian = 'analytic';  o.chunkSize = 37;
[xc, cc] = fitBatchedLM(vModel, x0v, g1, lbv, ubv, o);
ok = max(abs(cc - ca)./ca) < 1e-9;
allOk = report('2c. chunkSize = 37 vs 1024 (voxels are independent): max relative cost difference = %.1e', ok, max(abs(cc - ca)./ca)) && allOk;

%% 3. models the solver has never seen, vs. a per-voxel lsqnonlin loop
lsqo = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
% 3a. real-valued: y(t) = A exp(-t/T) + C
nE = 300; nt = 40; tt = (0:nt-1).'/nt;
tru = [0.5 + 2*rand(nE,1), 0.1 + 0.9*rand(nE,1), 0.3*(2*rand(nE,1) - 1)];
Ye = tru(:,1).*exp(-tt.'./tru(:,2)) + tru(:,3) + 0.02*randn(nE, nt);
mExp = @(X, cols) X(1,:).*exp(-tt(cols)./X(2,:)) + X(3,:);
lbE = [0 0.02 -1];  ubE = [5 5 1];  x0E = [1 0.5 0];
[xE, cE, ~, cvE] = fitBatchedLM(mExp, x0E, Ye, lbE, ubE);          % defaults: forward differences, automatic xscale
[xEa, cEa] = fitBatchedLM(@expModelJ, x0E, Ye, lbE, ubE, struct('jacobian', 'analytic'));
cL = zeros(nE,1); xL = zeros(nE,3);
for q = 1:nE
    [xL(q,:), cL(q)] = lsqnonlin(@(p) p(1)*exp(-tt/p(2)) + p(3) - Ye(q,:).', x0E, lbE, ubE, lsqo);
end
rE = cE./cL;  rEa = cEa./cE;
ok = mean(cvE) >= 0.99 && median(rE) < 1 + 1e-6 && sum(rE > 1.01) <= 0.05*nE;
allOk = report('3a. real exponential model, finite differences: converged %.1f%%, cost/lsqnonlin cost: median %.8f, worse by >1%% in %d of %d voxels; median |A - truth|: %.3f (lsqnonlin %.3f)', ...
    ok, 100*mean(cvE), median(rE), sum(rE > 1.01), nE, median(abs(xE(:,1) - tru(:,1))), median(abs(xL(:,1) - tru(:,1)))) && allOk;
ok = median(abs(rEa - 1)) < 1e-6 && max(abs(rEa - 1)) < 1e-2;
allOk = report('3b. same model with an analytic Jacobian agrees with finite differences: median |cost ratio - 1| = %.1e, max = %.1e', ok, median(abs(rEa - 1)), max(abs(rEa - 1))) && allOk;

% 3c. complex-valued with a per-voxel start point: y(t) = A exp((-r + 2 pi i f) t) + DC
nC = 200; tC = (0:29).'/100;
trC = [0.6 + 0.4*rand(nC,1), 1 + 3*rand(nC,1), -3 + 6*rand(nC,1), 0.05*rand(nC,1)];
Yc = trC(:,1).*exp(complex(-trC(:,2).*tC.', 2*pi*trC(:,3).*tC.')) + trC(:,4) + 0.02*(randn(nC,30) + 1i*randn(nC,30))/sqrt(2);
mCx = @(X, cols) X(1,:).*exp(complex(-X(2,:).*tC(cols), 2*pi*X(3,:).*tC(cols))) + X(4,:);
x0C = [0.8*ones(nC,1), 2*ones(nC,1), trC(:,3) + 0.3*randn(nC,1), zeros(nC,1)];
lbC = [0 0 -10 0];  ubC = [2 10 10 1];
[xC, cC, ~, cvC] = fitBatchedLM(mCx, x0C, Yc, lbC, ubC);
cLc = zeros(nC,1);
for q = 1:nC
    cpx = @(p) p(1)*exp(complex(-p(2)*tC, 2*pi*p(3)*tC)) + p(4) - Yc(q,:).';
    [~, cLc(q)] = lsqnonlin(@(p) [real(cpx(p)); imag(cpx(p))], x0C(q,:), lbC, ubC, lsqo);
end
rC = cC./cLc;
ok = mean(cvC) >= 0.99 && median(rC) < 1 + 1e-6 && sum(rC > 1.01) <= 0.05*nC;
allOk = report('3c. complex damped-oscillation model, per-voxel x0: converged %.1f%%, cost/lsqnonlin cost: median %.8f, worse by >1%% in %d of %d voxels; median |f - truth|: %.3f Hz', ...
    ok, 100*mean(cvC), median(rC), sum(rC > 1.01), nC, median(abs(xC(:,3) - trC(:,3)))) && allOk;

%% 4. bounds, start points, windows, weights, model arity (real exponential model)
[xB0, cB0] = fitBatchedLM(mExp, x0E, Ye, lbE, ubE);
% 4a. fixed parameter (lb = ub), partial NaN override, per-voxel bounds
lbF = [NaN NaN 0.1];  ubF = [NaN NaN 0.1];                      % C fixed at 0.1; A, T unbounded
xF = fitBatchedLM(mExp, x0E, Ye, lbF, ubF);
lbP = nan(1,3); lbP(2) = 0.5;                                    % only T bounded below
xP = fitBatchedLM(mExp, x0E, Ye, lbP, []);
hw = 0.05 + 0.2*rand(nE,1);  lbV = nan(nE,3);  ubV = nan(nE,3);  lbV(:,2) = tru(:,2) - hw;  ubV(:,2) = tru(:,2) + hw;  ubV(:,1) = 1.5;
xV = fitBatchedLM(mExp, x0E, Ye, lbV, ubV);
ok = all(abs(xF(:,3) - 0.1) < 1e-14) && all(xP(:,2) >= 0.5 - 1e-12) && all(xV(:,2) >= tru(:,2) - hw - 1e-12 & xV(:,2) <= tru(:,2) + hw + 1e-12) && all(xV(:,1) <= 1.5 + 1e-12);
allOk = report('4a. fixed parameter (lb = ub): all C == 0.1 is %d; partial override lb(T) = 0.5: min T = %.4f; per-voxel bounds and a ub(A) cap honored', ok, all(abs(xF(:,3) - 0.1) < 1e-14), min(xP(:,2))) && allOk;

% 4b. [1,nP] start point == the same start repeated per voxel; window as last-index vector == the equivalent logical mask
xR = fitBatchedLM(mExp, repmat(x0E, nE, 1), Ye, lbE, ubE);
lastE = randi([15 nt], nE, 1);
[xW1, cW1] = fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('window', lastE, 't1i', 3));
maskE = (1:nt) >= 3 & (1:nt) <= lastE;
[xW2, cW2] = fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('window', maskE));
ok = isequal(xR, xB0) && isequal(xW1, xW2) && isequal(cW1, cW2);
allOk = report('4b. [1,nP] vs per-voxel x0 identical = %d; window as last-index vector (t1i = 3) vs equivalent logical mask identical = %d', ok, isequal(xR, xB0), isequal(xW1, xW2)) && allOk;

% 4c. a mask with holes == zero weights at the holes; weights as [1,nT], [nT,1] and [nVox,nT]
maskH = maskE;  maskH(:, 20:22) = false;  maskH(1:50, 8) = false;
wZero = double(maskH);  wUnit = ones(1, nt);
[xH1, cH1] = fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('window', maskH));
[xH2, cH2] = fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('window', maskE, 'weights', wZero));
[xU1, cU1] = fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('weights', wUnit));
[xU2, cU2] = fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('weights', wUnit.'));
[xU3, cU3] = fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('weights', repmat(wUnit, nE, 1)));
ok = isequal(xH1, xH2) && isequal(cH1, cH2) && isequal(xU1, xB0) && isequal(xU2, xB0) && isequal(xU3, xB0);
allOk = report('4c. logical mask with holes == zero weights at the holes: %d; unit weights as [1,nT], [nT,1], [nVox,nT] == no weights: %d', ok, isequal(xH1, xH2) && isequal(cH1, cH2), isequal(xU1, xB0) && isequal(xU2, xB0) && isequal(xU3, xB0)) && allOk;

% 4d. non-uniform weights vs. a weighted per-voxel lsqnonlin, and the weighted cost is what is reported
wLin = linspace(1, 0.1, nt);
[xWt, cWt] = fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('weights', wLin));
nq = 60; cLw = zeros(nq,1);
for q = 1:nq
    [~, cLw(q)] = lsqnonlin(@(p) wLin(:).*(p(1)*exp(-tt/p(2)) + p(3) - Ye(q,:).'), x0E, lbE, ubE, lsqo);
end
rw = cWt(1:nq)./cLw;
ok = median(rw) < 1 + 1e-6 && sum(rw > 1.01) <= 0.05*nq;
allOk = report('4d. residual weights vs weighted lsqnonlin (%d voxels): cost ratio median %.8f, worse by >1%% in %d', ok, nq, median(rw), sum(rw > 1.01)) && allOk;

% 4e. model handles with fewer inputs: @(X) (full window only) and @(X, cols)
[xM1, cM1] = fitBatchedLM(@(X) X(1,:).*exp(-tt./X(2,:)) + X(3,:), x0E, Ye, lbE, ubE);
ok = isequal(xM1, xB0) && isequal(cM1, cB0);
allOk = report('4e. model handle @(X) with a single input gives the same fit as @(X, cols): %d', ok, ok) && allOk;

%% 5. bad input and non-finite models
ok5 = true;
ok5 = ok5 && throws(@() fitBatchedLM(mExp, x0E, Ye, ubE, lbE), 'badBounds');                                              % lb > ub
ok5 = ok5 && throws(@() fitBatchedLM(mExp, x0E, Ye, ones(7,3), []), 'badBounds');                                        % wrong-size bounds
ok5 = ok5 && throws(@() fitBatchedLM(mExp, x0E, Ye, struct('tol', 1e-3)), 'optsMoved');                                   % opts where lb belongs
ok5 = ok5 && throws(@() fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('nNodes', 24)), 'badOption');                        % unknown option
ok5 = ok5 && throws(@() fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('stepCap', [1 2])), 'badStepCap');                   % wrong-size stepCap
ok5 = ok5 && throws(@() fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('xscale', [1 2])), 'badXscale');                     % wrong-size xscale
ok5 = ok5 && throws(@() fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('weights', -ones(1, nt))), 'badWeights');            % negative weights
ok5 = ok5 && throws(@() fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('window', true(3, 3))), 'badWindow');                % wrong-size mask
ok5 = ok5 && throws(@() fitBatchedLM(mExp, [1 NaN 0], Ye, lbE, ubE), 'badX0');                                            % non-finite x0
ok5 = ok5 && throws(@() fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('jacobian', 'magic')), 'MATLAB:unrecognizedStringChoice'); % unknown Jacobian mode
ok5 = ok5 && throws(@() fitBatchedLM(@(X, cols) X(1:2,:), x0E, Ye, lbE, ubE), 'badModelOutput');                         % wrong-size model output
ok5 = ok5 && throws(@() fitBatchedLM(@(X, cols) X(1,:).*exp(1i*tt(cols)), x0E, Ye, lbE, ubE), 'complexModel');            % complex model, real data
ok5 = ok5 && throws(@() fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('jacobian', 'analytic')), 'noJacobian');           % analytic requested, model has one output
allOk = report('5a. lb > ub, wrong-size bounds/stepCap/xscale/mask, opts as lb, unknown option, negative weights, NaN x0, unknown Jacobian mode, wrong-size or complex model output are rejected', ok5) && allOk;

% 5b. a model that is not finite at some start points: those voxels are skipped, the rest is fitted, no crash
x0N = repmat(x0E, nE, 1);  x0N(1:5, 1) = 0.01;                                          % A = 0.01 -> model returns NaN there
mNaN = @(X, cols) mExp(X, cols) + 0./(X(1,:) > 0.05);                                    % 0/0 = NaN where A <= 0.05 (start only)
lastwarn('', '');
[xN, cN, ~, cvN] = fitBatchedLM(mNaN, x0N, Ye, lbE, ubE);
[~, wid] = lastwarn;
ok = strcmp(wid, 'fitBatchedLM:badStart') && all(isnan(cN(1:5))) && ~any(cvN(1:5)) && max(abs(xN(1:5,:) - x0N(1:5,:)), [], 'all') < 1e-12 && all(isfinite(cN(6:end))) && mean(cvN(6:end)) >= 0.99;
allOk = report('5b. model NaN at the start of 5 voxels: warns = %d, those return x0 with cost NaN, the other %d voxels are fitted (converged %.1f%%)', ok, strcmp(wid, 'fitBatchedLM:badStart'), nE - 5, 100*mean(cvN(6:end))) && allOk;

% 5c. empty fit windows: skipped with a warning
lastE0 = lastE;  lastE0(7) = 0;
lastwarn('', '');
[xZ, cZ, ~, cvZ] = fitBatchedLM(mExp, x0E, Ye, lbE, ubE, struct('window', lastE0, 't1i', 3));
[~, wid] = lastwarn;
others = [1:6 8:nE];  sameOthers = max(abs(xZ(others,:) - xW1(others,:)), [], 'all') < 1e-9;
ok = strcmp(wid, 'fitBatchedLM:emptyWindow') && isnan(cZ(7)) && ~cvZ(7) && isequal(xZ(7,:), x0E) && sameOthers;
allOk = report('5c. voxel with an empty window: warns = %d, returns x0, the others are unaffected = %d', ok, strcmp(wid, 'fitBatchedLM:emptyWindow'), sameOthers) && allOk;

%% 6. per-voxel model through the adapter; single precision
nS = 20;  idS = 1:nS;
perVox   = @(x, cols) vUS_3D_quad_vec(x, tau(cols), k0, sigma, s, w);
perVoxJ  = @(x, cols) vUS_3D_quad_complex_Jac(x, tau(cols), k0, sigma, s, w);
oS = struct('window', tdi(idS), 't1i', t1i, 'xscale', xscale);
[xS0, cS0] = fitBatchedLM(vModel, x0v(idS,:), g1(idS,:), lbv(idS,:), ubv(idS,:), oS);
[xS1, cS1] = fitBatchedLM(batchModelFromSingle(perVox), x0v(idS,:), g1(idS,:), lbv(idS,:), ubv(idS,:), oS);
oS.jacobian = 'analytic';
[xS2, cS2] = fitBatchedLM(batchModelFromSingle(perVoxJ), x0v(idS,:), g1(idS,:), lbv(idS,:), ubv(idS,:), oS);
[xS3, cS3] = fitBatchedLM(vModel, x0v(idS,:), g1(idS,:), lbv(idS,:), ubv(idS,:), oS);
ok = max(abs(cS1./cS0 - 1)) < 1e-6 && max(abs(cS2./cS3 - 1)) < 1e-9;
allOk = report('6a. per-voxel model via batchModelFromSingle (%d voxels): finite-difference cost vs batched model max |ratio - 1| = %.1e; per-voxel analytic Jacobian vs batched analytic: %.1e', ok, nS, max(abs(cS1./cS0 - 1)), max(abs(cS2./cS3 - 1))) && allOk;

o = vOpts;  o.jacobian = 'analytic';  o.precision = 'single';
[xSg, cSg, ~, cvSg] = fitBatchedLM(vModel, x0v, g1, lbv, ubv, o);
rSg = cSg./ca;
ok = all(isfinite(xSg(:))) && mean(cvSg) >= 0.95 && median(rSg) < 1.001 && sum(rSg > 1.01) <= 0.10*nV;
allOk = report('6b. single precision: converged %.1f%%, cost vs double: median %.6f, worse by >1%% in %d of %d voxels', ok, 100*mean(cvSg), median(rSg), sum(rSg > 1.01), nV) && allOk;

if allOk, fprintf('\nAll checks passed.\n'); else, warning('One or more checks failed -- inspect before use.'); end

%% 7. timing (informational)
nBig = 2000;  ib = mod((0:nBig-1).', nV) + 1;
gB = g1(ib,:);  tdB = tdi(ib);  x0B = x0v(ib,:);  lbB = lbv(ib,:);  ubB = ubv(ib,:);
oT = struct('window', tdB, 't1i', t1i, 'xscale', xscale);
oW = oT;  oW.window = tdB(1:50);
fitBatchedLM(vModel, x0B(1:50,:), gB(1:50,:), lbB(1:50,:), ubB(1:50,:), oW);   % warm-up
tic; vUS_3D_quad_fitBatched(gB, tdB, Vz0(ib), tau, k0, sigma); tOld = toc;
fprintf('\nTiming on %d voxels (vUS model, %d threads):\n  vUS_3D_quad_fitBatched (specialised)          %6.2f s\n', nBig, maxNumCompThreads, tOld);
for mode = {'analytic', 'forward', 'central'}
    oT.jacobian = mode{1};
    tic; fitBatchedLM(vModel, x0B, gB, lbB, ubB, oT); tFit = toc;
    fprintf('  fitBatchedLM, %-8s Jacobian               %6.2f s   (%.2fx the specialised solver)\n', mode{1}, tFit, tFit/tOld);
end
nL = 100; tic
for q = 1:nL
    inds = t1i:tdB(q); gv = gB(q,:).';
    lsqnonlin(@(xs) vUS_3D_quad_residJac(xs, tau(inds), k0, sigma, s, w, real(gv(inds)), imag(gv(inds)), xscale), x0B(q,:)./xscale, lbB(q,:)./xscale, ubB(q,:)./xscale, ...
        optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true, 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10));
end
tL = toc*nBig/nL;
fprintf('  serial lsqnonlin loop (analytic Jacobian)      %6.2f s   (extrapolated from %d voxels)\n', tL, nL);

function tf = report(fmt, ok, varargin)
% print one check line with PASS/FAIL; ok may be a logical or (when varargin holds the values) the verdict
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

function [Y, J] = expModelJ(X, cols)
% y(t) = A exp(-t/T) + C with its analytic Jacobian (test model for the analytic path)
    tt = (0:39).'/40;  t = tt(cols);
    A = X(1,:);  T = X(2,:);  C = X(3,:);
    e = exp(-t./T);
    Y = A.*e + C;
    if nargout > 1, J = permute(cat(3, e, A.*e.*t./T.^2, ones(size(e))), [1 3 2]); end
end
