%% Description:
%   Self-contained check of the batched fitting code (no data files needed):
%     1. vUS_3D_quad_batchEval.m (all voxels at once) vs. the per-voxel, symbolically
%        derived vUS_3D_quad_residJac.m -- residual and Jacobian must agree to round-off
%     2. vUS_3D_quad_fitBatched.m vs. a per-voxel lsqnonlin loop on synthetic voxels
%        generated from the model plus noise: final cost, and recovery of v_zgp
%   Prints PASS/FAIL per check. Requires the Optimization Toolbox (for lsqnonlin).

rng(1);
sigma = [368.1124 368.1124 126.4505]*1e-6;     % RC15gV PSF (5 x 2 angles, -6 to 6 deg)
k0 = 2*pi/0.000113235294117647;                % 13.6 MHz, c = 1540 m/s
frameRate = 4900; nTau = 60; tau = (0:nTau-1)'/frameRate; t1i = 2;
[s, w] = gaussLegendre01(48); xscale = [1e-2 1e-2 1 1 1 1];

%% Synthetic voxels: model + complex noise, random fit windows
nV = 200; win = randi([8 40], nV, 1); tdi = t1i + win - 1;
truth = [(2 + 18*rand(nV,1))*1e-3, (-25 + 50*rand(nV,1))*1e-3, 0.5 + 0.35*rand(nV,1), 0.15*rand(nV,1), 2 + rand(nV,1), 0.3 + 0.7*rand(nV,1)];
g1 = zeros(nV, nTau);
for i = 1:nV
    g1(i,:) = (vUS_3D_quad_vec(truth(i,:), tau, k0, sigma, s, w) + 0.03*(randn(nTau,1) + 1i*randn(nTau,1))/sqrt(2)).';
end
Vz0 = truth(:,2) + 3e-3*randn(nV,1);

%% 1. batched model + Jacobian vs. per-voxel residJac
c = (1:nV).'; P = vUS_3D_quad_batchPack(c, g1, tdi, t1i, tau, sigma, k0, s, w, xscale, 'double');
Xt = [truth(:,1:2)./xscale(1:2), 0.7*ones(nV,1), 0.1*ones(nV,1), 2.4*ones(nV,1), 0.6*ones(nV,1)].';   % generic interior point (scaled)
[r, J] = vUS_3D_quad_batchEval(Xt, 1:nV, P); L = size(P.YR, 1);
maxdr = 0; maxdJ = 0;
for q = 1:nV
    inds = t1i:tdi(q); gv = g1(q,:).'; m = numel(inds);
    [r2, J2] = vUS_3D_quad_residJac(Xt(:,q).', tau(inds), k0, sigma, s, w, real(gv(inds)), imag(gv(inds)), xscale);
    ri = [1:m, L + (1:m)]; rq = r(:,:,q); Jq = J(:,:,q);
    maxdr = max(maxdr, max(abs(rq(ri) - r2)));  maxdJ = max(maxdJ, max(abs(Jq(ri,:) - J2), [], 'all'));
end
ok1 = maxdr < 1e-12 && maxdJ < 1e-9;
fprintf('1. batched model/Jacobian vs vUS_3D_quad_residJac: max|dr| = %.1e, max|dJ| = %.1e -> %s\n', maxdr, maxdJ, passfail(ok1));

%% 2. batched fit vs. per-voxel lsqnonlin (analytic Jacobian, scaled, tight tolerance)
[x, cost, iters, conv] = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma);
o = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true, 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
cL = zeros(nV,1); xL = zeros(nV,6);
for q = 1:nV
    inds = t1i:tdi(q); gv = g1(q,:).';
    x0 = [5e-3, Vz0(q), 1, 0, 2.5, 1]; lb = [0, Vz0(q)-0.01, 0, 0, 2, 0]; ub = [sqrt(2)*30e-3, Vz0(q)+0.01, 1, 1, 3, 1];
    [xs, cL(q)] = lsqnonlin(@(xs) vUS_3D_quad_residJac(xs, tau(inds), k0, sigma, s, w, real(gv(inds)), imag(gv(inds)), xscale), x0./xscale, lb./xscale, ub./xscale, o);
    xL(q,:) = xs.*xscale;
end
ratio = cost./cL;
fprintf('2. batched fit: converged %.1f%% | cost/lsqnonlin cost: median %.6f, batched worse by >1%% in %d of %d voxels, better by >1%% in %d\n', 100*mean(conv), median(ratio), sum(ratio > 1.01), nV, sum(ratio < 0.99));
fprintf('   median |v_zgp - truth|: batched %.2f mm/s, lsqnonlin %.2f mm/s | median iterations: batched %d\n', 1e3*median(abs(x(:,2) - truth(:,2))), 1e3*median(abs(xL(:,2) - truth(:,2))), median(iters));
ok2 = mean(conv) >= 0.99 && median(ratio) < 1 + 1e-6 && sum(ratio > 1.01) <= 0.10*nV;
fprintf('   -> %s (needs: >=99%% converged, median cost ratio ~1, <=10%% of voxels >1%% worse)\n', passfail(ok2));

%% 3. lb / ub inputs (vUS_3D_quad_fitBatched(..., sigma, lb, ub))
% 3a. passing the DEFAULT bounds explicitly ([nVox,6] matrices) must reproduce the default call exactly
lbD = [zeros(nV,1), Vz0 - 0.01, zeros(nV,1), zeros(nV,1), 2*ones(nV,1), zeros(nV,1)];
ubD = [sqrt(2)*30e-3*ones(nV,1), Vz0 + 0.01, ones(nV,1), ones(nV,1), 3*ones(nV,1), ones(nV,1)];
xE = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, lbD, ubD);
ok3a = max(abs(xE(:) - x(:))) == 0;
fprintf('3a. explicit default bounds reproduce the default fit: max|dx| = %.1e -> %s\n', max(abs(xE(:) - x(:))), passfail(ok3a));

% 3b. the bounds used in vUS_3D_newmodel.m's serial loop: -Inf/Inf entries, one [1,6] row for every voxel
lbI = [0, -Inf, 0, 0, 2, 0];  ubI = [Inf, Inf, 1, 1, 3, 1];
[xI, costI, ~, convI] = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, lbI, ubI);
inI = all(xI >= lbI - 1e-12 & xI <= ubI + 1e-12, 'all') && all(isfinite(xI(:)));
rI = costI./cost;
ok3b = inI && mean(convI) >= 0.99 && sum(rI > 1.01) <= 0.10*nV;
fprintf('3b. infinite bounds [0 -Inf 0 0 2 0]/[Inf Inf 1 1 3 1]: finite and inside bounds = %d, converged %.1f%%, cost vs default-bounds fit: median %.6f, >1%% worse in %d voxels -> %s\n', inI, 100*mean(convI), median(rI), sum(rI > 1.01), passfail(ok3b));

% 3c. fixing a parameter (lb = ub) and a partial NaN override (only k bounded below; everything else default)
lbF = nan(1,6); ubF = nan(1,6); lbF(5) = 2; ubF(5) = 2;
xF = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, lbF, ubF);
lbP = nan(1,6); lbP(5) = 2.2;
xP = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, lbP, []);
ok3c = all(xF(:,5) == 2) && all(xP(:,5) >= 2.2 - 1e-12) && all(xP(:,2) >= Vz0 - 0.01 - 1e-12 & xP(:,2) <= Vz0 + 0.01 + 1e-12);
fprintf('3c. k fixed by lb = ub = 2: all k == 2 is %d; partial override lb(k) = 2.2: min k = %.4f, v_zgp still inside its default window -> %s\n', all(xF(:,5) == 2), min(xP(:,5)), passfail(ok3c));

% 3d. per-voxel bounds honored ([nVox,6]): random v_zgp window half-width per voxel, F capped at 0.9
hw = 0.001 + 0.004*rand(nV,1);
lbV = nan(nV,6);  lbV(:,2) = Vz0 - hw;
ubV = nan(nV,6);  ubV(:,2) = Vz0 + hw;  ubV(:,3) = 0.9;
xV =vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, lbV, ubV);
ok3d = all(xV(:,2) >= Vz0 - hw - 1e-12 & xV(:,2) <= Vz0 + hw + 1e-12) && all(xV(:,3) <= 0.9 + 1e-12);
fprintf('3d. per-voxel bounds honored (v_zgp window %.1f-%.1f mm, F <= 0.9) -> %s\n', 1e3*min(hw), 1e3*max(hw), passfail(ok3d));

% 3e. bad input is rejected with clear errors
ok3e = throws(@() vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, ubD, lbD), 'badBounds') ...              % lb > ub
    && throws(@() vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, ones(3,6), []), 'badBounds') ...          % wrong size
    && throws(@() vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, struct('nNodes', 24)), 'optsMoved') ...   % old call style (opts as 7th input)
    && throws(@() vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, [], [], struct('stepCap', [1 2 3])), 'badStepCap');   % wrong-size stepCap
fprintf('3e. lb > ub, wrong-size bounds, the old opts-as-7th-input call and a wrong-size stepCap are rejected -> %s\n', passfail(ok3e));

% 3f. per-parameter stepCap is accepted; with an unbounded v_tgp the automatic default is the same as passing it explicitly
xA = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, lbI, ubI);
xC = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, lbI, ubI, struct('stepCap', [10 0.5 0.5 0.5 0.5 0.5]));
ok3f = all(isfinite(xC(:))) && max(abs(xA(:) - xC(:))) == 0;
fprintf('3f. automatic stepCap for an unbounded v_tgp equals stepCap = [10 .5 .5 .5 .5 .5]: max|dx| = %.1e -> %s\n', max(abs(xA(:) - xC(:))), passfail(ok3f));

% 3g. REGRESSION: a huge finite ub(v_tgp) (50e3 m/s, a slip for 50e-3) once pushed every start point to 0.1% of the box
% (v_tgp = 50 m/s), where the model is flat: 100% of voxels hit the iteration cap and the fit was ~20x worse.
% The start-point nudge is now capped, so it must behave like an unbounded ub(v_tgp), and it must warn about the units.
lbW = [0 NaN 0 0 2 0];
lastwarn('', '');
[xW, costW, itW, convW] = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, lbW, [50e3 NaN 1 1 3 1]);
[~, warnId] = lastwarn;
[~, costU, itU] = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, lbW, [Inf NaN 1 1 3 1]);       % same, but ub(v_tgp) = Inf
rW = costW./costU;
ok3g = mean(convW) >= 0.99 && mean(itW) < 100 && abs(median(rW) - 1) < 1e-3 && strcmp(warnId, 'vUS_3D_quad_fitBatched:unphysicalVtBound') && all(isfinite(xW(:)));
fprintf('3g. ub(v_tgp) = 50e3: converged %.1f%%, mean iterations %.1f (bug: 300), cost vs ub = Inf: median %.6f, warns about units = %d -> %s\n', ...
    100*mean(convW), mean(itW), median(rW), strcmp(warnId, 'vUS_3D_quad_fitBatched:unphysicalVtBound'), passfail(ok3g));

ok3 = ok3a && ok3b && ok3c && ok3d && ok3e && ok3f && ok3g;
if ok1 && ok2 && ok3, fprintf('\nAll checks passed.\n'); else, warning('One or more checks failed -- inspect before use.'); end

function s = passfail(tf)
    if tf, s = 'PASS'; else, s = 'FAIL'; end
end

function tf = throws(fh, idSuffix)
% true if fh() errors with an identifier ending in idSuffix
    try
        fh(); tf = false;
    catch e
        tf = endsWith(e.identifier, idSuffix);
    end
end
