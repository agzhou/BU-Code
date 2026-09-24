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

if ok1 && ok2, fprintf('\nAll checks passed.\n'); else, warning('One or more checks failed -- inspect before use.'); end

function s = passfail(tf)
    if tf, s = 'PASS'; else, s = 'FAIL'; end
end
