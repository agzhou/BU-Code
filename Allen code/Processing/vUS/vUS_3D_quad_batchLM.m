function [X, cost, iters, conv] = vUS_3D_quad_batchLM(P, X0, LB, UB, opt)
%% Description:
%   Batched box-constrained Levenberg-Marquardt for the corrected
%   Tangelder-profile g1 model: every voxel of a chunk advances in lockstep,
%   so the interpreter overhead of one solver step is paid once per chunk
%   instead of once per voxel. Voxels stay fully independent (own damping,
%   own bounds, own convergence); they only share loop iterations.
%
%   Per iteration, for all active voxels at once:
%     1. JtJ [6x6xB] and J'r [6x1xB] via pagemtimes
%     2. Marquardt-damped step  (JtJ + lambda_b*diag(JtJ)) d = -J'r  for all B
%        systems with one pagemldivide; variables sitting on a bound with the
%        gradient pushing outward are frozen (active-set flag), the step is
%        capped elementwise (crude trust region), and the trial point is
%        clipped to the box
%     3. one batched model + Jacobian evaluation at the trial points
%        (vUS_3D_quad_batchEval.m), accept where the cost decreased
%     4. lambda_b <- lambda_b/3 (accepted) or *4 (rejected), per voxel
%     5. converged voxels are written out and removed from all arrays
%
%   Different from lsqnonlin's trust-region-reflective algorithm, so a
%   minority of voxels (~1% each way on the test data) end in a different
%   local minimum. Class-agnostic: double, single or gpuArray inputs (the
%   gpuArray path is UNTESTED -- no GPU was available when this was written).
%
% Inputs:
%   P: chunk data/constants struct from vUS_3D_quad_batchPack.m
%   X0, LB, UB: [6, B] SCALED start point and bounds (parameters ./ xscale), order
%       [v_tgp, v_zgp, F, DC, k, a]; X0 must lie inside [LB, UB]
%   opt: struct with fields
%       maxIter: iteration cap (300)
%       tolF: converged when an accepted step lowers the cost by <= tolF relative
%       tolX: converged when the largest scaled parameter step <= tolX
%       lam0: initial damping (1; 1e-3 starts more aggressively and lands in worse basins more often)
%       stepCap: max |step| per parameter per iteration in scaled units (0.5; Inf disables)
%
% Outputs:
%   X: [6, B] SCALED solution;  cost: [1, B] final sum of squared residuals
%   iters: [1, B] iterations used;  conv: [1, B] false where maxIter was hit

    B = size(X0, 2);  cls = 'like';
    X  = X0;  LBc = LB;  UBc = UB;  sel = 1:B;  ids = 1:B;
    Xo = X0;  co = zeros(1, B, cls, X0);  iters = zeros(1, B);  conv = false(1, B);
    [r, J, cost] = vUS_3D_quad_batchEval(X, sel, P);
    lam = opt.lam0*ones(1, B, cls, X0);
    I6  = eye(6, cls, X0);
    for it = 1:opt.maxIter
        Bc  = size(X, 2);
        JtJ = pagemtimes(J, 'transpose', J, 'none');                        % 6x6xBc
        g   = reshape(pagemtimes(J, 'transpose', r, 'none'), 6, Bc);        % 6xBc
        dJ  = max(reshape(sum(J.^2, 1), 6, Bc), 1e-12);
        free = ~((X <= LBc) & (g > 0)) & ~((X >= UBc) & (g < 0));           % active-set: freeze variables pushing into a bound
        Fm  = reshape(double(free), 6, 1, Bc);  Fm = cast(Fm, 'like', X0);
        A   = JtJ + I6.*reshape(lam, 1, 1, Bc).*reshape(dJ, 6, 1, Bc) + I6.*reshape(1e-10*max(dJ, [], 1), 1, 1, Bc);   % Marquardt damping + tiny ridge
        A   = A.*(Fm.*permute(Fm, [2 1 3])) + I6.*(1 - Fm);
        d   = -reshape(pagemldivide(A, reshape(g.*cast(free,'like',X0), 6, 1, Bc)), 6, Bc);
        if isfinite(opt.stepCap), d = max(min(d, opt.stepCap), -opt.stepCap); end          % crude trust-region cap (scaled units)
        Xn  = min(max(X + d, LBc), UBc);
        [rn, Jn, cn] = vUS_3D_quad_batchEval(Xn, sel, P);
        acc = cn < cost;
        rel = (cost - cn)./max(cost, realmin('like', X0));
        stp = max(abs(Xn - X), [], 1);
        X(:, acc) = Xn(:, acc);  r(:,:,acc) = rn(:,:,acc);  J(:,:,acc) = Jn(:,:,acc);  cost(acc) = cn(acc);
        lam = min(max(lam.*(4.^double(~acc)).*(1/3).^double(acc), 1e-12), 1e12);
        done = (acc & (rel <= opt.tolF | stp <= opt.tolX)) | (~acc & (stp <= opt.tolX | lam >= 1e10));
        done = gather(done);                                                % host logical: it indexes the host bookkeeping arrays (no-op on CPU)
        if it == opt.maxIter, done(:) = true; end
        if any(done)
            gi = ids(done);  Xo(:, gi) = X(:, done);  co(gi) = cost(done);  iters(gi) = it;  conv(gi) = ~(it == opt.maxIter);
            keep = ~done;
            if ~any(keep), break; end
            X = X(:, keep); LBc = LBc(:, keep); UBc = UBc(:, keep); r = r(:,:,keep); J = J(:,:,keep);
            cost = cost(keep); lam = lam(keep); sel = sel(keep); ids = ids(keep);
        end
    end
    X = Xo;  cost = co;
end
