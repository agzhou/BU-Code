function [X, cost, iters, conv] = batchLM(ev, X0, LB, UB, opt)
%% Description:
%   Model-agnostic batched box-constrained Levenberg-Marquardt: every problem
%   (voxel) of a chunk advances in lockstep, so the interpreter overhead of one
%   solver step is paid once per chunk instead of once per voxel. Problems stay
%   fully independent (own damping, own bounds, own convergence); they only
%   share loop iterations. This is the generic version of
%   vUS_3D_quad_batchLM.m -- same algorithm, but the model enters only through
%   the residual/Jacobian callbacks in ev, so any nP and any model work.
%
%   Per iteration, for all active problems at once:
%     1. JtJ [nPxnPxB] and J'r [nPx1xB] via pagemtimes
%     2. Marquardt-damped step (JtJ + lambda_b*diag(JtJ)) d = -J'r for all B
%        systems with one pagemldivide; variables sitting on a bound with the
%        gradient pushing outward are frozen (active-set flag), the step is
%        capped elementwise (crude trust region), and the trial point is
%        clipped to the box
%     3. one batched residual evaluation at the trial points, accept where the
%        cost decreased
%     4. lambda_b <- lambda_b/3 (accepted) or *4 (rejected), per problem
%     5. the Jacobian is refreshed only for problems whose step was accepted
%        (with a finite-difference Jacobian that saves the nP extra model
%        evaluations for every rejected trial)
%     6. converged problems are written out and removed from all arrays
%
%   Class-agnostic: double, single or gpuArray inputs (the gpuArray path is
%   UNTESTED -- no GPU was available when this was written).
%
% Inputs:
%   ev: struct of function handles (built by fitBatchedLM.m); X is [nP, Bc] SCALED
%       parameters and sel [1, Bc] the indices of the active problems within the chunk
%       ev.analytic: true if the Jacobian comes together with the residual
%       ev.res(X, sel): residual R [Lr, Bc]  (used when ~ev.analytic)
%       ev.jac(X, sel, R): Jacobian J [Lr, nP, Bc] at X, given the residual R there (used when ~ev.analytic)
%       ev.resJac(X, sel): [R, J] together (used when ev.analytic)
%   X0, LB, UB: [nP, B] SCALED start point and bounds; X0 must lie inside [LB, UB]
%   opt: struct with fields
%       maxIter: iteration cap
%       tolF: converged when an accepted step lowers the cost by <= tolF relative
%       tolX: converged when the largest scaled parameter step <= tolX
%       lam0: initial damping (1; 1e-3 starts more aggressively and lands in worse basins more often)
%       stepCap: max |step| per parameter per iteration in scaled units, [nP, 1]; Inf disables the cap
%
% Outputs:
%   X: [nP, B] SCALED solution;  cost: [1, B] final sum of squared residuals
%       (NaN for problems whose residual was not finite at the start point; X0 is returned for those)
%   iters: [1, B] iterations used;  conv: [1, B] false where maxIter was hit or the start was not finite

    nP = size(X0, 1);  B = size(X0, 2);
    X  = X0;  LBc = LB;  UBc = UB;  sel = 1:B;
    Xo = X0;  co = nan(1, B, 'like', X0);  iters = zeros(1, B);  conv = false(1, B);
    cap = cast(opt.stepCap(:), 'like', X0);
    I   = eye(nP, 'like', X0);

    [R, J] = residJac(ev, X, sel);
    cost = sum(R.^2, 1);
    bad = gather(~isfinite(cost));                                          % non-finite model at the start point: give up on those problems
    if any(bad)
        keep = ~bad;
        if ~any(keep), X = Xo; cost = co; return, end
        X = X(:, keep); LBc = LBc(:, keep); UBc = UBc(:, keep); R = R(:, keep); J = J(:, :, keep);
        cost = cost(keep); sel = sel(keep);
    end
    lam = opt.lam0*ones(1, numel(sel), 'like', X0);

    for it = 1:opt.maxIter
        Bc  = size(X, 2);  Lr = size(R, 1);
        JtJ = pagemtimes(J, 'transpose', J, 'none');                        % nP x nP x Bc
        g   = reshape(pagemtimes(J, 'transpose', reshape(R, Lr, 1, Bc), 'none'), nP, Bc);   % nP x Bc
        dJ  = max(reshape(sum(J.^2, 1), nP, Bc), 1e-12);
        free = ~((X <= LBc) & (g > 0)) & ~((X >= UBc) & (g < 0));           % active-set: freeze variables pushing into a bound
        Fm  = reshape(cast(free, 'like', X0), nP, 1, Bc);
        A   = JtJ + I.*reshape(lam, 1, 1, Bc).*reshape(dJ, nP, 1, Bc) + I.*reshape(1e-10*max(dJ, [], 1), 1, 1, Bc);   % Marquardt damping + tiny ridge
        A   = A.*(Fm.*permute(Fm, [2 1 3])) + I.*(1 - Fm);
        d   = -reshape(pagemldivide(A, reshape(g.*cast(free, 'like', X0), nP, 1, Bc)), nP, Bc);
        d   = max(min(d, cap), -cap);                                       % crude trust-region cap (scaled units); Inf = no cap
        Xn  = min(max(X + d, LBc), UBc);

        if ev.analytic, [Rn, Jn] = ev.resJac(Xn, sel); else, Rn = ev.res(Xn, sel); end
        cn  = sum(Rn.^2, 1);
        acc = cn < cost;                                                    % a non-finite trial cost is never accepted
        rel = (cost - cn)./max(cost, realmin('like', X0));
        stp = max(abs(Xn - X), [], 1);
        X(:, acc) = Xn(:, acc);  R(:, acc) = Rn(:, acc);  cost(acc) = cn(acc);
        lam = min(max(lam.*(4.^double(~acc)).*(1/3).^double(acc), 1e-12), 1e12);
        doneNat = (acc & (rel <= opt.tolF | stp <= opt.tolX)) | (~acc & (stp <= opt.tolX | lam >= 1e10));
        doneNat = gather(doneNat);                                          % host logicals: they index the host bookkeeping arrays (no-op on CPU)
        acc = gather(acc);
        done = doneNat;
        if it == opt.maxIter, done(:) = true; end

        upd = acc & ~done;                                                  % Jacobian at the new point, only where it will be used
        if any(upd)
            if ev.analytic
                J(:, :, upd) = Jn(:, :, upd);
            else
                J(:, :, upd) = ev.jac(X(:, upd), sel(upd), R(:, upd));
            end
        end

        if any(done)
            gi = sel(done);  Xo(:, gi) = X(:, done);  co(gi) = cost(done);  iters(gi) = it;  conv(gi) = doneNat(done);
            keep = ~done;
            if ~any(keep), break; end
            X = X(:, keep); LBc = LBc(:, keep); UBc = UBc(:, keep); R = R(:, keep); J = J(:, :, keep);
            cost = cost(keep); lam = lam(keep); sel = sel(keep);
        end
    end
    X = Xo;  cost = co;
end

function [R, J] = residJac(ev, X, sel)
% Residual and Jacobian together (one call if the model supplies its own Jacobian, else res + jac)
    if ev.analytic
        [R, J] = ev.resJac(X, sel);
    else
        R = ev.res(X, sel);
        J = ev.jac(X, sel, R);
    end
end
