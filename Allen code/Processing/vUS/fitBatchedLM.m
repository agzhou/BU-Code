function [x, cost, iters, converged] = fitBatchedLM(model, x0, data, lb, ub, opts)
%% Description:
%   Fit an ARBITRARY model to MANY voxels at once with the batched
%   Levenberg-Marquardt solver (batchLM.m). This is the model-agnostic version
%   of vUS_3D_quad_fitBatched.m: the model, start point, bounds, fit window and
%   weights are all inputs, so a new model needs no new solver code. Same
%   algorithm, chunking, scaling and bound conventions as vUS_3D_quad_fitBatched.
%
%   The model is a function handle that is evaluated for ALL voxels of a chunk in
%   one call (that is where the speed comes from):
%
%       Y = model(X, cols, ids)            or   [Y, J] = model(X, cols, ids)   (opts.jacobian = 'analytic')
%
%       X:    [nP, Bc] parameters in PHYSICAL units, one column per voxel
%       cols: [L, 1] indices of the data samples (columns of data) to evaluate, e.g. tau(cols)
%       ids:  [1, Bc] index of each column's voxel (row of data / x0 / lb / ub), for per-voxel constants
%       Y:    [L, Bc] model prediction (or [L, 1, Bc]); complex or real, column b belongs to voxel ids(b)
%       J:    [L, nP, Bc] analytic dY/dX in physical units (only for opts.jacobian = 'analytic')
%
%   Handles with fewer inputs are fine: @(X) ..., @(X, cols) ... or @(X, cols, ids) ....
%   The residual is (model - data), stacked [real; imag] when data is complex, so the
%   solver minimises the same sum of squares as lsqnonlin on a split complex residual.
%   Nothing needs to be differentiated by hand: without a Jacobian the solver builds one
%   from batched finite differences (nP extra model calls per iteration, each over the
%   whole chunk, and only for the voxels whose last step was accepted).
%
%   Example -- the corrected Tangelder-profile g1 model (replaces vUS_3D_quad_fitBatched; for the 2D
%   [v_xgp, v_zgp] model of vUS_2D_newmodel.m use vUS_2D_quad_batchModel with sigma = [sigma_x, sigma_z]):
%       [s, w] = gaussLegendre01(48);
%       model  = @(X, cols) vUS_3D_quad_batchModel(X, tau(cols), PP.k0, sigma, s, w);     % [Y, J] with analytic Jacobian
%       vi     = find(maskToUse);
%       nv     = numel(vi);
%       x0     = [5e-3*ones(nv,1), Vz0(vi), ones(nv,1), zeros(nv,1), 2.5*ones(nv,1), ones(nv,1)];   % or one [1,6] row for all voxels
%       lb     = [0, -50e-3, 0, 0, 2, 0];   ub = [250e-3, 50e-3, 1, 1, 3, 1];                        % [1,6] or [nVox,6]; Inf/-Inf allowed
%       opts   = struct('window', tau_decayed_ind(vi), 't1i', 2, 'jacobian', 'analytic', 'xscale', [1e-2 1e-2 1 1 1 1]);
%       [x, cost] = fitBatchedLM(model, x0, g1_exp{j}(vi, :), lb, ub, opts);
%   With opts.jacobian = 'forward' (the default) and a model that returns only Y, e.g.
%       model  = @(X, cols) vUS_3D_quad_batchModel(X, tau(cols), PP.k0, sigma, s, w);   % called with nargout = 1 -> only Y
%   the same fit needs no analytic Jacobian. A model written for ONE voxel at a time (e.g. the
%   integral()-based vUS_3D_num_wrapper.m) can be used through batchModelFromSingle.m; nothing
%   is vectorised then, so it is only as fast as the per-voxel model.
%   Fixing a parameter: set lb = ub for it.
%
%   Measured with the vUS model on 2,000 synthetic voxels, 8 threads (test_fitBatchedLM.m section 7; range over
%   three runs):
%       vUS_3D_quad_fitBatched 1.05-1.09 s | fitBatchedLM analytic 1.03-1.06 s, forward 3.2-3.4 s, central 6.5-6.9 s |
%       serial lsqnonlin loop with an analytic Jacobian 32-39 s.
%   So the analytic path costs nothing over the specialised solver, and finite differences (no derivation
%   needed) cost about 3x (forward) / 6x (central) that but are still ~10x / ~5x faster than the serial loop.
%   With a Jacobian-free model the fits reached the same cost as the analytic ones (no voxel >1% worse than
%   vUS_3D_quad_fitBatched in 200; 1 of 200 better with forward, 5 with central differences).
%
% Inputs:
%   model: function handle, see above
%   x0: [1, nP] (same start point for every voxel) or [nVox, nP] initial guess in physical units; nP is
%       size(x0, 2). Clipped into [lb, ub] (nudged 0.1% into the box, never by more than 1e-3 scaled units)
%   data: [nVox, nT] observations (complex or real), one row per voxel
%   lb, ub: (optional) lower/upper bounds on the nP parameters in physical units, each [] (unbounded),
%       [1, nP] (same for every voxel) or [nVox, nP] (per voxel). Inf/-Inf allowed. NaN in any entry means
%       "unbounded" for that entry, so partial specifications are easy. lb must be <= ub for every entry.
%   opts: (optional) struct; any field may be omitted
%       window: which samples of each voxel enter the fit. [] (default) = all nT samples of every voxel;
%           [nVox, 1] numeric = index of the LAST sample per voxel, the fit window is opts.t1i:window
%           (findTauDecayed's tau_decayed_ind); [nVox, nT] logical = arbitrary per-voxel mask.
%           Voxels with an empty window are not fitted (x0 is returned for them, with a warning)
%       t1i: first sample of the window when window is a vector of last indices (1; use 2 to skip tau = 0)
%       weights: [] (default), [1, nT] / [nT, 1] (same per voxel) or [nVox, nT], non-negative. Multiplies the
%           RESIDUALS (not the squared residuals) sample by sample, as OF_weight in vUS_3D_num_OF.m does
%       jacobian: 'forward' (default) forward finite differences, 'central' central differences (2x the
%           model calls, more accurate near the optimum), or 'analytic' (the model's second output).
%           Finite-difference steps stay inside the bounds, and parameters fixed by lb = ub for every
%           voxel of a chunk are not perturbed at all
%       fdStep: relative finite-difference step in scaled units (sqrt(eps) forward, eps^(1/3) central)
%       xscale: [1, nP] parameter scale; parameters are fitted as x./xscale so they are all O(1). Default
%           (empty = auto): median |x0| per parameter, else the median finite box width, else 1. The step cap
%           and the step tolerance act in these scaled units, so scaling matters: on 600 synthetic vUS voxels with
%           lb/ub = [0 -50e-3 0 0 2 0]/[250e-3 50e-3 1 1 3 1], xscale = ones left 3% of voxels unconverged and 51
%           voxels >1% worse than the tuned [1e-2 1e-2 1 1 1 1] (what vUS_3D_quad_fitBatched used); the auto
%           scale ([5e-3 1.3e-2 1 1 2.5 1] there) converged 100% with 1 such voxel
%       chunkSize: voxels per batch (1024; the best size depends on the model's memory footprint)
%       precision: 'double' (default) or 'single' (the model receives single-precision parameters; in the vUS
%           test 2 of 200 voxels ended >1% worse than in double)
%       tol: convergence tolerance (1e-10 relative cost change; scaled step tolerance is 10*tol);
%           default 1e-5 for 'single'
%       maxIter: iteration cap per voxel (300)
%       lam0: initial Levenberg-Marquardt damping (1)
%       stepCap: max scaled step per parameter per iteration; scalar or [1, nP], Inf allowed (0.5). Raise it for
%           a parameter that is effectively unbounded and runs away along a flat direction
%       useGPU: run the chunks on a GPU via gpuArray (false). UNTESTED -- no GPU was available when written;
%           the model must accept gpuArray parameters
%
% Outputs:
%   x: [nVox, nP] fitted parameters in physical units
%   cost: [nVox, 1] final sum of squared (weighted) residuals over the fit window (lsqnonlin's resnorm);
%       NaN where the model was not finite at the start point
%   iters: [nVox, 1] iterations used
%   converged: [nVox, 1] false where the iteration cap was hit, the window was empty or the start was not finite

    if nargin < 4, lb = []; end
    if nargin < 5, ub = []; end
    if nargin < 6, opts = []; end
    if ~isa(model, 'function_handle')
        error('fitBatchedLM:badModel', 'model must be a function handle: Y = model(X, cols, ids).')
    end
    if isstruct(lb) || isstruct(ub)
        error('fitBatchedLM:optsMoved', 'Call as fitBatchedLM(model, x0, data, lb, ub, opts); use [] for lb/ub to leave them unbounded.')
    end

    o = struct('window', [], 't1i', 1, 'weights', [], 'jacobian', 'forward', 'fdStep', [], 'xscale', [], ...
        'chunkSize', 1024, 'precision', 'double', 'tol', [], 'maxIter', 300, 'lam0', 1, 'stepCap', 0.5, 'useGPU', false);
    if ~isempty(opts)
        fn = fieldnames(opts);
        for i = 1:numel(fn)
            if ~isfield(o, fn{i}), error('fitBatchedLM:badOption', 'Unknown option "%s".', fn{i}); end
            o.(fn{i}) = opts.(fn{i});
        end
    end
    mode = validatestring(o.jacobian, {'forward', 'central', 'analytic'});
    cls  = validatestring(o.precision, {'double', 'single'});
    if isempty(o.tol), o.tol = 1e-10; if strcmp(cls, 'single'), o.tol = 1e-5; end, end
    if isempty(o.fdStep)
        if strcmp(mode, 'central'), o.fdStep = eps(cls)^(1/3); else, o.fdStep = sqrt(eps(cls)); end
    end

    if ~isnumeric(data) || ~ismatrix(data), error('fitBatchedLM:badData', 'data must be a numeric [nVox x nT] matrix.'), end
    [n, nT] = size(data);
    nP = size(x0, 2);
    if ~isnumeric(x0) || ~ismatrix(x0) || nP < 1 || ~(size(x0, 1) == 1 || size(x0, 1) == n) || any(~isfinite(x0(:)))
        error('fitBatchedLM:badX0', 'x0 must be a finite [1 x nP] or [nVox x nP] matrix (nVox = %d), got size %s.', n, mat2str(size(x0)))
    end
    x0M = double(x0);  if size(x0M, 1) == 1, x0M = repmat(x0M, n, 1); end

    % Bounds per voxel [n, nP]; NaN entries mean unbounded
    lbM = resolveBounds(lb, n, nP, 'lb', -Inf);
    ubM = resolveBounds(ub, n, nP, 'ub', Inf);
    if any(lbM(:) > ubM(:))
        [bv, bp] = find(lbM > ubM, 1);
        error('fitBatchedLM:badBounds', 'lb > ub for voxel %d, parameter %d (lb = %g, ub = %g).', bv, bp, lbM(bv,bp), ubM(bv,bp))
    end
    if any(lbM(:) == Inf) || any(ubM(:) == -Inf), error('fitBatchedLM:badBounds', 'lb must be < Inf and ub > -Inf.'), end

    % Parameter scale and step cap
    if isempty(o.xscale)
        xscale = ones(1, nP);
        for j = 1:nP
            a0 = abs(x0M(:, j));  a0 = a0(a0 > 0);
            wd = ubM(:, j) - lbM(:, j);  wd = wd(isfinite(wd) & wd > 0);
            if ~isempty(a0), xscale(j) = median(a0); elseif ~isempty(wd), xscale(j) = median(wd); end
        end
    else
        xscale = double(o.xscale(:)).';
        if numel(xscale) ~= nP || any(~isfinite(xscale)) || any(xscale <= 0)
            error('fitBatchedLM:badXscale', 'opts.xscale must be a positive finite [1 x %d] vector.', nP)
        end
    end
    if ~any(numel(o.stepCap) == [1 nP]) || any(o.stepCap(:) <= 0) || any(isnan(o.stepCap(:)))
        error('fitBatchedLM:badStepCap', 'opts.stepCap must be a positive scalar or a positive [1 x %d] vector (Inf allowed).', nP)
    end
    stepCap = double(o.stepCap(:)) .* ones(nP, 1);
    optLM = struct('maxIter', o.maxIter, 'tolF', o.tol, 'tolX', 10*o.tol, 'lam0', o.lam0, 'stepCap', stepCap);

    % Fit window per voxel: first/last sample (and an explicit mask, if given)
    [first, last, mask] = parseWindow(o.window, o.t1i, n, nT);
    empty = last < first;
    if any(empty)
        warning('fitBatchedLM:emptyWindow', '%d voxel(s) have an empty fit window and are not fitted (x0 returned).', nnz(empty))
    end
    wts = o.weights;  wtsVec = false;
    if ~isempty(wts)
        wtsVec = isvector(wts) && numel(wts) == nT;            % one weight per sample, shared by every voxel
        if ~isnumeric(wts) || any(~isfinite(wts(:))) || any(wts(:) < 0) || ~(wtsVec || isequal(size(wts), [n nT]))
            error('fitBatchedLM:badWeights', 'weights must be finite, non-negative and [1 x nT], [nT x 1] or [nVox x nT] (nT = %d).', nT)
        end
        if wtsVec, wts = reshape(wts, [], 1); end
    end
    cplx = ~isreal(data);

    nIn = nargin(model);                                       % let the model take fewer than 3 inputs
    if nIn == 1, mdl = @(X, cols, ids) model(X);
    elseif nIn == 2, mdl = @(X, cols, ids) model(X, cols);
    elseif nIn == 0, error('fitBatchedLM:badModel', 'model must take at least the parameters X.')
    else, mdl = model;
    end

    x = x0M;  cost = nan(n, 1);  iters = zeros(n, 1);  converged = false(n, 1);
    valid = find(~empty);
    [~, ord] = sortrows([last(valid), first(valid)]);          % length-sorted chunks keep zero-padding small
    order = valid(ord);
    xs = xscale(:);
    for a0 = 1:o.chunkSize:numel(order)
        c = order(a0:min(a0 + o.chunkSize - 1, numel(order)));
        cols = (min(first(c)):max(last(c))).';

        Wc = (cols >= first(c).') & (cols <= last(c).');       % [L, B] fit window (zero-padded outside)
        if ~isempty(mask), Wc = Wc & mask(c, cols).'; end
        Wc = double(Wc);
        if ~isempty(wts)
            if wtsVec, Wc = Wc .* wts(cols); else, Wc = Wc .* wts(c, cols).'; end
        end
        D = data(c, cols).';                                   % [L, B]; plain transpose: complex data must not be conjugated
        Z = (Wc == 0);  D(Z) = 0;                              % samples outside the window: no data, no residual (and no NaN leaking in)
        if cplx, Z = repmat(Z, 2, 1); end                      % [real; imag] residual rows

        LB = cast(lbM(c, :).' ./ xs, cls);  UB = cast(ubM(c, :).' ./ xs, cls);  X0 = cast(x0M(c, :).' ./ xs, cls);
        d = min(1e-3*(UB - LB), 1e-3);                         % start strictly inside the box: nudge = 0.1% of the box width, but never more than 1e-3
        X0 = min(max(X0, LB + d), UB - d);                     % (scaled units); an uncapped nudge would push starts to 0.1% of a huge ub

        Q = struct('model', mdl, 'cols', cols, 'ids', c.', 'D', cast(D, cls), 'W', cast(Wc, cls), 'Z', Z, 'hasZero', any(Z(:)), ...
            'cplx', cplx, 'xscale', cast(xs, cls), 'LB', LB, 'UB', UB, 'central', strcmp(mode, 'central'), ...
            'fdStep', cast(o.fdStep, cls), 'active', find(any(UB > LB, 2)).');
        if o.useGPU
            for f = {'D', 'W', 'Z', 'xscale', 'LB', 'UB'}, Q.(f{1}) = gpuArray(Q.(f{1})); end
            X0 = gpuArray(X0);  LB = gpuArray(LB);  UB = gpuArray(UB);
        end
        ev = struct('analytic', strcmp(mode, 'analytic'), ...
            'res',    @(Xs, sel) evalRes(Q, Xs, sel), ...
            'jac',    @(Xs, sel, R0) evalJacFD(Q, Xs, sel, R0), ...
            'resJac', @(Xs, sel) evalResJac(Q, Xs, sel));

        [Xc, cc, ic, cv] = batchLM(ev, X0, LB, UB, optLM);

        x(c, :) = double(gather(Xc .* cast(xs, 'like', Xc))).';
        cost(c) = double(gather(cc)).';  iters(c) = double(gather(ic)).';  converged(c) = logical(gather(cv)).';
    end
    nBad = nnz(isnan(cost(valid)));
    if nBad > 0
        warning('fitBatchedLM:badStart', '%d voxel(s) gave a non-finite model residual at the start point and were not fitted (x0 returned, cost = NaN).', nBad)
    end
end

function M = resolveBounds(u, n, nP, name, def)
% Expand a user bound input ([], [1,nP] or [n,nP]; NaN = unbounded) to an [n,nP] matrix
    if isempty(u), M = repmat(def, n, nP); return, end
    if ~isnumeric(u) || ~ismatrix(u) || size(u, 2) ~= nP || ~(size(u, 1) == 1 || size(u, 1) == n)
        error('fitBatchedLM:badBounds', '%s must be [], [1 x %d] or [nVox x %d] (nVox = %d), got size %s.', name, nP, nP, n, mat2str(size(u)))
    end
    M = double(u);
    if size(M, 1) == 1, M = repmat(M, n, 1); end
    M(isnan(M)) = def;
end

function [first, last, mask] = parseWindow(win, t1i, n, nT)
% First/last sample index per voxel (and the explicit mask, if the window was given as one)
    mask = [];
    if islogical(win) && ~isempty(win)
        if ~isequal(size(win), [n nT]), error('fitBatchedLM:badWindow', 'A logical window must be [nVox x nT] = [%d x %d].', n, nT), end
        mask = win;
        [hasAny, first] = max(win, [], 2);                     % first true
        [~, lr] = max(fliplr(win), [], 2);  last = nT - lr + 1;
        first(~hasAny) = 1;  last(~hasAny) = 0;                % empty windows
        return
    end
    if ~(isscalar(t1i) && t1i >= 1 && t1i <= nT && t1i == floor(t1i))
        error('fitBatchedLM:badWindow', 'opts.t1i must be an integer between 1 and nT = %d.', nT)
    end
    first = t1i*ones(n, 1);
    if isempty(win)
        last = nT*ones(n, 1);
    else
        if ~isnumeric(win) || numel(win) ~= n || any(~isfinite(win(:)))
            error('fitBatchedLM:badWindow', 'A numeric window must hold one finite last-sample index per voxel (nVox = %d).', n)
        end
        last = min(floor(double(win(:))), nT);
    end
end

function Y = callModel(Q, Xs, sel)
% Model prediction [L, Bc] at scaled parameters Xs for the active voxels sel
    Y = checkY(Q, Q.model(Xs .* Q.xscale, Q.cols, Q.ids(sel)), numel(sel));
end

function Y = checkY(Q, Y, Bc)
    L = numel(Q.cols);  sz = size(Y);
    if numel(sz) == 3 && sz(2) == 1, Y = reshape(Y, sz(1), sz(3)); sz = size(Y); end   % [L,1,Bc] -> [L,Bc]
    if ~isequal(sz, [L Bc])
        error('fitBatchedLM:badModelOutput', 'The model must return Y as [L x Bc] = [%d x %d] (one column per voxel), got size %s.', L, Bc, mat2str(size(Y)))
    end
end

function R = toResidual(Q, Y, sel)
% Weighted, masked residual (model - data): [Lr, Bc], with Lr = 2L for complex data ([real; imag])
    if ~Q.cplx && ~isreal(Y)
        if any(imag(Y(:)) ~= 0), error('fitBatchedLM:complexModel', 'The model returned complex values but data is real; fit complex data or return a real model.'), end
        Y = real(Y);
    end
    W = Q.W(:, sel);  E = Y - Q.D(:, sel);
    if Q.cplx, R = [real(E); imag(E)] .* [W; W]; else, R = E .* W; end
    if Q.hasZero, R(Q.Z(:, sel)) = 0; end
end

function R = evalRes(Q, Xs, sel)
    R = toResidual(Q, callModel(Q, Xs, sel), sel);
end

function [R, J] = evalResJac(Q, Xs, sel)
% Residual and the model's own analytic Jacobian (returned w.r.t. the SCALED parameters)
    Bc = numel(sel);  nP = size(Xs, 1);  L = numel(Q.cols);
    try
        [Y, Jm] = Q.model(Xs .* Q.xscale, Q.cols, Q.ids(sel));
    catch e
        if any(strcmp(e.identifier, {'MATLAB:maxlhs', 'MATLAB:TooManyOutputs'}))
            error('fitBatchedLM:noJacobian', 'opts.jacobian = ''analytic'' needs a model with two outputs, [Y, J] = model(X, cols, ids); use ''forward'' or ''central'' otherwise.')
        end
        rethrow(e)
    end
    R = toResidual(Q, checkY(Q, Y, Bc), sel);
    if ~isequal(size(Jm, 1:3), [L nP Bc])                      % size(., 1:3): a one-voxel batch has no trailing singleton
        error('fitBatchedLM:badModelOutput', 'The analytic Jacobian must be [L x nP x Bc] = [%d x %d x %d], got size %s.', L, nP, Bc, mat2str(size(Jm)))
    end
    Jm = Jm .* reshape(Q.xscale, 1, nP);                       % chain rule for the parameter scaling
    W = Q.W(:, sel);
    if Q.cplx, J = [real(Jm); imag(Jm)]; W = [W; W]; else, J = real(Jm); end
    J = J .* reshape(W, [], 1, Bc);
    J(~isfinite(J)) = 0;
end

function J = evalJacFD(Q, Xs, sel, R0)
% Batched finite-difference Jacobian of the residual w.r.t. the SCALED parameters, [Lr, nP, Bc]. R0 is the residual at Xs.
% One model call over all Bc voxels per perturbed parameter (two for central differences). Steps stay inside the box.
    [nP, Bc] = size(Xs);  Lr = size(R0, 1);
    J  = zeros(Lr, nP, Bc, 'like', R0);
    LB = Q.LB(:, sel);  UB = Q.UB(:, sel);
    for j = Q.active                                           % parameters fixed for every voxel of the chunk (lb = ub) are skipped
        xj = Xs(j, :);
        h  = Q.fdStep .* max(1, abs(xj));                      % [1, Bc], scaled units
        up = UB(j, :) - xj;  dn = xj - LB(j, :);               % room to each bound
        if Q.central
            hp = min(h, up);  hm = min(h, dn);                 % symmetric where possible, one-sided at a bound
            fx = (hp + hm) <= 0;  hp(fx) = h(fx);              % fixed in this voxel: its column is frozen by the box anyway; any step avoids 0/0
            xp = xj + hp;  xm = xj - hm;  hp = xp - xj;  hm = xj - xm;   % the steps actually taken (exactly representable)
            Rp = perturbed(Q, Xs, sel, j, xp);  Rm = perturbed(Q, Xs, sel, j, xm);
            dR = (Rp - Rm) ./ (hp + hm);
        else
            sg = ones(1, Bc, 'like', xj);  sg(up < h & dn > up) = -1;   % step down only where up does not fit and down has more room
            h  = min(h, max(up, dn));                          % shrink into a box narrower than the step
            h(h <= 0) = Q.fdStep;                              % fixed in this voxel (see above)
            xp = xj + sg.*h;  hs = xp - xj;
            dR = (perturbed(Q, Xs, sel, j, xp) - R0) ./ hs;
        end
        J(:, j, :) = reshape(dR, Lr, 1, Bc);
    end
    J(~isfinite(J)) = 0;                                       % a model that fails at a perturbed point gives "no information" for that column
end

function R = perturbed(Q, Xs, sel, j, xj)
% Residual with parameter j replaced by xj (a [1, Bc] row of scaled values)
    Xs(j, :) = xj;
    R = toResidual(Q, callModel(Q, Xs, sel), sel);
end
