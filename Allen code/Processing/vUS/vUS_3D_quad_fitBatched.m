function [x, cost, iters, converged] = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, lb, ub, opts)
%% Description:
%   Fit the corrected Tangelder-profile g1 model to MANY voxels at once with
%   the batched Levenberg-Marquardt solver (vUS_3D_quad_batchLM.m). Drop-in
%   alternative to the per-voxel lsqnonlin loop in vUS_3D_newmodel.m section 5:
%   same model, same parameters, same default start point and bounds, but all
%   voxels advance together. On the 5,206 masked voxels of a real 3D dataset
%   this was ~12x faster than serial lsqnonlin (8 threads; ~4.6x on 1 thread).
%
%   Example (replaces the "for vi = 1:num_voxels ... lsqnonlin ..." loop):
%       vi = find(maskToUse);                     % voxels to fit
%       [x, cost] = vUS_3D_quad_fitBatched(g1_exp{j}(vi, :), tau_decayed_ind(vi), Vz0(vi), tau, PP.k0, sigma);
%       v_tgp_stacked(vi) = x(:,1);  v_zgp_stacked(vi) = x(:,2);  F_stacked(vi) = x(:,3);
%       DC_stacked(vi)    = x(:,4);  k_stacked(vi)     = x(:,5);  a_stacked(vi) = x(:,6);
%
%   With explicit bounds, same convention as lsqnonlin (Inf/-Inf allowed):
%       lb = [0, -Inf, 0, 0, 2, 0];   ub = [Inf, Inf, 1, 1, 3, 1];
%       [x, cost] = vUS_3D_quad_fitBatched(g1_exp{j}(vi, :), tau_decayed_ind(vi), Vz0(vi), tau, PP.k0, sigma, lb, ub);
%   Per-voxel bounds ([nVox, 6] matrices) and partial overrides work too:
%       lb = nan(nnz(vi), 6);  lb(:,5) = 2.2;   % only bound k from below; everything else stays at its default
%       [x, cost] = vUS_3D_quad_fitBatched(..., sigma, lb, []);
%   Fixing a parameter: set lb = ub for it (e.g. lb(5) = ub(5) = 2 fits the k = 2 profile).
%
% Inputs:
%   g1: [nVox, nTau] complex g1 of the voxels to fit, tau starting at 0
%   tdi: [nVox, 1] index of the last lag to fit per voxel (findTauDecayed's tau_decayed_ind);
%       the fit window is opts.t1i:tdi
%   Vz0: [nVox, 1] initial axial group velocity [m/s] (findVzPhaseDiff), used for the default start
%       point and the default v_zgp bounds
%   tau: [nTau, 1] time lags [s];  k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_y, sigma_z] [m]; sigma_x must equal sigma_y
%   lb, ub: (optional) lower/upper bounds on [v_tgp, v_zgp, F, DC, k, a] in physical units, each
%       [] (all defaults), [1, 6] (same bounds for every voxel) or [nVox, 6] (per voxel).
%       Inf/-Inf are allowed. NaN in any entry means "use the default for that entry":
%           default lb = [0,          Vz0 - opts.vzHalfWidth, 0, 0, 2, 0]
%           default ub = [sqrt(2)*30e-3, Vz0 + opts.vzHalfWidth, 1, 1, 3, 1]
%       lb must be <= ub for every entry. Velocities are in m/s (50 mm/s = 50e-3). A finite ub(v_tgp) above 1 m/s
%       (e.g. 50e3, a common unit slip for 50e-3) triggers the warning 'vUS_3D_quad_fitBatched:unphysicalVtBound',
%       because v_tgp is then effectively unbounded and its fitted values are unreliable.
%   opts: (optional) struct; any field may be omitted
%       t1i: first lag index to fit (2)
%       nNodes: Gauss-Legendre nodes (48; 24 was accurate to 4e-9 on real fit windows)
%       chunkSize: voxels per batch (1024 was fastest on CPU)
%       precision: 'double' (default) or 'single' (~2x faster on CPU; slightly more voxels in worse minima)
%       tol: convergence tolerance (1e-10 relative cost change; scaled step tolerance is 10*tol; use ~1e-5 for single)
%       maxIter: iteration cap per voxel (300)
%       lam0: initial Levenberg-Marquardt damping (1)
%       stepCap: max scaled step per parameter per iteration; scalar or [1,6]. Default (NaN = auto): 0.5, except
%           [10 0.5 0.5 0.5 0.5 0.5] when the v_tgp upper bound is above 1 m/s or infinite (effectively unbounded).
%           An unbounded v_tgp runs away along a flat direction and a 0.5 cap makes it crawl. On real data with
%           lb/ub = [0 -Inf 0 0 2 0]/[Inf Inf 1 1 3 1] this took convergence from 92% to 99.8%, run time 2.4x lower, and
%           voxels ending >1% worse than serial lsqnonlin from 689 to ~250 (of 5,206). Larger v_tgp caps (30, 100) were
%           worse. NOTE: with v_tgp unbounded the fit is degenerate -- serial lsqnonlin also returned v_tgp > 100 mm/s
%           in 85% of those voxels.
%       xscale: [1,6] parameter scale (1e-2 1e-2 1 1 1 1: v_tgp, v_zgp in units of 10 mm/s)
%       x0: [1,6] start point, NaN in position 2 = use Vz0 ([5e-3 NaN 1 0 2.5 1]); clipped into [lb, ub]
%       vzHalfWidth: half-width of the DEFAULT v_zgp bound around Vz0 [m/s] (0.01)
%       useGPU: run the chunks on a GPU via gpuArray (false). UNTESTED -- no GPU was available when written
%
% Outputs:
%   x: [nVox, 6] fitted [v_tgp, v_zgp, F, DC, k, a] in physical units (m/s for the velocities)
%   cost: [nVox, 1] final sum of squared residuals over the fit window (lsqnonlin's resnorm)
%   iters: [nVox, 1] iterations used
%   converged: [nVox, 1] false where the iteration cap was hit

    if nargin < 7, lb = []; end
    if nargin < 8, ub = []; end
    if nargin < 9, opts = []; end
    if isstruct(lb)
        error('vUS_3D_quad_fitBatched:optsMoved', ...
            'The 7th input is now lb (then ub, then opts): call as (g1, tdi, Vz0, tau, k0, sigma, lb, ub, opts). Use [] for lb/ub to keep the defaults.')
    end

    o = struct('t1i', 2, 'nNodes', 48, 'chunkSize', 1024, 'precision', 'double', 'tol', 1e-10, 'maxIter', 300, ...
        'lam0', 1, 'stepCap', NaN, 'xscale', [1e-2 1e-2 1 1 1 1], 'x0', [5e-3 NaN 1 0 2.5 1], ...
        'vzHalfWidth', 0.01, 'useGPU', false);
    if ~isempty(opts)
        fn = fieldnames(opts);
        for i = 1:numel(fn)
            if any(strcmp(fn{i}, {'lb', 'ub'}))
                error('vUS_3D_quad_fitBatched:optsMoved', 'lb and ub are now inputs 7 and 8, not fields of opts.')
            end
            if ~isfield(o, fn{i}), error('vUS_3D_quad_fitBatched:badOption', 'Unknown option "%s".', fn{i}); end
            o.(fn{i}) = opts.(fn{i});
        end
    end
    if sigma(1) ~= sigma(2)
        error('vUS_3D_quad_fitBatched:sigmaMismatch', 'For the transverse-combined model, sigma_x and sigma_y must be equal.')
    end

    n = numel(tdi);  tdi = tdi(:);  Vz0 = double(Vz0(:));  tau = double(tau(:));
    if size(g1, 1) ~= n || numel(Vz0) ~= n, error('vUS_3D_quad_fitBatched:sizeMismatch', 'g1, tdi and Vz0 must have one row/element per voxel.'), end
    [s, w] = gaussLegendre01(o.nNodes);
    x0row = o.x0;

    % Bounds per voxel [n, 6]: defaults, overridden by whatever lb/ub specify (NaN = keep the default)
    lbDef = repmat([0 NaN 0 0 2 0], n, 1);                lbDef(:,2) = Vz0 - o.vzHalfWidth;
    ubDef = repmat([sqrt(2)*30e-3 NaN 1 1 3 1], n, 1);    ubDef(:,2) = Vz0 + o.vzHalfWidth;
    lbM = resolveBounds(lb, lbDef, n, 'lb');
    ubM = resolveBounds(ub, ubDef, n, 'ub');
    if any(lbM(:) > ubM(:))
        [bv, bp] = find(lbM > ubM, 1);
        error('vUS_3D_quad_fitBatched:badBounds', 'lb > ub for voxel %d, parameter %d (lb = %g, ub = %g).', bv, bp, lbM(bv,bp), ubM(bv,bp))
    end

    if any(isfinite(ubM(:,1)) & ubM(:,1) > 1)                 % velocities are in m/s: a bound above 1 m/s is almost certainly a unit slip (50 mm/s = 50e-3)
        warning('vUS_3D_quad_fitBatched:unphysicalVtBound', ...
            ['ub(v_tgp) = %g m/s is far above any physical flow speed, so v_tgp is effectively unbounded and its fitted values will be unreliable. ', ...
             'Velocities are in m/s (50 mm/s = 50e-3). Use Inf if you really want it unbounded.'], max(ubM(:,1)))
    end
    if isscalar(o.stepCap) && isnan(o.stepCap)                % auto
        o.stepCap = 0.5;
        if any(ubM(:,1) > 1), o.stepCap = [10 0.5 0.5 0.5 0.5 0.5]; end      % v_tgp effectively unbounded (ub > 1 m/s or Inf)
    end
    if ~any(numel(o.stepCap) == [1 6]) || any(o.stepCap(:) <= 0) || any(isnan(o.stepCap(:)))
        error('vUS_3D_quad_fitBatched:badStepCap', 'opts.stepCap must be a positive scalar or a positive [1 x 6] vector (Inf allowed).')
    end
    optLM = struct('maxIter', o.maxIter, 'tolF', o.tol, 'tolX', 10*o.tol, 'lam0', o.lam0, 'stepCap', o.stepCap(:));

    x = zeros(n, 6);  cost = zeros(n, 1);  iters = zeros(n, 1);  converged = false(n, 1);
    [~, order] = sort(tdi);                                   % length-sorted chunks keep zero-padding small
    for a0 = 1:o.chunkSize:n
        c = order(a0:min(a0 + o.chunkSize - 1, n));  B = numel(c);
        P = vUS_3D_quad_batchPack(c, g1, tdi, o.t1i, tau, sigma, k0, s, w, o.xscale, o.precision);

        x0 = repmat(x0row(:), 1, B);  x0(2,:) = Vz0(c).';
        LB = cast(lbM(c,:).' ./ o.xscale(:), o.precision);
        UB = cast(ubM(c,:).' ./ o.xscale(:), o.precision);
        X0 = cast(x0 ./ o.xscale(:), o.precision);
        d = min(1e-3*(UB - LB), 1e-3);                            % start strictly inside the box: nudge = 0.1% of the box width, but never more than 1e-3
        X0 = min(max(X0, LB + d), UB - d);                        % (scaled units). An uncapped nudge pushed every start point to 0.1% of a huge ub (50 m/s for ub = 50e3)

        if o.useGPU
            fn = setdiff(fieldnames(P), {'sigma', 'k0'});
            for i = 1:numel(fn), P.(fn{i}) = gpuArray(P.(fn{i})); end
            X0 = gpuArray(X0);  LB = gpuArray(LB);  UB = gpuArray(UB);
        end

        [Xc, cc, ic, cv] = vUS_3D_quad_batchLM(P, X0, LB, UB, optLM);

        x(c, :) = double(gather(Xc .* cast(o.xscale(:), 'like', Xc))).';
        cost(c) = double(gather(cc)).';  iters(c) = double(ic).';  converged(c) = logical(cv).';
    end
end

function M = resolveBounds(u, def, n, name)
% Expand a user bound input ([], [1,6] or [n,6]; NaN = default) to an [n,6] matrix
    if isempty(u), M = def; return, end
    if ~isnumeric(u) || ~ismatrix(u) || size(u, 2) ~= 6 || ~(size(u, 1) == 1 || size(u, 1) == n)
        error('vUS_3D_quad_fitBatched:badBounds', '%s must be [], [1 x 6] or [nVox x 6] (nVox = %d), got size %s.', name, n, mat2str(size(u)))
    end
    M = double(u);
    if size(M, 1) == 1, M = repmat(M, n, 1); end
    isdef = isnan(M);
    M(isdef) = def(isdef);
end
