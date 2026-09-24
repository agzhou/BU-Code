function [x, cost, iters, converged] = vUS_3D_quad_fitBatched(g1, tdi, Vz0, tau, k0, sigma, opts)
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
% Inputs:
%   g1: [nVox, nTau] complex g1 of the voxels to fit, tau starting at 0
%   tdi: [nVox, 1] index of the last lag to fit per voxel (findTauDecayed's tau_decayed_ind);
%       the fit window is opts.t1i:tdi
%   Vz0: [nVox, 1] initial axial group velocity [m/s] (findVzPhaseDiff); v_zgp is bounded to Vz0 +/- opts.vzHalfWidth
%   tau: [nTau, 1] time lags [s];  k0: wavenumber [rad/m]
%   sigma: [sigma_x, sigma_y, sigma_z] [m]; sigma_x must equal sigma_y
%   opts: (optional) struct; any field may be omitted
%       t1i: first lag index to fit (2)
%       nNodes: Gauss-Legendre nodes (48; 24 was accurate to 4e-9 on real fit windows)
%       chunkSize: voxels per batch (1024 was fastest on CPU)
%       precision: 'double' (default) or 'single' (~2x faster on CPU; slightly more voxels in worse minima)
%       tol: convergence tolerance (1e-10 relative cost change; scaled step tolerance is 10*tol; use ~1e-5 for single)
%       maxIter: iteration cap per voxel (300)
%       lam0: initial Levenberg-Marquardt damping (1);  stepCap: max scaled step per iteration (0.5)
%       xscale: [1,6] parameter scale (1e-2 1e-2 1 1 1 1: v_tgp, v_zgp in units of 10 mm/s)
%       x0: [1,6] start point, NaN in position 2 = use Vz0 ([5e-3 NaN 1 0 2.5 1])
%       lb, ub: [1,6] bounds, NaN in position 2 = Vz0 -/+ vzHalfWidth ([0 NaN 0 0 2 0], [sqrt(2)*30e-3 NaN 1 1 3 1])
%       vzHalfWidth: half-width of the v_zgp bound around Vz0 [m/s] (0.01)
%       useGPU: run the chunks on a GPU via gpuArray (false). UNTESTED -- no GPU was available when written
%
% Outputs:
%   x: [nVox, 6] fitted [v_tgp, v_zgp, F, DC, k, a] in physical units (m/s for the velocities)
%   cost: [nVox, 1] final sum of squared residuals over the fit window (lsqnonlin's resnorm)
%   iters: [nVox, 1] iterations used
%   converged: [nVox, 1] false where the iteration cap was hit

    o = struct('t1i', 2, 'nNodes', 48, 'chunkSize', 1024, 'precision', 'double', 'tol', 1e-10, 'maxIter', 300, ...
        'lam0', 1, 'stepCap', 0.5, 'xscale', [1e-2 1e-2 1 1 1 1], 'x0', [5e-3 NaN 1 0 2.5 1], ...
        'lb', [0 NaN 0 0 2 0], 'ub', [sqrt(2)*30e-3 NaN 1 1 3 1], 'vzHalfWidth', 0.01, 'useGPU', false);
    if nargin > 6 && ~isempty(opts)
        fn = fieldnames(opts);
        for i = 1:numel(fn)
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
    lbrow = o.lb; ubrow = o.ub; x0row = o.x0;
    optLM = struct('maxIter', o.maxIter, 'tolF', o.tol, 'tolX', 10*o.tol, 'lam0', o.lam0, 'stepCap', o.stepCap);

    x = zeros(n, 6);  cost = zeros(n, 1);  iters = zeros(n, 1);  converged = false(n, 1);
    [~, order] = sort(tdi);                                   % length-sorted chunks keep zero-padding small
    for a0 = 1:o.chunkSize:n
        c = order(a0:min(a0 + o.chunkSize - 1, n));  B = numel(c);
        P = vUS_3D_quad_batchPack(c, g1, tdi, o.t1i, tau, sigma, k0, s, w, o.xscale, o.precision);

        vz = Vz0(c).';
        x0 = repmat(x0row(:), 1, B);  x0(2,:) = vz;
        lb = repmat(lbrow(:), 1, B);  lb(2,:) = vz - o.vzHalfWidth;
        ub = repmat(ubrow(:), 1, B);  ub(2,:) = vz + o.vzHalfWidth;
        LB = cast(lb ./ o.xscale(:), o.precision);  UB = cast(ub ./ o.xscale(:), o.precision);  X0 = cast(x0 ./ o.xscale(:), o.precision);
        d = 1e-3*(UB - LB);  X0 = min(max(X0, LB + d), UB - d);   % start strictly inside the box

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
