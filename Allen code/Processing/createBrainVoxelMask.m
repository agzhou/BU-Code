%% Description: automatically segment a 3D volume into brain-tissue voxels vs.
%   everything above the brain surface (coupling medium, skull/dura, cranial
%   window), by detecting the shallow specular reflection at that interface.
%
%   IMPORTANT: use the coherent sum of RAW IQ data as the anatomical reference
%   (e.g. sum(abs(IQ), 4), summed BEFORE any SVD clutter filtering/declutter).
%   The clutter filter intentionally suppresses static reflectors like the
%   skull/coupling interface, so a post-clutter-filter image (PDI/CDI) will NOT
%   show the surface clearly and should not be used here.
%
%   Intended use: mask out non-brain voxels (by zeroing them in IQ) BEFORE
%   calling getSVs2D/applySVs2D_accel, so the SVD's covariance matrix (and
%   therefore the clutter/blood-flow eigenbasis) is built only from brain-tissue
%   dynamics, not dominated by the much higher-energy coupling-medium/skull
%   reflection. Also apply the same mask to the final PDI/CDI/activation maps as
%   a belt-and-suspenders spatial exclusion.
%
% Inputs:
%   anatVol: anatomical reference volume. One of:
%       - (x, y, z) real/complex array: used directly (abs() is taken)
%       - (x, y, z, frames) real/complex array: reduced via sum(abs(.), 4)
%       - cell array of either of the above (e.g. a {PDI/IQ per component} cell,
%         such as PDI{3} / IQ{3}): the LAST cell element is used
%   Name-value inputs:
%       - 'ZDim' (default 3): the depth dimension of anatVol
%       - 'SearchRange' [z1, z2] (default [1, round(0.6*nz)]): z-index range to
%         search for the specular skull/coupling interface, so a deep bright
%         vessel isn't mistaken for the surface
%       - 'SmoothSpatialSigma' (default 2): Gaussian sigma [voxels] for smoothing
%         the per-(x,y) surface depth map across x,y (robust to probe tilt /
%         noisy single-column detections)
%       - 'MarginVoxels' (default 4): extra voxels to exclude below the detected
%         surface, since the specular reflection has some axial ringdown/PSF width
%       - 'Method' ('gradient' default | 'threshold'):
%           'gradient':  per (x,y) column, find the steepest rising edge (max
%                        positive first difference) within SearchRange -- robust
%                        to the absolute gain/TGC curve, since it looks for the
%                        edge, not an absolute brightness level
%           'threshold': per (x,y) column, find the first z within SearchRange
%                        where intensity crosses ThresholdPercentile of the
%                        volume's intensity distribution
%       - 'ThresholdPercentile' (default 75): for 'threshold' method only
%       - 'MinEdgeConfidence' (default 0.15): columns whose strongest edge is
%         weaker than this fraction of that column's peak intensity are treated
%         as unreliable and replaced with the (spatially smoothed) median surface
%         instead of their own (noisy) detection
%       - 'Visualize' (default true): produce a diagnostic figure (coronal +
%         sagittal MIPs with the detected surface overlaid, and the resulting mask)
%
% Outputs:
%   brainMask: logical (x, y, z) array; true = brain tissue to keep
%   info: struct with
%       - surfaceZ: (x, y) map of the detected surface z-index (pre-margin)
%       - confident: (x, y) logical map of which columns had a confident direct
%         detection vs. a filled-in value
%       - anatMIP: struct with coronal/sagittal max-intensity projections, for
%         your own plotting/QC
%       - params: the resolved options used
%
% Example (recommended usage, masking BEFORE the clutter filter):
%   load([IQpath, IQfilenameStructure, '1'])
%   IQ = squeeze(IData + 1i .* QData);
%   anatVol = sum(abs(IQ), 4);
%   [brainMask, info] = createBrainVoxelMask(anatVol);
%   IQ_masked = IQ .* brainMask;  % zero out non-brain voxels before SVD
%   [CM, EVs, V] = getSVs2D(IQ_masked(xrange, yrange, zrange, :));
%   ... % continue as in IQ2g1T_3D.m, and also apply brainMask(xrange,yrange,zrange)
%       % to the resulting PDI/CDI/g1 before downstream stats

function [brainMask, info] = createBrainVoxelMask(anatVol, opts)
    arguments
        anatVol
        opts.ZDim (1, 1) double = 3
        opts.SearchRange (1, 2) double = [NaN, NaN]
        opts.SmoothSpatialSigma (1, 1) double = 2
        opts.MarginVoxels (1, 1) double {mustBeNonnegative} = 4
        opts.Method (1, :) char = 'gradient'
        opts.ThresholdPercentile (1, 1) double = 75
        opts.MinEdgeConfidence (1, 1) double = 0.15
        opts.Visualize (1, 1) logical = true
    end

    %% Resolve the anatomical reference to a real-valued (x, y, z) volume
    if iscell(anatVol)
        anatVol = anatVol{end};
    end
    anatVol = abs(anatVol);
    if ndims(anatVol) == 4
        anatVol = sum(anatVol, 4);
    elseif ndims(anatVol) > 4
        error('createBrainVoxelMask:badInput', 'anatVol must be (x,y,z) or (x,y,z,frames).')
    end

    % Move the depth dimension to the end for convenience, then move back at the end
    nd = ndims(anatVol);
    order = [setdiff(1:nd, opts.ZDim), opts.ZDim];
    A = permute(anatVol, order); % now (x, y, ..., z)
    sz = size(A);
    nz = sz(end);
    xy = A(:, :, :); % collapse any leading dims into one "column" dimension if nd > 3
    xy = reshape(A, [], nz); % (nCols, nz)

    searchRange = opts.SearchRange;
    if any(isnan(searchRange))
        searchRange = [1, max(round(0.6 * nz), 2)];
    end
    z1 = max(1, round(searchRange(1)));
    z2 = min(nz, round(searchRange(2)));
    if z2 - z1 < 1
        error('createBrainVoxelMask:badRange', 'SearchRange must span at least 2 z-indices within 1:%d.', nz)
    end

    %% Per-column surface detection
    nCols = size(xy, 1);
    surfaceZ = nan(nCols, 1);
    confidence = zeros(nCols, 1);

    switch validatestring(opts.Method, {'gradient', 'threshold'})
        case 'gradient'
            seg = xy(:, z1:z2);
            d = diff(seg, 1, 2); % (nCols, z2-z1) first difference along z
            [edgeVal, edgeIdx] = max(d, [], 2);
            surfaceZ = z1 + edgeIdx - 1 + 0.5; % edge sits between sample edgeIdx and edgeIdx+1
            colPeak = max(seg, [], 2);
            colPeak(colPeak == 0) = eps;
            confidence = edgeVal ./ colPeak; % normalized edge strength, ~0-1
        case 'threshold'
            thresh = prctile(xy(:), opts.ThresholdPercentile);
            seg = xy(:, z1:z2) >= thresh;
            [hasCrossing, firstIdx] = max(seg, [], 2); % first true index (or 1 if none found)
            surfaceZ = z1 + firstIdx - 1;
            surfaceZ(~hasCrossing) = NaN;
            confidence = double(hasCrossing);
    end

    %% Fill in / smooth unreliable columns using a spatially-smoothed consensus
    spatialSz = sz(1:end-1);
    surfaceZimg = reshape(surfaceZ, spatialSz);
    confImg = reshape(confidence, spatialSz);

    reliable = confImg >= opts.MinEdgeConfidence & ~isnan(surfaceZimg);
    if ~any(reliable, 'all')
        error('createBrainVoxelMask:noReliableSurface', ...
            'No columns had a confident surface detection. Try a different SearchRange, Method, or lower MinEdgeConfidence.')
    end

    % Spatially smooth the reliable detections, then fill unreliable columns by
    % nearest-neighbor extrapolation from the smoothed reliable map
    filled = surfaceZimg;
    filled(~reliable) = NaN;
    if opts.SmoothSpatialSigma > 0 && numel(spatialSz) >= 2
        filled = fillmissing2(filled, 'movmedian', max(3, round(4 * opts.SmoothSpatialSigma)));
        filled = imgaussfilt(filled, opts.SmoothSpatialSigma, 'FilterDomain', 'spatial');
    else
        filled = fillmissing(filled, 'linear', 'EndValues', 'nearest');
    end
    surfaceZfinal = filled;

    %% Build the mask: brain = below (deeper than) the detected surface + margin
    cutoffZ = surfaceZfinal + opts.MarginVoxels;
    zGrid = reshape(1:nz, [ones(1, numel(spatialSz)), nz]);
    maskPerm = zGrid > cutoffZ; % broadcast (spatialSz x nz) > (spatialSz x 1)

    brainMask = ipermute(maskPerm, order);

    %% Info / diagnostics
    info.surfaceZ = surfaceZfinal;
    info.rawSurfaceZ = surfaceZimg;
    info.confident = reliable;
    info.cutoffZ = cutoffZ;
    info.params = opts;
    info.params.SearchRange = searchRange;

    if numel(spatialSz) >= 2 % only meaningful for true 3D (x,y,z) input
        info.anatMIP.coronal = squeeze(max(A, [], 1)); % (y, z)
        info.anatMIP.sagittal = squeeze(max(A, [], 2)); % (x, z)
    end

    %% Visualization
    if opts.Visualize && numel(spatialSz) >= 2
        f = figure('Name', 'createBrainVoxelMask diagnostic');
        subplot(2, 2, 1)
        imagesc(info.anatMIP.coronal.' .^ 0.5); axis image; colormap(gca, 'hot'); colorbar
        hold on; plot(1:spatialSz(2), surfaceZfinal(round(spatialSz(1)/2), :), 'c-', 'LineWidth', 1.5); hold off
        title('Coronal MIP (raw IQ) with detected surface'); xlabel('y'); ylabel('z')

        subplot(2, 2, 2)
        imagesc(info.anatMIP.sagittal.' .^ 0.5); axis image; colormap(gca, 'hot'); colorbar
        hold on; plot(1:spatialSz(1), surfaceZfinal(:, round(spatialSz(2)/2)), 'c-', 'LineWidth', 1.5); hold off
        title('Sagittal MIP (raw IQ) with detected surface'); xlabel('x'); ylabel('z')

        subplot(2, 2, 3)
        imagesc(squeeze(max(maskPerm, [], 1)).'); axis image; colormap(gca, 'gray'); colorbar
        title('Brain mask (coronal MIP)'); xlabel('y'); ylabel('z')

        subplot(2, 2, 4)
        imagesc(reliable.'); axis image; colormap(gca, 'gray'); colorbar
        title(sprintf('Confident direct detections (%.0f%% of columns)', 100*mean(reliable(:))))
        xlabel('x'); ylabel('y')
    end
end
