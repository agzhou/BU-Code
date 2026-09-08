% Description: Estimate a dense 2D tissue-displacement field between two
% image frames (e.g. consecutive |IQ| amplitude or power-Doppler frames)
% via block matching (normalized cross-correlation speckle tracking).
%
% Why block matching instead of optical flow: ultrasound speckle is a
% granular, spatially-correlated pattern rather than a piecewise-smooth
% brightness field, so brightness-constancy-based optical flow tends to
% be noisier here than exploiting the speckle correlation directly. This
% is the same principle as the 1D lateral speckle tracking in
% findVxSpeckleTracking.m, extended to a 2D (axial + lateral) grid of
% blocks instead of a single velocity per depth row.
%
% Operates on ONE frame pair -- loop this over consecutive (or lagged)
% frame pairs externally to build a motion time series. Use
% (xCenters, zCenters, dxMap, dzMap) directly as quiver() inputs for
% visualization.
%
% Inputs:
%   img0: [nz, nx] real reference frame (e.g. abs(IQ(:,:,t)) or a PDI/B-mode frame)
%   img1: [nz, nx] real frame to track relative to img0 (e.g. abs(IQ(:,:,t+lag))),
%       same size as img0
%   dz, dx: axial and lateral pixel spacing [m]
%   blockSize: (optional, default [21, 21]) [bz, bx] template block size [pixels]
%   searchMargin: (optional, default [10, 10]) [mz, mx] max searched
%       displacement on each side of the block [pixels]
%   step: (optional, default = blockSize, i.e. non-overlapping blocks)
%       [sz, sx] spacing between block centers [pixels]
%   minPeakCorr: (optional, default 0.3) minimum normalized cross-correlation
%       peak to trust a block. Blocks below this still get a displacement
%       estimate (returned) but are flagged in `failed`
%
% Outputs:
%   dzMap, dxMap: [nBlocksZ, nBlocksX] estimated axial/lateral displacement
%       of img1 relative to img0, sub-pixel refined [m]. Positive dz = motion
%       toward greater depth; positive dx = motion toward increasing column
%       index. NaN where tracking failed outright (see `failed`).
%   zCenters, xCenters: [nBlocksZ, 1] / [1, nBlocksX] pixel position (row/col
%       index into img0/img1) of each block's center, for use as quiver plot
%       coordinates (multiply by dz/dx for physical units)
%   peakCorr: [nBlocksZ, nBlocksX] normalized cross-correlation peak value at
%       each block (~1 = confident match, low = unreliable/decorrelated)
%   failed: [nBlocksZ, nBlocksX] logical. True where the correlation peak fell
%       on the edge of the search window or the block was degenerate (flat) --
%       dzMap/dxMap are NaN in that case -- OR where peakCorr < minPeakCorr
%       (displacement estimate still returned, just low-confidence)
function [dzMap, dxMap, zCenters, xCenters, peakCorr, failed] = findMotionField2D(img0, img1, dz, dx, blockSize, searchMargin, step, minPeakCorr)
    if nargin < 5 || isempty(blockSize), blockSize = [21, 21]; end
    if nargin < 6 || isempty(searchMargin), searchMargin = [10, 10]; end
    if nargin < 7 || isempty(step), step = blockSize; end
    if nargin < 8 || isempty(minPeakCorr), minPeakCorr = 0.3; end

    [nz, nx] = size(img0);
    bz = blockSize(1); bx = blockSize(2);
    mz = searchMargin(1); mx = searchMargin(2);
    sz = step(1); sx = step(2);

    z0range = (mz+1) : sz : (nz-bz+1-mz);
    x0range = (mx+1) : sx : (nx-bx+1-mx);
    if isempty(z0range) || isempty(x0range)
        error('findMotionField2D:imageTooSmall', 'Image is too small for the given blockSize/searchMargin -- need nz > 2*mz+bz and nx > 2*mx+bx.')
    end

    nBlocksZ = numel(z0range); nBlocksX = numel(x0range);
    dzMap = nan(nBlocksZ, nBlocksX);
    dxMap = nan(nBlocksZ, nBlocksX);
    peakCorr = nan(nBlocksZ, nBlocksX);
    failed = true(nBlocksZ, nBlocksX);
    zCenters = zeros(nBlocksZ, 1);
    xCenters = zeros(1, nBlocksX);

    for bi = 1:nBlocksZ
        z0 = z0range(bi);
        zCenters(bi) = z0 + floor((bz-1)/2);

        for bj = 1:nBlocksX
            x0 = x0range(bj);
            if bi == 1
                xCenters(bj) = x0 + floor((bx-1)/2);
            end

            T = img0(z0:z0+bz-1, x0:x0+bx-1);
            S = img1(z0-mz:z0+bz-1+mz, x0-mx:x0+bx-1+mx);

            if std(T(:)) == 0 || std(S(:)) == 0
                continue % degenerate (flat) block -- leave as NaN/failed
            end

            C = normxcorr2(T, S);
            % Crop to the fully-overlapping ("valid") region only, so the
            % peak location directly indexes displacement in [-mz,mz] x [-mx,mx]
            Cvalid = C(bz:bz+2*mz, bx:bx+2*mx);

            [pk, idx] = max(Cvalid(:));
            [ypeak, xpeak] = ind2sub(size(Cvalid), idx);
            peakCorr(bi, bj) = pk;

            onEdge = (ypeak==1 || ypeak==size(Cvalid,1) || xpeak==1 || xpeak==size(Cvalid,2));
            if onEdge || isnan(pk)
                continue % true displacement likely >= search window; leave as NaN/failed
            end

            % Sub-pixel refine via separable 3-point parabolic fit around the peak
            y1 = Cvalid(ypeak-1, xpeak); y2 = Cvalid(ypeak, xpeak); y3 = Cvalid(ypeak+1, xpeak);
            deltaZ = 0.5*(y1-y3)/(y1-2*y2+y3+eps);
            x1 = Cvalid(ypeak, xpeak-1); x2 = Cvalid(ypeak, xpeak); x3 = Cvalid(ypeak, xpeak+1);
            deltaX = 0.5*(x1-x3)/(x1-2*x2+x3+eps);

            dzMap(bi, bj) = ((ypeak - (mz+1)) + deltaZ) * dz;
            dxMap(bi, bj) = ((xpeak - (mx+1)) + deltaX) * dx;
            failed(bi, bj) = pk < minPeakCorr; % geometrically valid match, but flag low-confidence (e.g. low-texture) blocks
        end
    end
end
