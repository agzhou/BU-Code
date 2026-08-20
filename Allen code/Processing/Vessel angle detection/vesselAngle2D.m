% Description: Extract the anatomical orientation of each vessel in a 2D
% structural fUS/vUS image (e.g. the power Doppler image, PDI). This is
% purely image-based (Hessian ridge detection of the vessel walls) and is
% independent of any flow-velocity fit -- it answers "which way does this
% vessel run", not "which way is the blood moving".
%
% Pipeline: vesselnessAngle2D.m (Jerman-filter vesselness + per-pixel
% tangent angle, see that file) -> threshold to a binary vessel mask ->
% skeletonize (bwskel) -> cut the skeleton at branch points -> for each
% resulting vessel segment, take the vesselness-weighted circular mean of
% the pixelwise tangent angle (correctly handling the mod-180-degree
% wraparound of an undirected line orientation).
%
% Inputs:
%   I: 2D structural image [dim1, dim2] (e.g. PDI = squeeze(mean(abs(IQf_HPF).^2, 3))
%      from vUS_2D.m). Vessels should be the brightest structures in I.
%   sigmas: vector of vessel radii to search over [pixels, or spacing units
%       if spacing ~= 1]. Should roughly span the range of vessel calibers
%       present in I.
%   spacing: (optional) [dim1; dim2] pixel spacing, passed through to
%       vesselnessAngle2D.m for anisotropic pixels. Default [1; 1].
%   tau: (optional) Jerman filter response-uniformity parameter, in
%       (0.5, 1); lower = more intense response. Default 1.
%   brightondark: (optional) true if vessels are brighter than the
%       background (true for PDI). Default true.
%   vesselnessThreshold: (optional) threshold on the normalized [0, 1]
%       vesselness map used to build the binary vessel mask. Default 0.1.
%   minBranchLengthPix: (optional) minimum skeleton branch length [pixels]
%       passed to bwskel for spur pruning. Default 5.
%   minSegLengthPix: (optional) minimum number of skeleton pixels for a
%       connected piece between branch points to be reported as a vessel
%       segment. Default 5.
%
% Outputs:
%   vessels: 1 x N struct array, one entry per vessel segment (a piece of
%       skeleton between branch/end points), with fields:
%       .PixelIdxList: linear indices (into size(I)) of the segment's skeleton pixels
%       .coordsDim1Dim2: [# pixels, 2], (dim1, dim2) coordinates of those pixels
%       .lengthPix: number of skeleton pixels in the segment (NOT
%           corrected for diagonal steps -- a rough proxy for length, not
%           a true arclength)
%       .angleDeg: vesselness-weighted circular mean of the pixelwise
%           tangent angle along the segment [degrees], wrapped to
%           [-90, 90) (see vesselnessAngle2D.m for the angle convention)
%       .straightness: resultant length in [0, 1] of the circular mean
%           (1 = every pixel along the segment agrees on orientation, i.e.
%           dead straight; closer to 0 = the segment's local orientation
%           varies a lot, e.g. it's curved or the estimate is noisy)
%       .meanVesselness: mean vesselness value along the segment
%   angleMap: [size(I)], pixelwise tangent angle [degrees] from
%       vesselnessAngle2D.m (NaN outside detected vessel structure)
%   vesselness: [size(I)], normalized [0, 1] vesselness map
%   vesselMask: [size(I)], logical, vesselness > vesselnessThreshold
%   segLabel: [size(I)], integer label image; segLabel == k marks the
%       pixels of vessels(k). 0 = not part of any reported segment.
%   skel: [size(I)], logical, the full pruned skeleton (before cutting at
%       branch points)
%
% Requires vesselnessAngle2D.m (Allen code/Processing/Jerman Enhancement
% Filter/) and Image Processing Toolbox (bwskel, bwmorph, bwconncomp) on
% the MATLAB path.
%
% Example:
%   PDI = squeeze(mean(abs(IQf_HPF).^2, 3));
%   [vessels, angleMap] = vesselAngle2D(PDI, 1:0.5:4);
function [vessels, angleMap, vesselness, vesselMask, segLabel, skel] = vesselAngle2D(I, sigmas, spacing, tau, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix)

    if nargin < 3 || isempty(spacing), spacing = [1; 1]; end
    if nargin < 4 || isempty(tau), tau = 1; end
    if nargin < 5 || isempty(brightondark), brightondark = true; end
    if nargin < 6 || isempty(vesselnessThreshold), vesselnessThreshold = 0.1; end
    if nargin < 7 || isempty(minBranchLengthPix), minBranchLengthPix = 5; end
    if nargin < 8 || isempty(minSegLengthPix), minSegLengthPix = 5; end

    % ===== 1. Ridge (vesselness) detection + per-pixel tangent angle ===== %
    [vesselness, angleMap] = vesselnessAngle2D(I, sigmas, spacing, tau, brightondark);

    % ===== 2. Threshold + skeletonize ===== %
    vesselMask = vesselness > vesselnessThreshold;
    skel = bwskel(vesselMask, 'MinBranchLength', minBranchLengthPix);

    % Cut the skeleton at branch points so each connected piece left over
    % is a single vessel segment between two branch/end points
    branchPts = bwmorph(skel, 'branchpoints');
    segSkel = skel & ~imdilate(branchPts, strel('square', 3));

    % ===== 3. Per-segment angle: vesselness-weighted circular mean ===== %
    % Orientation is a line (undirected, mod 180 deg), so the standard
    % circular mean would wrap incorrectly (e.g. -89 deg and +89 deg are
    % nearly parallel, not opposite). Fix: double the angle before taking
    % the circular mean, then halve the result -- this is the standard
    % trick for averaging axial (mod-180) data.
    CC = bwconncomp(segSkel);
    segLabel = zeros(size(I), 'like', double(1));

    vessels = struct('PixelIdxList', {}, 'coordsDim1Dim2', {}, 'lengthPix', {}, 'angleDeg', {}, 'straightness', {}, 'meanVesselness', {});
    nSeg = 0;
    for k = 1:CC.NumObjects
        idx = CC.PixelIdxList{k};
        if numel(idx) < minSegLengthPix
            continue
        end
        nSeg = nSeg + 1;
        segLabel(idx) = nSeg;

        [d1, d2] = ind2sub(size(I), idx);
        theta = angleMap(idx); % degrees, in [-90, 90)
        w = vesselness(idx);
        if sum(w) == 0
            w = ones(size(w)); % guard against a zero-vesselness segment (shouldn't normally happen post-threshold)
        end

        phi = 2 * theta * pi/180; % double the angle to remove the mod-180 ambiguity
        Cx = sum(w .* cos(phi)) / sum(w);
        Cy = sum(w .* sin(phi)) / sum(w);
        meanTheta = mod(atan2d(Cy, Cx)/2 + 90, 180) - 90; % undo the doubling, re-wrap to [-90, 90)
        straightness = sqrt(Cx^2 + Cy^2); % resultant length in [0, 1]

        vessels(nSeg).PixelIdxList = idx;
        vessels(nSeg).coordsDim1Dim2 = [d1, d2];
        vessels(nSeg).lengthPix = numel(idx);
        vessels(nSeg).angleDeg = meanTheta;
        vessels(nSeg).straightness = straightness;
        vessels(nSeg).meanVesselness = mean(w);
    end
end
