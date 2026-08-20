% Description: Extract the anatomical orientation of each vessel in a 3D
% structural fUS/vUS volume (e.g. a 3D power Doppler volume). This is the
% 3D counterpart of vesselAngle2D.m -- same pipeline, same reasoning,
% generalized from a 2D line orientation (one angle, mod 180 deg) to a 3D
% line orientation (a unit axis, mod a sign flip). It is purely
% image-based (Hessian ridge detection of the vessel walls) and
% independent of any flow-velocity fit.
%
% Pipeline: vesselnessAngle3D.m (Jerman-filter vesselness + per-voxel
% tangent direction, see that file) -> threshold to a binary vessel mask
% -> skeletonize (bwskel) -> cut the skeleton at branch points -> for
% each resulting vessel segment, average the pixelwise tangent direction
% into one representative axis.
%
% Averaging a 3D line direction is NOT a plain vector mean: each voxel's
% tangent is only defined up to a sign flip (dir ~ -dir), and a vector
% mean of, say, [1,0,0] and [-1,0,0] would wrongly cancel to zero even
% though both describe the same axis. Instead this uses the standard
% orientation-tensor trick (the same one behind DTI principal-direction
% estimation): average the SIGN-INVARIANT outer product dir*dir' over the
% segment's voxels, weighted by vesselness, then take the eigenvector of
% the largest eigenvalue of that averaged 3x3 tensor as the mean axis.
% "Straightness" is Westin's tensor linearity measure c_l = lambda1 -
% lambda2 (eigenvalues sorted descending, summing to 1 since every input
% vector is unit length): 1 for a segment where every voxel agrees on one
% axis, 0 where the tensor is isotropic (no consistent axis at all).
%
% Inputs:
%   I: 3D structural image [dim1, dim2, dim3] (e.g. a 3D PDI volume).
%      Vessels should be the brightest structures in I.
%   sigmas: vector of vessel radii to search over [voxels, or spacing
%       units if spacing ~= 1]. Should roughly span the range of vessel
%       calibers present in I.
%   spacing: (optional) [dim1; dim2; dim3] voxel spacing, passed through
%       to vesselnessAngle3D.m for anisotropic voxels. Default [1; 1; 1].
%   tau: (optional) Jerman filter response-uniformity parameter, in
%       (0.5, 1); lower = more intense response. Default 1.
%   brightondark: (optional) true if vessels are brighter than the
%       background (true for PDI). Default true.
%   vesselnessThreshold: (optional) threshold on the normalized [0, 1]
%       vesselness map used to build the binary vessel mask. Default 0.1.
%   minBranchLengthPix: (optional) minimum skeleton branch length [voxels]
%       passed to bwskel for spur pruning. Default 5.
%   minSegLengthPix: (optional) minimum number of skeleton voxels for a
%       connected piece between branch points to be reported as a vessel
%       segment. Default 5.
%
% Outputs:
%   vessels: 1 x N struct array, one entry per vessel segment (a piece of
%       skeleton between branch/end points), with fields:
%       .PixelIdxList: linear indices (into size(I)) of the segment's skeleton voxels
%       .coordsDim123: [# voxels, 3], (dim1, dim2, dim3) coordinates of those voxels
%       .lengthPix: number of skeleton voxels in the segment (NOT
%           corrected for diagonal steps -- a rough proxy for length, not
%           a true arclength)
%       .direction: [1, 3] unit vector, the segment's mean tangent axis
%           (dim1, dim2, dim3 components). Sign is canonicalized so its
%           largest-magnitude component is positive, for a reproducible
%           (if still ultimately arbitrary) reporting convention -- it is
%           still a LINE direction, not a flow direction.
%       .azimuthDeg: atan2d(direction(2), direction(1)) -- angle in the
%           dim1/dim2 plane, [-180, 180]. A convenience scalar for
%           reporting/plotting only; use .direction for computation.
%       .elevationDeg: asind(direction(3)) -- angle out of the dim1/dim2
%           plane toward dim3, [-90, 90].
%       .straightness: Westin linearity c_l = lambda1 - lambda2 of the
%           segment's orientation tensor, in [0, 1] (1 = every voxel
%           along the segment agrees on orientation; lower = curved
%           and/or noisy)
%       .meanVesselness: mean vesselness value along the segment
%   dir1, dir2, dir3: [size(I)], pixelwise tangent direction components
%       from vesselnessAngle3D.m (NaN outside detected vessel structure)
%   vesselness: [size(I)], normalized [0, 1] vesselness map
%   vesselMask: [size(I)], logical, vesselness > vesselnessThreshold
%   segLabel: [size(I)], integer label volume; segLabel == k marks the
%       voxels of vessels(k). 0 = not part of any reported segment.
%   skel: [size(I)], logical, the full pruned skeleton (before cutting at
%       branch points)
%
% Requires vesselnessAngle3D.m (Allen code/Processing/Jerman Enhancement
% Filter/, which itself requires the compiled eig3volume.c/.mexw64 in
% that same folder) and Image Processing Toolbox (bwskel, bwconncomp,
% imdilate) on the MATLAB path.
%
% Example:
%   PDI3D = squeeze(mean(abs(IQf_HPF).^2, 4));
%   [vessels, dir1, dir2, dir3] = vesselAngle3D(PDI3D, 1:0.5:4);
function [vessels, dir1, dir2, dir3, vesselness, vesselMask, segLabel, skel] = vesselAngle3D(I, sigmas, spacing, tau, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix)

    if nargin < 3 || isempty(spacing), spacing = [1; 1; 1]; end
    if nargin < 4 || isempty(tau), tau = 1; end
    if nargin < 5 || isempty(brightondark), brightondark = true; end
    if nargin < 6 || isempty(vesselnessThreshold), vesselnessThreshold = 0.1; end
    if nargin < 7 || isempty(minBranchLengthPix), minBranchLengthPix = 5; end
    if nargin < 8 || isempty(minSegLengthPix), minSegLengthPix = 5; end

    % ===== 1. Ridge (vesselness) detection + per-voxel tangent direction ===== %
    [vesselness, dir1, dir2, dir3] = vesselnessAngle3D(I, sigmas, spacing, tau, brightondark);

    % ===== 2. Threshold + skeletonize ===== %
    vesselMask = vesselness > vesselnessThreshold;
    skel = bwskel(vesselMask, 'MinBranchLength', minBranchLengthPix);

    % Cut the skeleton at branch points so each connected piece left over
    % is a single vessel segment between two branch/end points.
    % bwmorph (used for this in 2D) has no 3D support, so branch points
    % are found directly: a skeleton voxel with >= 3 skeleton neighbors
    % in its 26-connected neighborhood is a branch point.
    nbrCount = convn(double(skel), ones(3,3,3), 'same') - double(skel);
    branchPts = skel & (nbrCount >= 3);
    segSkel = skel & ~imdilate(branchPts, ones(3,3,3));

    % ===== 3. Per-segment direction: vesselness-weighted orientation tensor ===== %
    CC = bwconncomp(segSkel);
    segLabel = zeros(size(I));

    vessels = struct('PixelIdxList', {}, 'coordsDim123', {}, 'lengthPix', {}, 'direction', {}, 'azimuthDeg', {}, 'elevationDeg', {}, 'straightness', {}, 'meanVesselness', {});
    nSeg = 0;
    for k = 1:CC.NumObjects
        idx = CC.PixelIdxList{k};
        if numel(idx) < minSegLengthPix
            continue
        end
        nSeg = nSeg + 1;
        segLabel(idx) = nSeg;

        [d1, d2, d3] = ind2sub(size(I), idx);
        v = [dir1(idx), dir2(idx), dir3(idx)]; % [# voxels, 3]
        w = vesselness(idx);
        if sum(w) == 0
            w = ones(size(w)); % guard against a zero-vesselness segment (shouldn't normally happen post-threshold)
        end
        w = w / sum(w);

        % Sign-invariant orientation tensor: weighted mean of dir*dir'
        T = (v .* w)' * v; % 3x3, symmetric since it's a weighted sum of outer products

        [V, D] = eig(T, 'vector');
        [D, order] = sort(D, 'descend'); % lambda1 >= lambda2 >= lambda3, sum(D) = 1
        V = V(:, order);

        meanDir = V(:, 1)';
        [~, signIdx] = max(abs(meanDir));
        if meanDir(signIdx) < 0
            meanDir = -meanDir; % canonicalize sign for reproducible reporting
        end

        vessels(nSeg).PixelIdxList = idx;
        vessels(nSeg).coordsDim123 = [d1, d2, d3];
        vessels(nSeg).lengthPix = numel(idx);
        vessels(nSeg).direction = meanDir;
        vessels(nSeg).azimuthDeg = atan2d(meanDir(2), meanDir(1));
        vessels(nSeg).elevationDeg = asind(meanDir(3));
        vessels(nSeg).straightness = D(1) - D(2);
        vessels(nSeg).meanVesselness = mean(vesselness(idx));
    end
end
