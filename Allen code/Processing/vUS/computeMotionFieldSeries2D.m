% Description: Loop findMotionField2D.m over a stack of frames to build a
% 2D dense tissue-motion field time series.
%
% Inputs:
%   ampStack: [nz, nx, nf] real amplitude images to track (e.g.
%       squeeze(abs(IQ)) -- use the RAW/unfiltered amplitude, not a
%       clutter-filtered blood-only signal like IQf in
%       continuous_2D_testing.m, since the tissue speckle pattern being
%       tracked here is exactly what SVD clutter filtering is designed to
%       remove)
%   dz, dx: axial/lateral pixel spacing [m] (e.g. PData.PDelta(3)*P.wl,
%       PData.PDelta(1)*P.wl)
%   lag: (optional, default 1) frame separation to track [frames] --
%       frame t is tracked against frame t+lag. Larger lag gives larger,
%       easier-to-measure displacements but risks speckle decorrelation
%       between the two tracked frames
%   blockSize, searchMargin, step, minPeakCorr: (optional) passed through
%       to findMotionField2D.m -- see that function for defaults/meaning
%
% Outputs:
%   dzSeries, dxSeries: [nBlocksZ, nBlocksX, nf-lag] displacement field
%       time series [m]. dzSeries(:,:,t)/dxSeries(:,:,t) is the motion
%       from frame t to frame t+lag
%   peakCorrSeries, failedSeries: [nBlocksZ, nBlocksX, nf-lag] diagnostics,
%       see findMotionField2D.m
%   zCenters, xCenters: block-center pixel coordinates (constant across
%       time, since the block grid geometry doesn't change frame to
%       frame) -- see findMotionField2D.m
function [dzSeries, dxSeries, peakCorrSeries, failedSeries, zCenters, xCenters] = computeMotionFieldSeries2D(ampStack, dz, dx, lag, blockSize, searchMargin, step, minPeakCorr)
    if nargin < 4 || isempty(lag), lag = 1; end
    if nargin < 5, blockSize = []; end
    if nargin < 6, searchMargin = []; end
    if nargin < 7, step = []; end
    if nargin < 8, minPeakCorr = []; end

    nf = size(ampStack, 3);
    nPairs = nf - lag;
    if nPairs < 1
        error('computeMotionFieldSeries2D:notEnoughFrames', 'ampStack has %d frames, which is not enough for lag = %d.', nf, lag)
    end

    % Run the first pair to size the outputs from the actual block grid
    [dz1, dx1, zCenters, xCenters, pk1, fail1] = findMotionField2D(ampStack(:,:,1), ampStack(:,:,1+lag), dz, dx, blockSize, searchMargin, step, minPeakCorr);
    [nBlocksZ, nBlocksX] = size(dz1);

    dzSeries = nan(nBlocksZ, nBlocksX, nPairs);
    dxSeries = nan(nBlocksZ, nBlocksX, nPairs);
    peakCorrSeries = nan(nBlocksZ, nBlocksX, nPairs);
    failedSeries = true(nBlocksZ, nBlocksX, nPairs);

    dzSeries(:,:,1) = dz1; dxSeries(:,:,1) = dx1;
    peakCorrSeries(:,:,1) = pk1; failedSeries(:,:,1) = fail1;

    for t = 2:nPairs
        [dzSeries(:,:,t), dxSeries(:,:,t), ~, ~, peakCorrSeries(:,:,t), failedSeries(:,:,t)] = ...
            findMotionField2D(ampStack(:,:,t), ampStack(:,:,t+lag), dz, dx, blockSize, searchMargin, step, minPeakCorr);
    end
end
