% Description: Render a 2D dense tissue-motion field time series (from
% computeMotionFieldSeries2D.m) as a video -- each frame shows a
% grayscale anatomical background image with the block-wise motion
% vectors overlaid as a quiver plot.
%
% Inputs:
%   outFilename: output video file path (e.g. 'motion.mp4'); extension
%       selects the VideoWriter profile ('.mp4' -> 'MPEG-4', otherwise
%       'Motion JPEG AVI')
%   bgStack: [nz, nx, nPairs] real background image to display per
%       frame (e.g. squeeze(abs(IQ(:,:,1:end-lag))) or a PDI stack) --
%       same [nz, nx] pixel grid the motion was tracked on, and the same
%       number of frames as dzSeries/dxSeries (nPairs = nf - lag)
%   zCenters, xCenters, dzSeries, dxSeries, failedSeries: outputs of
%       computeMotionFieldSeries2D.m (dz/dxSeries in meters)
%   dz, dx: axial/lateral pixel spacing [m], to convert dzSeries/dxSeries
%       back into pixel units matching bgStack's pixel grid
%   fps: (optional, default 15) video PLAYBACK frame rate [Hz] -- this is
%       independent of the acquisition frame rate, so N seconds of e.g.
%       800 Hz data will play back much faster than real time unless you
%       subsample the series first (or set fps to something like
%       P.frameRate/lag for real-time playback)
%   vectorScale: (optional, default 20) multiplies displacement before
%       drawing, purely for on-screen visibility -- sub-pixel tissue
%       motion is otherwise too small to see as an arrow. A red scale-bar
%       arrow of known physical length is drawn each frame so the
%       exaggeration doesn't mislead
%   clim: (optional, default = [0, 99th percentile] of bgStack) color
%       limits for the background image
%
% Outputs: none (writes outFilename to disk)
function writeMotionFieldVideo2D(outFilename, bgStack, zCenters, xCenters, dzSeries, dxSeries, failedSeries, dz, dx, fps, vectorScale, clim)
    if nargin < 10 || isempty(fps), fps = 15; end
    if nargin < 11 || isempty(vectorScale), vectorScale = 20; end
    if nargin < 12 || isempty(clim), clim = [0, prctile(bgStack(:), 99)]; end

    nPairs = size(dzSeries, 3);
    if size(bgStack, 3) ~= nPairs
        error('writeMotionFieldVideo2D:sizeMismatch', 'bgStack has %d frames but dzSeries/dxSeries have %d -- they must match (bgStack(:,:,t) should be the reference frame for dzSeries(:,:,t)).', size(bgStack, 3), nPairs)
    end

    [~, ~, ext] = fileparts(outFilename);
    if strcmpi(ext, '.mp4')
        vidProfile = 'MPEG-4';
    else
        vidProfile = 'Motion JPEG AVI';
    end
    vw = VideoWriter(outFilename, vidProfile);
    vw.FrameRate = fps;
    open(vw)

    % Pick a round physical scale-bar length near the typical displacement magnitude
    dispMag = sqrt(dzSeries.^2 + dxSeries.^2);
    typicalDispMag = median(dispMag(~failedSeries), 'omitnan');
    if isempty(typicalDispMag) || isnan(typicalDispMag) || typicalDispMag == 0
        scaleBarLengthM = 10e-6; % fallback: 10 um
    else
        scaleBarLengthM = 10^floor(log10(typicalDispMag));
    end

    [X, Z] = meshgrid(xCenters, zCenters);
    zRange = max(zCenters) - min(zCenters);
    sbX0 = min(xCenters);
    sbZ0 = max(zCenters) + 0.08*zRange;

    fig = figure('Visible', 'off', 'Color', 'w');
    ax = axes(fig);

    for t = 1:nPairs
        cla(ax)
        imagesc(ax, bgStack(:,:,t), clim); colormap(ax, 'gray'); axis(ax, 'image'); hold(ax, 'on')

        dzPix = dzSeries(:,:,t) / dz; dxPix = dxSeries(:,:,t) / dx;
        dzPix(failedSeries(:,:,t)) = NaN; dxPix(failedSeries(:,:,t)) = NaN;

        quiver(ax, X, Z, dxPix*vectorScale, dzPix*vectorScale, 0, 'y', 'LineWidth', 1.2, 'MaxHeadSize', 0.5)

        % Scale-bar arrow, drawn in the same (pixel * vectorScale) units as the motion vectors
        sbLenPix = scaleBarLengthM/dx * vectorScale;
        quiver(ax, sbX0, sbZ0, sbLenPix, 0, 0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 1)
        text(ax, sbX0, sbZ0 + 0.06*zRange, sprintf('%.1f \\mum (%dx)', scaleBarLengthM*1e6, vectorScale), 'Color', 'r', 'FontSize', 9)

        ylim(ax, [min(zCenters) - 0.05*zRange, sbZ0 + 0.12*zRange])
        title(ax, sprintf('Frame pair %d / %d', t, nPairs))
        hold(ax, 'off')
        drawnow

        writeVideo(vw, getframe(fig))
    end

    close(vw)
    close(fig)
end
