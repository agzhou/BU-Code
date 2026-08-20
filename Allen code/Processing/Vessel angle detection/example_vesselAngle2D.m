%% Demonstrate/validate vesselAngle2D.m on synthetic test images
% Builds a handful of synthetic line images with known geometry (a
% straight line at an arbitrary angle, a straight line near the +-90 deg
% wrap boundary, a bent line with no true junction, and a Y-junction) and
% runs vesselAngle2D.m on each, plotting the input image, the vesselness
% map, the pixelwise tangent-angle map, and the resulting per-segment
% angle estimates. Useful as a quick regression check after touching
% vesselAngle2D.m or vesselnessAngle2D.m (../Jerman Enhancement Filter/).
%
% Set exportPNGs = true (or predefine it as true in the base workspace
% before running this script) to also export each figure as a PNG next
% to this file.

if ~exist('exportPNGs', 'var')
    exportPNGs = false;
end

sigmas = 1:0.5:3; % Vessel radii to search over [pixels], matches the calibers used in the synthetic images below

%% Case 1: straight line at an arbitrary angle (30 deg)
N = 80;
trueAngleDeg = 30;
I1 = syntheticLine(N, trueAngleDeg, 2.0);
runAndPlotCase(I1, sigmas, 'straight_line_30deg', trueAngleDeg, exportPNGs)

%% Case 2: straight line near the +-90 deg wrap boundary (-80 deg)
trueAngleDeg = -80;
I2 = syntheticLine(N, trueAngleDeg, 2.0);
runAndPlotCase(I2, sigmas, 'straight_line_neg80deg', trueAngleDeg, exportPNGs)

%% Case 3: bent line, no true 3-way junction -- should stay ONE segment
N3 = 80;
I3 = zeros(N3, N3);
for r = 10:70
    I3(r, 50) = 100; % vertical trunk
end
for c = 50:90
    r = round(70 + 0.6*(c - 50));
    if r <= N3
        I3(r, c) = 100; % diagonal continuation (kink, not a branch)
    end
end
I3 = imgaussfilt(I3, 0.7);
runAndPlotCase(I3, sigmas, 'bent_line_no_branch', [], exportPNGs)

%% Case 4: Y-junction (vertical trunk + 45 deg diagonal branch) -- should split into 3 segments
N4 = 120;
I4 = zeros(N4, N4);
for r = 10:100
    I4(r, 60) = 100; % vertical trunk + straight continuation
end
for t = 0:40
    r = 60 + t; c = 60 + t;
    if r <= N4 && c <= N4
        I4(r, c) = 100; % diagonal branch peeling off at the midpoint
    end
end
I4 = imgaussfilt(I4, 0.7);
runAndPlotCase(I4, 1:0.5:4, 'y_junction', [], exportPNGs)

findfigs

%% Helper functions

% Build a synthetic image of a single straight bright line at angleDeg
% (same convention as vesselnessAngle2D.m: 0 deg = along dim 1/rows,
% +-90 deg = along dim 2/columns), width ~2*halfWidthPix, on a dark
% background, through the center of an N x N image.
function I = syntheticLine(N, angleDeg, halfWidthPix)
    [X, Y] = meshgrid(1:N, 1:N); % X = dim2 (columns), Y = dim1 (rows)
    theta = angleDeg * pi/180;
    dirRow = cos(theta); dirCol = sin(theta);
    normalRow = -dirCol; normalCol = dirRow; % unit normal to the line
    cx = N/2; cy = N/2;
    dist = (Y - cy).*normalRow + (X - cx).*normalCol;
    I = double(abs(dist) < halfWidthPix) * 100;
    I = imgaussfilt(I, 0.5);
end

% Run vesselAngle2D on a synthetic image and plot: input image,
% vesselness map, pixelwise angle map, and the detected vessel segments
% colored by their per-segment angle (with measured vs. expected angle
% and straightness reported in the title).
function runAndPlotCase(I, sigmas, caseTitle, trueAngleDeg, exportPNGs)
    [vessels, angleMap, vesselness, ~, segLabel, ~] = vesselAngle2D(I, sigmas, [1;1], 1, true, 0.1, 5, 5);

    segAngleImg = nan(size(segLabel));
    for k = 1:numel(vessels)
        segAngleImg(vessels(k).PixelIdxList) = vessels(k).angleDeg;
    end

    fig = figure('Name', caseTitle, 'Position', [100 100 800 800]);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile(tl)
    imagesc(I); colormap(gca, gray); axis image; axis off; title('Input image')

    nexttile(tl)
    imagesc(vesselness); colormap(gca, gray); axis image; axis off; title('Vesselness')

    nexttile(tl)
    imagesc(angleMap, 'AlphaData', ~isnan(angleMap)); axis image; axis off; set(gca, 'Color', [0.15 0.15 0.15])
    colormap(gca, hsv); clim(gca, [-90, 90]); colorbar
    title('Pixelwise tangent angle [deg]')

    nexttile(tl)
    imagesc(segAngleImg, 'AlphaData', ~isnan(segAngleImg)); axis image; axis off; set(gca, 'Color', [0.15 0.15 0.15])
    colormap(gca, hsv); clim(gca, [-90, 90]); colorbar
    if isempty(trueAngleDeg)
        title(['Segments (' num2str(numel(vessels)) ')'])
    else
        title(sprintf('Segments (expected %.1f deg)', trueAngleDeg))
    end
    segStr = sprintf('seg %d: %.1f deg (straightness %.2f)\n', ...
        [1:numel(vessels); [vessels.angleDeg]; [vessels.straightness]]);
    xlabel(segStr, 'FontSize', 8, 'Color', 'k', 'Visible', 'on')

    sgtitle(tl, strrep(caseTitle, '_', ' '), 'Interpreter', 'none')

    if exportPNGs
        exportgraphics(fig, fullfile(fileparts(mfilename('fullpath')), [caseTitle '.png']));
    end
end
