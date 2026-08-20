%% Demonstrate/validate vesselAngle3D.m on synthetic test volumes
% Builds a handful of synthetic tube volumes with known geometry (a
% straight tube at an arbitrary 3D direction, a straight tube in a
% different octant to sanity check the sign-canonicalization convention,
% a bent tube with no true junction, and a Y-junction with a genuinely
% 3D-diverging branch) and runs vesselAngle3D.m on each, plotting a MIP
% of the input volume, a MIP of the vesselness volume, a 3D scatter of
% the detected skeleton segments, and a text summary of each segment's
% recovered azimuth/elevation/straightness. 3D counterpart of
% example_vesselAngle2D.m; see that file for the 2D version.
%
% Set exportPNGs = true (or predefine it as true in the base workspace
% before running this script) to also export each figure as a PNG next
% to this file.

if ~exist('exportPNGs', 'var')
    exportPNGs = false;
end

sigmas = 1:0.5:3; % Vessel radii to search over [voxels], matches the calibers used in the synthetic volumes below
N = 50; % Volume is N x N x N

%% Case 1: straight tube at an arbitrary 3D direction (az = 30 deg, el = 20 deg)
trueAz = 30; trueEl = 20;
I1 = syntheticTube3D(N, trueAz, trueEl, 2.0);
runAndPlotCase3D(I1, sigmas, 'straight_tube_az30_el20', trueAz, trueEl, exportPNGs)

%% Case 2: straight tube in a different octant (az = -100 deg, el = -35 deg)
trueAz = -100; trueEl = -35;
I2 = syntheticTube3D(N, trueAz, trueEl, 2.0);
runAndPlotCase3D(I2, sigmas, 'straight_tube_az-100_el-35', trueAz, trueEl, exportPNGs)

%% Case 3: bent tube, no true 3-way junction -- should stay ONE segment
I3 = zeros(N, N, N);
I3 = addTubeSegment3D(I3, [10 25 25], [25 25 25], 2.0, 100); % arm 1, along dim1
I3 = addTubeSegment3D(I3, [25 25 25], [42 38 35], 2.0, 100); % arm 2, kinks off diagonally
I3 = imgaussian3(I3, 0.7);
runAndPlotCase3D(I3, sigmas, 'bent_tube_no_branch', [], [], exportPNGs)

%% Case 4: Y-junction (trunk + straight continuation + a diverging 3D branch)
I4 = zeros(N, N, N);
I4 = addTubeSegment3D(I4, [8 25 25], [26 25 25], 2.0, 100); % trunk
I4 = addTubeSegment3D(I4, [26 25 25], [44 25 25], 2.0, 100); % straight continuation
I4 = addTubeSegment3D(I4, [26 25 25], [44 40 38], 2.0, 100); % branch, diverges in dim2 AND dim3
I4 = imgaussian3(I4, 0.7);
runAndPlotCase3D(I4, sigmas, 'y_junction_3d', [], [], exportPNGs)

findfigs

%% Helper functions

% Build a synthetic volume of a single straight bright tube at
% (azimuthDeg, elevationDeg) -- same convention as vesselAngle3D.m's
% .azimuthDeg/.elevationDeg: direction = [cosd(el)*cosd(az),
% cosd(el)*sind(az), sind(el)] -- through the center of an N x N x N
% volume.
function I = syntheticTube3D(N, azimuthDeg, elevationDeg, halfWidthVox)
    c = [N N N]/2;
    dir = [cosd(elevationDeg)*cosd(azimuthDeg), cosd(elevationDeg)*sind(azimuthDeg), sind(elevationDeg)];
    halfLen = 0.4*N;
    A = c - halfLen*dir;
    B = c + halfLen*dir;
    I = addTubeSegment3D(zeros(N,N,N), A, B, halfWidthVox, 100);
    I = imgaussian3(I, 0.5);
end

% Add a bright tube of the given half-width along the finite segment A-B
% to vol (taking the max with whatever's already there), using the
% vectorized point-to-segment distance over the whole grid.
function vol = addTubeSegment3D(vol, A, B, halfWidthVox, brightness)
    [N1, N2, N3] = size(vol);
    [X1, X2, X3] = ndgrid(1:N1, 1:N2, 1:N3);
    AB = B - A;
    ABsq = dot(AB, AB);
    t = ((X1-A(1))*AB(1) + (X2-A(2))*AB(2) + (X3-A(3))*AB(3)) / ABsq;
    t = min(max(t, 0), 1);
    cx = A(1) + t*AB(1); cy = A(2) + t*AB(2); cz = A(3) + t*AB(3);
    d = sqrt((X1-cx).^2 + (X2-cy).^2 + (X3-cz).^2);
    vol = max(vol, double(d < halfWidthVox) * brightness);
end

% Minimal separable 3D Gaussian smoothing (isotropic sigma, pixel units)
% for lightly blurring the synthetic binary tubes -- same purpose as the
% imgaussfilt(I,0.5)/(0.7) calls in example_vesselAngle2D.m, just a 3D
% equivalent using only base image-processing primitives.
function I = imgaussian3(I, sigma)
    r = ceil(3*sigma);
    x = -r:r;
    h = exp(-(x.^2)/(2*sigma^2)); h = h/sum(h);
    I = imfilter(imfilter(imfilter(I, reshape(h,[],1,1), 'same', 'replicate'), reshape(h,1,[],1), 'same', 'replicate'), reshape(h,1,1,[]), 'same', 'replicate');
end

% Run vesselAngle3D on a synthetic volume and plot: MIP of the input
% volume, MIP of the vesselness volume, a 3D scatter of the detected
% skeleton segments (one color per segment), and a text summary of each
% segment's recovered azimuth/elevation/straightness (with the expected
% values, for the single-tube cases).
function runAndPlotCase3D(I, sigmas, caseTitle, trueAz, trueEl, exportPNGs)
    [vessels, ~, ~, ~, vesselness, ~, ~, ~] = vesselAngle3D(I, sigmas, [1;1;1], 1, true, 0.1, 5, 5);

    fig = figure('Name', caseTitle, 'Position', [100 100 1000 850]);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile(tl)
    imagesc(max(I, [], 3)); colormap(gca, gray); axis image; axis off; title('Input volume (MIP over dim3)')

    nexttile(tl)
    imagesc(max(vesselness, [], 3)); colormap(gca, gray); axis image; axis off; title('Vesselness (MIP over dim3)')

    nexttile(tl)
    hold on
    colors = lines(max(numel(vessels), 1));
    for k = 1:numel(vessels)
        c = vessels(k).coordsDim123;
        plot3(c(:,1), c(:,2), c(:,3), '.', 'Color', colors(k,:), 'MarkerSize', 16);
    end
    hold off
    axis equal; grid on; view(35, 20)
    xlabel('dim1'); ylabel('dim2'); zlabel('dim3')
    title(['Segments (' num2str(numel(vessels)) ')'])

    nexttile(tl)
    axis off
    lines_ = {['\bf' strrep(caseTitle, '_', ' ')]};
    if ~isempty(trueAz)
        lines_{end+1} = sprintf('expected  az %.1f  el %.1f', trueAz, trueEl);
    end
    lines_{end+1} = '';
    for k = 1:numel(vessels)
        lines_{end+1} = sprintf('seg %d  az %6.1f  el %6.1f  str %.2f', ...
            k, vessels(k).azimuthDeg, vessels(k).elevationDeg, vessels(k).straightness); %#ok<AGROW>
    end
    text(0.04, 0.95, lines_, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontName', 'FixedWidth')

    sgtitle(tl, strrep(caseTitle, '_', ' '), 'Interpreter', 'none')

    if exportPNGs
        exportgraphics(fig, fullfile(fileparts(mfilename('fullpath')), [caseTitle '.png']));
    end
end
