%% Description

% Checks RCA_DAS_ffdas against a simulated point scatterer. The data is built directly from the RCA_DAS delay model
% (independent of ffdas), so this checks the things in the ffdas interface that the docs don't pin down:
%   - sample index base and units of the offsets (peak lands at the scatterer, within a pixel)
%   - sign of the phase rotation (compare the peak amplitude for wavenumSign +1 vs. -1; the wrong one is much lower)
%   - in-plane receive distance (exactRx true vs. false, scatterer placed off-center)
%   - out of range samples don't crash (the grid reaches beyond the acquired samples for the steered angles)
% Needs a CUDA GPU, the Parallel Computing Toolbox, and ffdas on the MATLAB path.

c = 1540; fc = 15.625e6; wl = c / fc; spw = 4; fs = fc * spw; dps = c / fs;
N = 128; pitch = 0.1e-3; nf = 2;
startDepth = 10; endDepth = 60; twPeak = 1; % wavelengths
% Number of angles and max angle (deg) can be set beforehand with numAngles and maxAngleDeg. The defaults have enough
% angles that the grating lobes (lambda/sin(dtheta)) fall outside the region
if exist('numAngles', 'var'), na = numAngles; else, na = 11; end
if exist('maxAngleDeg', 'var'), maxAng = maxAngleDeg; else, maxAng = 15; end
angles = linspace(-maxAng, maxAng, na) * pi/180;

xe = ((0:N-1) - (N-1)/2) * pitch; ye = xe;
x_end = max(xe); y_end = max(ye);

% Minimal P with the fields RCA_DAS_ffdas reads
P.Trans.frequency = fc / 1e6;
P.Trans.spacingMm = pitch * 1e3;
P.Trans.lensCorrection = 0;
P.Trans.Connector = (1 : 2*N)';
P.Trans.ElementPos = zeros(2*N, 4);
P.Trans.ElementPos(1:N, 1) = xe / wl;
P.Trans.ElementPos(N+1 : 2*N, 2) = ye / wl;
P.na = na; P.numElements = N; P.angles = angles;
P.startDepth = startDepth; P.endDepth = endDepth;
P.TW.peak = twPeak;
P.Resource.Parameters.speedOfSound = c;
P.numFramesPerBuffer = nf;

nS = 2 * (endDepth - startDepth) * spw; % samples per acquisition
for a = 1 : 2*na
    P.Receive(a).startSample = (a-1) * nS + 1;
    P.Receive(a).endSample = a * nS;
    P.Receive(a).samplesPerWave = spw;
end

%% Simulate a point scatterer (off-center, so the 3D vs. in-plane receive distance matters)
scat = [3.0e-3, -2.5e-3, 3.5e-3]; % x, y, z in m
j = (0 : nS-1)';
sigma = 1.5 * spw;
pulse = @(n0) exp(-0.5 * ((j - n0) / sigma).^2) .* cos(2*pi * (j - n0) / spw); % real RF, peak at (fractional) sample n0

rf = zeros(nS * 2*na, 2*N, nf);
for a = 1:na
    % Group 1: TX columns steered in x-z, RX rows (positioned along y)
    th = angles(a);
    tx = scat(3) * cos(th) + (scat(1) + sign(th) * x_end) * sin(th);
    rx = sqrt(scat(3)^2 + (scat(2) - ye).^2);
    n0 = (tx + rx + wl * twPeak - 2 * startDepth * wl) / dps; % 0-based sample index of the echo peak
    for ch = 1:N
        rf((a-1)*nS + (1:nS), N + ch, :) = repmat(pulse(n0(ch)), 1, 1, nf);
    end

    % Group 2: TX rows steered in y-z, RX columns (positioned along x). The angle order is flipped
    th = angles(na + 1 - a);
    tx = scat(3) * cos(th) + (sign(th) * y_end - scat(2)) * sin(th);
    rx = sqrt(scat(3)^2 + (scat(1) - xe).^2);
    n0 = (tx + rx + wl * twPeak - 2 * startDepth * wl) / dps;
    for ch = 1:N
        rf((na + a - 1)*nS + (1:nS), ch, :) = repmat(pulse(n0(ch)), 1, 1, nf);
    end
end

%% Reconstruct and locate the peak
cfgs = {'exactRx, wavenum +', struct('exactRx', true, 'wavenumSign', 1);
        'exactRx, wavenum -', struct('exactRx', true, 'wavenumSign', -1);
        '3D Rx distance    ', struct('exactRx', false, 'wavenumSign', 1)};
% cfgs = {'exactRx, wavenum -', struct('exactRx', true, 'wavenumSign', -1)};
if exist('useGPU', 'var') % can be set to false beforehand to run with a CPU stand-in for ffdas.das
    for k = 1 : size(cfgs, 1)
        cfgs{k, 2}.useGPU = useGPU;
    end
end

fprintf('Scatterer at x = %.0f, y = %.0f, z = %.0f um\n', scat * 1e6)
fprintf('%-20s %-6s %10s %8s %8s %8s\n', 'config', 'vol', 'peak', 'dx um', 'dy um', 'dz um')
if ~exist('pix', 'var') % pixel size, can be set beforehand for a coarser (faster) grid
    pix = wl / 2;
end
peaks = zeros(1, size(cfgs, 1));
posOK = true;
mags = cell(size(cfgs, 1), 2); % envelope of the first frame, for plotting
for k = 1 : size(cfgs, 1)
    [IQ_CR, IQ_RC, grid] = RCA_DAS_ffdas(rf, P, pix, pix, pix, cfgs{k, 2});
    vols = {IQ_CR, IQ_RC};
    names = {'C-R', 'R-C'};
    for v = 1:2
        mags{k, v} = abs(vols{v}(:, :, :, 1));
        [m, ind] = max(reshape(abs(vols{v}(:, :, :, 1)), [], 1));
        [ix, iy, iz] = ind2sub(size(vols{v}, 1:3), ind);
        err = [grid.x(ix) - scat(1), grid.y(iy) - scat(2), grid.z(iz) - scat(3)] * 1e6;
        fprintf('%-20s %-6s %10.3g %8.1f %8.1f %8.1f\n', cfgs{k, 1}, names{v}, m, err)
        if k == 1
            peaks(1) = max(peaks(1), m);
            posOK = posOK && all(abs(err) <= pix * 1e6);
        elseif k == 2
            peaks(2) = max(peaks(2), m);
        end
    end
end

fprintf('\nPosition within one pixel (%.0f um) with exactRx: %s\n', pix * 1e6, passFail(posOK))
fprintf('Peak with wavenum + vs. - (should be much higher for +): %.3g vs. %.3g -> %s\n', peaks(1), peaks(2), passFail(peaks(1) > 2 * peaks(2)))
fprintf('If the peak is a constant depth off (about %.0f um per sample), try Opts.sampleShift = +/-1.\n', dps / 2 * 1e6)

%% Plot the reconstructions against the true scatterer location
% Max intensity projections in dB, red circle = true location. Set saveFigDir beforehand to also save PNGs
for k = 1 : size(cfgs, 1)
    fig = figure('Name', strtrim(cfgs{k, 1}), 'Position', [100, 100, 1100, 650]);
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, strtrim(cfgs{k, 1}))
    for v = 1:2
        for p = 1:3
            nexttile
            plotProjection(mags{k, v}, grid, scat, p, names{v})
        end
    end
    cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'dB re. peak of this volume';
    if exist('saveFigDir', 'var')
        exportgraphics(fig, fullfile(saveFigDir, sprintf('test_RCA_DAS_ffdas_cfg%d.png', k)), 'Resolution', 150)
    end
end

function plotProjection(mag, grid, scat, p, name)
    dB = 20 * log10(mag / max(mag(:)));
    switch p
        case 1 % x-z, max over y
            im = squeeze(max(dB, [], 2))'; ax1 = grid.x; ax2 = grid.z; l = {'x (mm)', 'z (mm)'}; tp = scat([1, 3]);
        case 2 % y-z, max over x
            im = squeeze(max(dB, [], 1))'; ax1 = grid.y; ax2 = grid.z; l = {'y (mm)', 'z (mm)'}; tp = scat([2, 3]);
        case 3 % x-y, max over z
            im = squeeze(max(dB, [], 3))'; ax1 = grid.x; ax2 = grid.y; l = {'x (mm)', 'y (mm)'}; tp = scat([1, 2]);
    end
    imagesc(ax1 * 1e3, ax2 * 1e3, im, [-40, 0])
    hold on
    plot(tp(1) * 1e3, tp(2) * 1e3, 'ro', 'MarkerSize', 14, 'LineWidth', 1.5)
    hold off
    axis image
    xlabel(l{1}); ylabel(l{2}); title(name)
    colormap(gca, gray)
end

function s = passFail(tf)
    if tf
        s = 'PASS';
    else
        s = 'FAIL';
    end
end
