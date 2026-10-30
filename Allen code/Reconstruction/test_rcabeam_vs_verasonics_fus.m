%% Description

% Checks that RCA_DAS_rcabeam gives IQ data that leads to the same Power Doppler and g1 as the Verasonics reconstruction, on a
% real acquisition (RC15gV, undersampled BS67BW). The Verasonics recon here is the mid xz plane ('IQ - midplane', from
% RCA_VerasonicsRecon_loop_unstackedFrames_midXZ.m), so the same plane is reconstructed with rcabeam (RCA_DAS_rcabeam with
% yIndex, on the Verasonics PData grid) from the RF file of the same number. Both IQ sets then go through the processing of
% IQ2g1T_3D_loop.m (IQ2g1T_3D: SVD clutter filter, high pass filter, PDI, CDI, g1) and are compared:
%   1. the raw IQ, per voxel over time (complex correlation: the phase reference and the scale of the two recons differ by a
%      constant per voxel, which does not affect PDI or g1, so a correlation magnitude near 1 is a match)
%   2. the clutter filtered IQ (the blood signal, which is what PDI and g1 are made of), same measure
%   3. the PDI, CDI and g1 maps/curves
% Which PData row the Verasonics mid plane is (rows 40 and 41 are the same distance from y = 0) is found from the data.

% Set beforehand (all optional):
%   runPath      folder with params.mat, the RF files and 'IQ - midplane'
%   fileNum      RF/IQ file number (default 1)
%   nFrames      # of frames to compare (default 1000, all of a file). The CPU reference backend is slow, 100-200 is faster
%   backend      'rcabeam' (GPU, R2026b) or 'reference' (CPU port of the kernel, default when the GPU path isn't available)
%   fNumber      receive f-number of the rcabeam recon (default 0.731 = 1/(2*tan(0.6)), the Verasonics Recon.senscutoff of 0.6 rad)
%   svRange      [lower, upper] singular value cut for the clutter filter (default [30, nFrames], like IQ2g1T_3D_loop.m)
%   zRange       z voxels to process (default [20, 124])
%   cacheFile    .mat file to keep the rcabeam IQ in, so a rerun skips the reconstruction
%   saveDir      folder to save the figures to as PNG
%   pythonExe    python.exe with cupy and rcabeam, for backend 'rcabeam'

if exist('runPath', 'var') ~= 1
    runPath = 'H:\Ultrasound\10-08-2026 AZ07 air puff\run 2 - RC15gV';
end
if exist('fileNum', 'var') ~= 1, fileNum = 1; end
if exist('nFrames', 'var') ~= 1, nFrames = 1000; end
if exist('backend', 'var') ~= 1
    backend = 'reference';
    try %#ok<TRYNC>
        if ~isMATLABReleaseOlderThan("R2026b"), backend = 'rcabeam'; end
    end
end
if exist('fNumber', 'var') ~= 1, fNumber = 1 / (2 * tan(0.6)); end
if exist('svRange', 'var') ~= 1, svRange = [30, nFrames]; end
if exist('zRange', 'var') ~= 1, zRange = [20, 124]; end
if exist('cacheFile', 'var') ~= 1, cacheFile = ''; end
if exist('saveDir', 'var') ~= 1, saveDir = ''; end
if exist('pythonExe', 'var') ~= 1, pythonExe = ''; end

addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Processing'))) % IQ2g1T_3D and what it calls

S = load(fullfile(runPath, 'params.mat'), 'P'); P = S.P;
wl = P.Resource.Parameters.speedOfSound / (P.Trans.frequency * 1e6);
pix = [P.Trans.spacing, P.Trans.spacing, 0.5] * wl; % the Verasonics PDelta of the fUS runs, [x, y, z] in m
nS = P.Receive(1).endSample; nspf = nS * P.na * 2; % samples per frame
fprintf('%s, file %d, %d frames, sampleMode %s (demodFrequency %.4g MHz, Trans.frequency %.4g MHz), rcabeam backend "%s"\n', ...
    runPath, fileNum, nFrames, P.Receive(1).sampleMode, P.Receive(1).demodFrequency, P.Trans.frequency, backend)

%% Verasonics mid plane IQ: [x, z, frames]
vFile = fullfile(runPath, 'IQ - midplane', sprintf('IQ-%d-%d-%d-%d-1-%d.mat', round(P.maxAngle), P.na, P.frameRate, P.numFramesPerBuffer, fileNum));
V = load(vFile, 'IQ');
IQv = reshape(V.IQ(:, :, 1 : nFrames), size(V.IQ, 1), 1, size(V.IQ, 2), nFrames); % [x, 1, z, frames]
clear V
PD = load(fullfile(runPath, 'IQ - midplane', 'PData.mat')); PD = PD.PData;
midRow = round((PD.Origin(2) - PD.Region(2).Shape.oPAIntersect) / PD.PDelta(2)) + 1;
candRows = [midRow - 1, midRow];

%% rcabeam IQ of the same plane
if ~isempty(cacheFile) && exist(cacheFile, 'file') && ~isempty(who('-file', cacheFile, 'IQr'))
    C = load(cacheFile);
    if isequal(C.backend, backend) && C.nFrames == nFrames && C.fileNum == fileNum && isequal(C.fNumber, fNumber)
        IQr = C.IQr; row = C.row; rho = C.rho;
        fprintf('Loaded the rcabeam IQ from %s\n', cacheFile)
    end
end
if ~exist('IQr', 'var')
    fprintf('Loading the RF data...\n')
    m = matfile(fullfile(runPath, sprintf('RF-%d-%d-%d-%d-1-%d.mat', round(P.maxAngle), P.na, P.frameRate, P.numFramesPerBuffer, fileNum)));
    r = permute(reshape(m.RcvData(1 : nspf * nFrames, :), nspf, nFrames, []), [1, 3, 2]); % unstacked: [samples per frame, channels, frames]
    base = struct('backend', backend, 'fNumber', fNumber, 'pythonExe', pythonExe, 'sumOnly', true, 'frameChunk', 16);

    % Which of the two rows next to y = 0 is the Verasonics plane: reconstruct both for a few frames
    nProbe = min(nFrames, 40);
    o = base; o.yIndex = candRows;
    tic
    [probe, ~, grid] = RCA_DAS_rcabeam(r(:, :, 1 : nProbe), P, pix(2), pix(1), pix(3), o);
    fprintf('Probe (%d frames, 2 rows): %.0f s\n', nProbe, toc)
    rho = zeros(1, 2);
    for q = 1:2
        rho(q) = median(voxelCorr(squeeze(probe(:, q, :, :)), squeeze(IQv(:, 1, :, 1 : nProbe)), 0.1));
    end
    [~, q] = max(rho); row = candRows(q);
    fprintf('Raw IQ correlation (median over the brightest voxels) with PData row %d: %.3f, row %d: %.3f -> using row %d (y = %.1f um)\n', ...
        candRows(1), rho(1), candRows(2), rho(2), row, grid.y(q) * 1e6) % grid.y has only the two rows

    o = base; o.yIndex = row;
    tic
    [IQr4, ~, grid] = RCA_DAS_rcabeam(r, P, pix(2), pix(1), pix(3), o);
    fprintf('Reconstructed %d frames with the %s backend: %.0f s\n', nFrames, backend, toc)
    IQr = IQr4; clear IQr4 r
    if ~isempty(cacheFile)
        save(cacheFile, 'IQr', 'row', 'rho', 'backend', 'nFrames', 'fileNum', 'fNumber', '-v7.3')
    end
end
IQr = reshape(IQr, size(IQr, 1), 1, size(IQr, 3), nFrames); % [x, 1, z, frames]
RCA_DAS_rcabeam('clear');

%% 1. Raw IQ
zr = zRange(1) : zRange(2);
rawRho = voxelCorr(squeeze(IQr(:, 1, zr, :)), squeeze(IQv(:, 1, zr, :)), 0.1);
fprintf('\n1. Raw IQ, per voxel complex correlation over time (brightest 90%% of voxels): median %.3f, 10th percentile %.3f\n', ...
    median(rawRho), prctile(rawRho, 10))
% amplitude relation (the Verasonics IQ is scaled differently): ratio of the temporal RMS
rmsr = sqrt(mean(abs(IQr(:, 1, zr, :)) .^ 2, 4)); rmsv = sqrt(mean(abs(IQv(:, 1, zr, :)) .^ 2, 4));
fprintf('   Image amplitude, rcabeam / Verasonics: median %.3g, spread (10-90th percentile) %.3g to %.3g\n', ...
    median(rmsr(:) ./ rmsv(:)), prctile(rmsr(:) ./ rmsv(:), 10), prctile(rmsr(:) ./ rmsv(:), 90))
fprintf('   Correlation of the (log) mean image: %.3f\n', corr(log10(rmsr(:)), log10(rmsv(:))))

% Phase of the two recons relative to each other per voxel: the reference phase is not the Verasonics one, so a phase offset
% that depends on depth but not on x is expected, while one that is random from voxel to voxel would be a timing mismatch
dphi = angle(sum(squeeze(IQr(:, 1, zr, :)) .* conj(squeeze(IQv(:, 1, zr, :))), 3)); % [x, z]
resid = angle(exp(1j * (dphi - angle(mean(exp(1j * dphi), 1))))); % remove the depth dependent offset (circular mean over x)
fprintf('   Phase offset between the recons: std over x at a depth (rad) median %.2f; ', median(std(resid, 0, 1)))
fprintf('circular mean over x runs from %.2f to %.2f rad over the depths\n', min(angle(mean(exp(1j * dphi), 1))), max(angle(mean(exp(1j * dphi), 1))))

%% 2. Clutter filtered IQ, high pass filtered (what IQ2g1T_3D does before PDI and g1)
frameRate = P.frameRate; % IQ2g1T_3D_loop.m uses checkFrameRate(P) (actual rate from the RF time tags, needs a file dialog). Both IQ sets get the same one here
nTau = ceil(20e-3 * frameRate);
tau = (0 : nTau - 1)' ./ frameRate;
HPF.fc = 25; HPF.fs = P.frameRate; HPF.order = 4;
[HPF.b, HPF.a] = butter(HPF.order, HPF.fc / (HPF.fs / 2), 'high');
voxelRange = {1 : size(IQr, 1), 1, zr};

filt = cell(1, 2);
for k = 1:2
    if k == 1, X = IQr(voxelRange{:}, :); else, X = IQv(voxelRange{:}, :); end
    [CM, EVs, Vs] = getSVs2D(X);
    Xf = applySVs2D_accel(X, CM, EVs, Vs, svRange(1), svRange(2));
    filt{k} = filter(HPF.b, HPF.a, Xf, [], 4);
end
clear X Xf CM EVs Vs
rhoF = voxelCorr(squeeze(filt{1}), squeeze(filt{2}), 0);
engF = sum(abs(squeeze(filt{2})) .^ 2, 3);
w = engF(:) > prctile(engF(:), 90); % the voxels with the most blood signal
fprintf('\n2. Clutter filtered IQ (SVD %d-%d, HPF %g Hz), complex correlation per voxel: median over all voxels %.3f, over the top 10%% by energy %.3f\n', ...
    svRange, HPF.fc, median(rhoF(:)), median(rhoF(w)))

%% 3. PDI, CDI and g1 with the same function as IQ2g1T_3D_loop.m
out = cell(1, 2);
names = {'rcabeam', 'Verasonics'};
for k = 1:2
    if k == 1, X = IQr; else, X = IQv; end
    tic
    [out{k}.PDI, out{k}.CDI, out{k}.g1] = IQ2g1T_3D(X, P, voxelRange, svRange(1), svRange(2), HPF, nTau);
    fprintf('   IQ2g1T_3D %s: %.1f s\n', names{k}, toc)
end
clear X
pdi = cellfun(@(c) squeeze(c(:, 1, :)), {out{1}.PDI{3}, out{2}.PDI{3}}, 'UniformOutput', false);
cdi = cellfun(@(c) squeeze(c(:, 1, :)), {out{1}.CDI{3}, out{2}.CDI{3}}, 'UniformOutput', false);
% The PDI grows with depth (noise gain), so the "vessel" voxels are the ones standing out from the PDI of their own depth:
% PDI divided by the median PDI of its depth row, top 2% of the Verasonics image
flatV = pdi{2} ./ median(pdi{2}, 1);
top = flatV > prctile(flatV(:), 98);
fprintf('\n3. Power Doppler (all frequencies), rcabeam vs. Verasonics, %d x %d voxels:\n', size(pdi{1}))
fprintf('   correlation of PDI %.3f, of log10 PDI %.3f; PDI sum ratio (rcabeam / Verasonics) %.3g\n', ...
    corr(pdi{1}(:), pdi{2}(:)), corr(log10(pdi{1}(:)), log10(pdi{2}(:))), sum(pdi{1}(:)) / sum(pdi{2}(:)))
flatR = pdi{1} ./ median(pdi{1}, 1);
fprintf('   depth flattened PDI (PDI / median of its depth row): correlation %.3f, of log10 %.3f\n', ...
    corr(flatR(:), flatV(:)), corr(log10(flatR(:)), log10(flatV(:))))
fprintf('   vessel voxels (top 2%% of the flattened Verasonics PDI, %d voxels): flattened PDI rcabeam / Verasonics, median %.3f (1 = the same contrast to the background); correlation of PDI %.3f\n', ...
    nnz(top), median(flatR(top) ./ flatV(top)), corr(pdi{1}(top), pdi{2}(top)))
fprintf('   Color Doppler (all frequencies), correlation over the vessel voxels: %.3f\n', corr(cdi{1}(top), cdi{2}(top)))

% g1: {negative, positive, all frequencies}, each [x, 1, z, lag] (or [x, z, lag])
g = cell(1, 2);
for k = 1:2
    G = out{k}.g1{3};
    G = reshape(G, [], size(G, ndims(G))); % [voxels, lags]
    g{k} = G;
end
if size(g{1}, 2) ~= nTau, warning('g1 has %d lags, expected %d', size(g{1}, 2), nTau), end
topV = top(:);
a1 = abs(g{1}(topV, :)); a2 = abs(g{2}(topV, :));
fprintf('\n   g1 (all frequencies) over the vessel voxels, %d lags (%.1f ms):\n', size(g{1}, 2), tau(end) * 1e3)
fprintf('   mean |g1| curve RMS difference %.3g (|g1| goes from 1 to about %.2f over the lags)\n', ...
    sqrt(mean((mean(a1, 1) - mean(a2, 1)) .^ 2)), mean(a2(:, end)))
fprintf('   per voxel |g1| curve RMS difference: median %.3g, 90th percentile %.3g\n', ...
    median(sqrt(mean((a1 - a2) .^ 2, 2))), prctile(sqrt(mean((a1 - a2) .^ 2, 2)), 90))
lagc = min(size(a1, 2), 6);
fprintf('   correlation over voxels of |g1| at lag %d (%.2f ms): %.3f\n', lagc - 1, tau(lagc) * 1e3, corr(a1(:, lagc), a2(:, lagc)))
% decay: lag where |g1| first drops below 1/e
decay = @(A) arrayfun(@(i) find([A(i, :), -1] < exp(-1), 1), 1 : size(A, 1)) - 1;
d1 = decay(a1); d2 = decay(a2);
fprintf('   lag (frames) at which |g1| < 1/e: median rcabeam %g, Verasonics %g; correlation over voxels %.3f\n', ...
    median(d1), median(d2), corr(d1(:), d2(:)))

%% Figures
x_mm = (0 : size(IQr, 1) - 1) * pix(1) * 1e3; z_mm = (voxelRange{3} - 1) * pix(3) * 1e3 + P.startDepth * wl * 1e3;
fig1 = figure('Name', 'rcabeam vs Verasonics: images', 'Position', [50, 50, 1300, 700], 'Color', 'w');
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s, file %d, mid plane (row %d), rcabeam backend %s', ' ', fileNum, row, backend), 'Interpreter', 'none')
for k = 1:2
    nexttile; imagesc(x_mm, z_mm, 20 * log10(squeeze(sqrt(mean(abs((k == 1) * IQr(voxelRange{:}, :) + (k == 2) * IQv(voxelRange{:}, :)) .^ 2, 4)))' ...
        ./ max(rmsr(:) * (k == 1) + rmsv(:) * (k == 2))), [-60, 0]); axis image; colormap(gca, gray); colorbar
    title([names{k}, ' mean |IQ| (dB)']); xlabel('x (mm)'); ylabel('z (mm)')
end
nexttile; imagesc(x_mm, z_mm, rhoF', [0, 1]); axis image; colormap(gca, parula); colorbar
title('clutter filtered IQ, correlation per voxel'); xlabel('x (mm)'); ylabel('z (mm)')
for k = 1:2
    nexttile; imagesc(x_mm, z_mm, (pdi{k} ./ max(pdi{k}(:)))' .^ 0.5); axis image; colormap(gca, hot); colorbar
    title([names{k}, ' PDI (sqrt, normalized)']); xlabel('x (mm)'); ylabel('z (mm)')
end
nexttile; loglog(pdi{2}(:) ./ mean(pdi{2}(:)), pdi{1}(:) ./ mean(pdi{1}(:)), '.', 'MarkerSize', 3); hold on
plot(xlim, xlim, 'k--'); grid on; xlabel('Verasonics PDI / mean'); ylabel('rcabeam PDI / mean'); title('PDI per voxel')

fig2 = figure('Name', 'rcabeam vs Verasonics: g1', 'Position', [50, 50, 1100, 420], 'Color', 'w');
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on
plot(tau * 1e3, mean(a2, 1), 'LineWidth', 1.6); plot(tau * 1e3, mean(a1, 1), '--', 'LineWidth', 1.6)
grid on; xlabel('lag (ms)'); ylabel('mean |g1|'); legend('Verasonics', 'rcabeam'); title('vessel voxels (top 2% flattened PDI)')
nexttile; plot(a2(:, lagc), a1(:, lagc), '.'); hold on; plot([0, 1], [0, 1], 'k--'); grid on; axis square
xlabel(sprintf('Verasonics |g1| at lag %d', lagc - 1)); ylabel('rcabeam |g1|'); title('per voxel')

if ~isempty(saveDir)
    exportgraphics(fig1, fullfile(saveDir, 'rcabeam_vs_Verasonics_images.png'), 'Resolution', 150)
    exportgraphics(fig2, fullfile(saveDir, 'rcabeam_vs_Verasonics_g1.png'), 'Resolution', 150)
end


%% Per voxel |complex correlation| over time between a and b, [x, z, frames]. Only the voxels above the given quantile of the
% energy of b (0: all). Without the quantile output zero-energy voxels are NaN
function rho = voxelCorr(a, b, q)

    num = abs(sum(a .* conj(b), 3));
    den = sqrt(sum(abs(a) .^ 2, 3) .* sum(abs(b) .^ 2, 3));
    rho = num ./ den;
    if q > 0
        e = sum(abs(b) .^ 2, 3);
        rho = rho(e > quantile(e(:), q));
    end

end
