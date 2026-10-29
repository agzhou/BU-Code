%% Description

% Plots the point spread function of a simulated point scatterer (RC15gV, one RF file) reconstructed with RCA_DAS_rcabeam,
% next to the Verasonics recon and the ffdas recon of the same data when they are in the folder (Half wl voxel spacing\IQ-*.mat,
% IQ Data - ffdas\IQ-*.mat). Slices through each volume's peak in dB, and line profiles through the peak with the -6 dB widths.

% Set beforehand (all optional):
%   psfPath    folder with params.mat, the RF file(s) and the other reconstructions
%   backend    'reference' (default, CPU port of rcabeam's kernel, a few minutes) or 'rcabeam' (GPU, R2026b)
%   cacheFile  .mat file to keep the reconstruction in, so replotting is instant
%   saveDir    folder to save the figures to as PNG
%   dBrange    display range (default 60 dB)

if exist('psfPath', 'var') ~= 1
    psfPath = 'G:\My Drive\Data\PSF Simulations\RC15gV PSF sim - 5 angles from -6 to 6 deg';
end
if exist('backend', 'var') ~= 1, backend = 'reference'; end
if exist('cacheFile', 'var') ~= 1, cacheFile = ''; end
if exist('saveDir', 'var') ~= 1, saveDir = ''; end
if exist('dBrange', 'var') ~= 1, dBrange = 60; end

S = load(fullfile(psfPath, 'params.mat'), 'P'); P = S.P;
if ~isfield(P, 'numFramesPerBuffer') && isfield(P, 'numFramesPerBuff'), P.numFramesPerBuffer = P.numFramesPerBuff; end
wl = P.Resource.Parameters.speedOfSound / (P.Trans.frequency * 1e6);
pix = 0.5 * wl; % the Half wl voxel spacing grid
rfFile = dir(fullfile(psfPath, 'RF-*.mat')); rfFile = rfFile(1).name;

%% Reconstruct (x, y, z) with the rcabeam port
if ~isempty(cacheFile) && exist(cacheFile, 'file')
    C = load(cacheFile); IQ_CR = C.IQ_CR; IQ_RC = C.IQ_RC; grid = C.grid;
    fprintf('Loaded the reconstruction from %s\n', cacheFile)
else
    R = load(fullfile(psfPath, rfFile), 'RcvData');
    r = unstackFrames(R.RcvData, P);
    clear R
    tic
    [IQ_CR, IQ_RC, grid] = RCA_DAS_rcabeam(r, P, pix, pix, pix, struct('backend', backend));
    fprintf('Reconstructed with the %s backend in %.0f s\n', backend, toc)
    if ~isempty(cacheFile)
        save(cacheFile, 'IQ_CR', 'IQ_RC', 'grid', '-v7.3')
    end
end
IQ_CR = IQ_CR(:, :, :, 1); IQ_RC = IQ_RC(:, :, :, 1);

% All volumes below are [y, x, z], like the Verasonics IQ
vols = struct('name', {}, 'iq', {}, 'x', {}, 'y', {}, 'z', {});
vols(1).name = 'rcabeam port (C-R + R-C)'; vols(1).iq = permute(IQ_CR + IQ_RC, [2, 1, 3]);
vols(1).x = grid.x; vols(1).y = grid.y; vols(1).z = grid.z;
vols(2).name = 'rcabeam port C-R'; vols(2).iq = permute(IQ_CR, [2, 1, 3]);
vols(2).x = grid.x; vols(2).y = grid.y; vols(2).z = grid.z;
vols(3).name = 'rcabeam port R-C'; vols(3).iq = permute(IQ_RC, [2, 1, 3]);
vols(3).x = grid.x; vols(3).y = grid.y; vols(3).z = grid.z;

vf = dir(fullfile(psfPath, 'Half wl voxel spacing', 'IQ-*.mat'));
if ~isempty(vf)
    V = load(fullfile(vf(1).folder, vf(1).name), 'IQ');
    PD = load(fullfile(vf(1).folder, 'PData.mat'), 'PData'); PD = PD.PData;
    k = numel(vols) + 1;
    vols(k).name = 'Verasonics recon'; vols(k).iq = V.IQ(:, :, :, 1);
    % Pixel positions from PData (the corner of the first pixel, in wavelengths); y direction is not checked
    vols(k).x = (PD.Origin(1) + (0 : size(V.IQ, 2) - 1) * PD.PDelta(1)) * wl;
    vols(k).y = (PD.Origin(2) - (0 : size(V.IQ, 1) - 1) * PD.PDelta(2)) * wl;
    vols(k).z = (PD.Origin(3) + (0 : size(V.IQ, 3) - 1) * PD.PDelta(3)) * wl;
end
ff = dir(fullfile(psfPath, 'IQ Data - ffdas', 'IQ-*.mat'));
if ~isempty(ff)
    F = load(fullfile(ff(1).folder, ff(1).name), 'IQ');
    FP = load(fullfile(ff(1).folder, 'ffdas_recon_params.mat'), 'grid');
    k = numel(vols) + 1;
    vols(k).name = 'ffdas recon'; vols(k).iq = F.IQ(:, :, :, 1);
    vols(k).x = FP.grid.x; vols(k).y = FP.grid.y; vols(k).z = FP.grid.z;
end

%% Peak, -6 dB widths, and the dB volumes
nv = numel(vols);
fprintf('\n%-26s %10s %9s %9s %9s %10s %10s %10s\n', 'volume', 'peak', 'x (um)', 'y (um)', 'z (um)', 'FWHM x', 'FWHM y', 'FWHM z')
for k = 1:nv
    A = abs(vols(k).iq);
    [vols(k).peak, i] = max(A(:));
    [iy, ix, iz] = ind2sub(size(A), i);
    vols(k).idx = [iy, ix, iz];
    vols(k).dB = 20 * log10(A / vols(k).peak);
    vols(k).prof = {squeeze(vols(k).dB(iy, :, iz)), squeeze(vols(k).dB(:, ix, iz))', squeeze(vols(k).dB(iy, ix, :))'};
    ax = {vols(k).x, vols(k).y, vols(k).z};
    for d = 1:3
        vols(k).fwhm(d) = fwhm(ax{d}, vols(k).prof{d}, [iy, ix, iz]);
    end
    fprintf('%-26s %10.3g %9.1f %9.1f %9.1f %9.0f um %9.0f um %9.0f um\n', vols(k).name, vols(k).peak, ...
        vols(k).x(ix) * 1e6, vols(k).y(iy) * 1e6, vols(k).z(iz) * 1e6, vols(k).fwhm * 1e6)
end
fprintf('(scatterer: x = 0, y = 0, z = %.0f um; wavelength %.0f um, voxel %.0f um)\n', P.Media.MP(3) * wl * 1e6, wl * 1e6, pix * 1e6)

%% Slices through each peak
show = [1, 4 : nv]; % OPW sum, Verasonics, ffdas
show = show(show <= nv);
fig1 = figure('Name', 'PSF slices', 'Position', [50, 50, 1250, 330 * numel(show)], 'Color', 'w');
tl = tiledlayout(numel(show), 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('PSF through the peak (dB re. peak), %s', fileparts(psfPath)), 'Interpreter', 'none')
for r = 1:numel(show)
    v = vols(show(r));
    iy = v.idx(1); ix = v.idx(2); iz = v.idx(3);
    sl = {squeeze(v.dB(iy, :, :))', squeeze(v.dB(:, ix, :))', squeeze(v.dB(:, :, iz))};
    a1 = {v.x, v.y, v.x}; a2 = {v.z, v.z, v.y};
    pk = {[v.x(ix), v.z(iz)], [v.y(iy), v.z(iz)], [v.x(ix), v.y(iy)]};
    lab = {{'x (mm)', 'z (mm)'}, {'y (mm)', 'z (mm)'}, {'x (mm)', 'y (mm)'}};
    for c = 1:3
        nexttile
        imagesc((a1{c} - pk{c}(1)) * 1e3, (a2{c} - pk{c}(2)) * 1e3, sl{c}, [-dBrange, 0])
        axis image; colormap(gca, gray); xlabel(lab{c}{1}); ylabel(lab{c}{2})
        xlim([-2, 2]); ylim([-2, 2]) % zoom on the PSF (mm from the peak)
        if c == 1, title(v.name, 'Interpreter', 'none'), end
    end
end
cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'dB re. peak of this volume';

%% Line profiles through the peaks
fig2 = figure('Name', 'PSF profiles', 'Position', [50, 50, 1250, 380], 'Color', 'w');
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Envelope profiles through the peak', 'Interpreter', 'none')
dn = {'x', 'y', 'z'}; col = lines(nv);
for d = 1:3
    nexttile; hold on
    for k = 1 : nv
        ax = {vols(k).x, vols(k).y, vols(k).z};
        ctr = ax{d}(vols(k).idx(c2i(d)));
        ls = '-'; if k == 2 || k == 3, ls = ':'; end
        plot((ax{d} - ctr) * 1e3, vols(k).prof{d}, ls, 'Color', col(k, :), 'LineWidth', 1.3, 'DisplayName', vols(k).name)
    end
    yline(-6, 'k--', 'HandleVisibility', 'off')
    hold off; grid on; xlim([-2, 2]); ylim([-dBrange, 3])
    xlabel(sprintf('%s from the peak (mm)', dn{d})); ylabel('dB re. peak')
    if d == 1, legend('Location', 'southoutside', 'Interpreter', 'none'), end
end

if ~isempty(saveDir)
    exportgraphics(fig1, fullfile(saveDir, 'PSF_slices_rcabeam.png'), 'Resolution', 150)
    exportgraphics(fig2, fullfile(saveDir, 'PSF_profiles_rcabeam.png'), 'Resolution', 150)
end


%% -6 dB width of a profile (m), by linear interpolation of the two crossings around the peak
function w = fwhm(ax, prof, ~)

    [~, ip] = max(prof);
    lo = NaN; hi = NaN;
    for i = ip : -1 : 2
        if prof(i-1) < -6
            lo = interp1(prof(i-1:i), ax(i-1:i), -6); break
        end
    end
    for i = ip : numel(prof) - 1
        if prof(i+1) < -6
            hi = interp1(prof(i:i+1), ax(i:i+1), -6); break
        end
    end
    w = abs(hi - lo);

end


% Index into [iy, ix, iz] for the x, y, z profile directions
function i = c2i(d)

    i = [2, 1, 3];
    i = i(d);

end
