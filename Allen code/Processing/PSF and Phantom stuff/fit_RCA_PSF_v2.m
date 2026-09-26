%% Description: fit the RCA-specific PSF model (sum of two "orthogonal" lateral Gaussians x an axial Gaussian) to simulated PSF data in 3D
clearvars

%% Add the Processing folder to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
AllenProcessingCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing\");
addpath(genpath(AllenProcessingCodePath))

%% Load simulation data and parameters (assuming Verasonics)
datapath = uigetdir('G:\', 'Select the data path');
datapath = [datapath, '\'];

if isfile([datapath, 'params.mat']) % params.mat is either next to PSF.mat or one folder up, depending on the dataset
    load([datapath, 'params.mat'])
else
    load([datapath, '..\params.mat'])
end
load([datapath, 'PData.mat'])
load([datapath, 'PSF.mat'])

%% Define coordinates
[x_mm, y_mm, z_mm] = getReconCoords3D(PData, P); % IQ is indexed as [x, y, z]

%% Plot the data
figure; imagesc(x_mm, z_mm, squeeze(max(abs(IQ), [], 2))'); title('xz MIP'); colormap gray; xlabel('x [mm]'); ylabel('z [mm]'); axis equal; axis tight
figure; imagesc(x_mm, y_mm, squeeze(max(abs(IQ), [], 3))'); title('xy MIP'); colormap gray; xlabel('x [mm]'); ylabel('y [mm]'); axis equal; axis tight

%% Get the position of the point scatterer - assuming there is only one scatterer
pos_wl = P.Media.MP(1, 1:3); % Position in wavelengths
pos_m = pos_wl .* P.wl; % [x, y, z] position in meters

% Known scatterer position = fixed center of the PSF
x0 = pos_m(1)*1e3; y0 = pos_m(2)*1e3; z0 = pos_m(3)*1e3; % [mm]

%% Fit window
% The model is fit to the envelope abs(IQ), so the sigmas are field (amplitude) sigmas, not intensity sigmas.
%
% The real RCA PSF has long, low-amplitude cross-artifact arms (~10% of the peak out to the edge of the volume in the
% 5-angle simulation, ~1% with 11 angles) and an axial spread that grows with lateral offset. The model can follow
% neither, and over the whole volume the arms can dominate the least-squares cost and pull the fit to a spurious
% solution (in the 5-angle case: sigma_wide of several mm, and the core amplitude underestimated ~5x). So the fit is
% restricted to a window around the scatterer, and the window is derived from the data in two stages:
%   1. Pilot fit in a small window seeded from the measured FWHM of the PSF core.
%   2. Final fit in a window with half-widths of window_sigma_mult * the pilot sigma (sigma_wide laterally,
%      sigma_axial axially). By 4 sigma the model has decayed to <0.03% of its peak (exp(-8)), so data further out
%      can only be arms/noise that the model can't represent, and would just bias sigma_wide upward.
% The window is deliberately not iterated to self-consistency: a larger window gives a larger sigma_wide, which gives
% a larger window, and beyond a few mm that runs away to the spurious solution described above. With 4 sigma the
% final fit didn't depend on the pilot window (tested with pilot windows of 1.5 - 5 x FWHM on the 5-angle data).
window_FWHM_mult = 3; % Pilot window half-width, in multiples of the measured FWHM
window_sigma_mult = 4; % Final window half-width, in multiples of the pilot sigma
min_halfwidth_vox = 5; % Floor on the window half-width [voxels]

% Measured FWHM of the PSF core, through the voxel nearest the scatterer
[~, ix0] = min(abs(x_mm - x0)); [~, iy0] = min(abs(y_mm - y0)); [~, iz0] = min(abs(z_mm - z0));
fwhm_lat_mm = mean([fwhm(x_mm, reshape(abs(IQ(:, iy0, iz0)), 1, [])), fwhm(y_mm, reshape(abs(IQ(ix0, :, iz0)), 1, []))]); % Lateral (mean of x and y)
fwhm_ax_mm = fwhm(z_mm, reshape(abs(IQ(ix0, iy0, :)), 1, [])); % Axial
dx = mean(diff(x_mm)); dz = mean(diff(z_mm)); % Voxel sizes [mm]

%% Define the model and the fit
PSF = @(x, y, z, x0, y0, z0, A, sigma_narrow, sigma_wide, sigma_axial) ...
        A.* ( exp(-(x - x0).^2 ./ (2*sigma_narrow^2) - (y - y0).^2 ./ (2*sigma_wide^2)) ...
        + exp(-(x - x0).^2 ./ (2*sigma_wide^2) - (y - y0).^2 ./ (2*sigma_narrow^2)) ) ...
        .* exp(-(z - z0).^2 ./ (2*sigma_axial^2));
% test = PSF(xr, yr, zr, x0, y0, z0, 1, 0.05, 0.25, 0.04); volumeViewer(test)

% Fit parameters: p = [A, sigma_narrow, sigma_wide/sigma_narrow, sigma_axial] (sigmas in mm)
%   The model is symmetric under sigma_narrow <-> sigma_wide, so fitting sigma_wide directly lets the optimizer
%   land on a mirrored "narrow > wide" solution. Fitting the ratio with a lower bound of 1 rules that out.
lb = [0,   0.005, 1,  0.005];
ub = [Inf, 1,     50, 1];

% Starting points [A, sigma_narrow, ratio, sigma_axial]. The first is the original guess (50 / 250 / 40 um); the
% rest are only there to check that the answer doesn't depend on the start.
p0_all = [0.5, 0.050, 5,  0.040;
          0.5, 0.100, 3,  0.030;
          0.5, 0.020, 20, 0.080;
          0.5, 0.150, 10, 0.050;
          0.3, 0.060, 30, 0.050];

opts = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-10, 'MaxIterations', 500, 'MaxFunctionEvaluations', 5000);

%% Fit: stage 1 = pilot window from the measured FWHM, stage 2 = final window from the pilot sigmas
stage_names = {'pilot', 'final'};
for stage = 1:2
    if stage == 1
        halfwidth_lat_mm = window_FWHM_mult * fwhm_lat_mm;
        halfwidth_ax_mm = window_FWHM_mult * fwhm_ax_mm;
    else
        halfwidth_lat_mm = window_sigma_mult * p_fit(2) * p_fit(3); % sigma_wide
        halfwidth_ax_mm = window_sigma_mult * p_fit(4); % sigma_axial
    end
    halfwidth_lat_mm = max(halfwidth_lat_mm, min_halfwidth_vox * dx);
    halfwidth_ax_mm = max(halfwidth_ax_mm, min_halfwidth_vox * dz);

    ix = find(abs(x_mm - x0) <= halfwidth_lat_mm);
    iy = find(abs(y_mm - y0) <= halfwidth_lat_mm);
    iz = find(abs(z_mm - z0) <= halfwidth_ax_mm);

    % Everything is fit in mm on peak-normalized data, so all fit parameters are O(0.01 - 1). Fitting the raw data
    % in meters (IQ ~ 1e10, sigmas ~ 1e-4) is badly scaled for finite-difference Jacobians.
    D = abs(IQ(ix, iy, iz)); % Data in the fit window, dimensions are [x, y, z] like IQ
    data_scale = max(D, [], 'all');
    D = D ./ data_scale;

    % Coordinate vectors oriented along their data dimension, so the model broadcasts to the [x, y, z] window
    xr = reshape(x_mm(ix), [], 1);
    yr = reshape(y_mm(iy), 1, []);
    zr = reshape(z_mm(iz), 1, 1, []);

    residual_fun = @(p) reshape(PSF(xr, yr, zr, x0, y0, z0, p(1), p(2), p(2)*p(3), p(4)) - D, [], 1);

    p_all = zeros(size(p0_all));
    cost_all = zeros(size(p0_all, 1), 1);
    for k = 1:size(p0_all, 1)
        [p_all(k, :), cost_all(k)] = lsqnonlin(residual_fun, p0_all(k, :), lb, ub, opts);
    end
    [cost_best, best] = min(cost_all);
    p_fit = p_all(best, :);
    n_agree = sum(abs(cost_all - cost_best) < 1e-6 * cost_best); % Number of start points that reached the best fit
    R2 = 1 - cost_best / sum((D - mean(D, 'all')).^2, 'all'); % In the fit window

    fprintf('%s fit: window +-%.2f mm lateral x +-%.2f mm axial (%d x %d x %d voxels), sigma_narrow / sigma_wide / sigma_axial = %.1f / %.1f / %.1f um, R^2 = %.3f, %d/%d start points agree\n', ...
        stage_names{stage}, halfwidth_lat_mm, halfwidth_ax_mm, numel(ix), numel(iy), numel(iz), ...
        p_fit(2)*1e3, p_fit(2)*p_fit(3)*1e3, p_fit(4)*1e3, R2, n_agree, size(p0_all, 1))

    if stage == 1
        p_pilot = p_fit;
    end
end

% The core sigmas shouldn't depend on the window. If they do, the window has started to include arms or noise.
if any(abs(p_fit([2 4]) ./ p_pilot([2 4]) - 1) > 0.2)
    warning('sigma_narrow or sigma_axial changed by more than 20%% between the pilot and final fits. The fit depends on the window, so check the plots.')
end

%% Results
A_fit = p_fit(1) * data_scale; % Amplitude in the original IQ units
sigma_field_mm = [p_fit(2), p_fit(2)*p_fit(3), p_fit(4)]; % [narrow, wide, axial]
sigma_field_um = sigma_field_mm .* 1e3;
sigma_intensity_um = sigma_field_um ./ sqrt(2); % Intensity sigmas, same convention as calcPSFShape3D

model_at = @(x, y, z) PSF(x, y, z, x0, y0, z0, p_fit(1), sigma_field_mm(1), sigma_field_mm(2), sigma_field_mm(3));

fprintf('Final: sigma_narrow / sigma_wide / sigma_axial = %.1f / %.1f / %.1f um (field), %.1f / %.1f / %.1f um (intensity)\n', ...
    sigma_field_um, sigma_intensity_um)

%% Plot the fit against the data (through the scatterer, full extent of the volume)
xPSF_data = abs(IQ(:, iy0, iz0)) ./ data_scale;
xPSF_fit = model_at(reshape(x_mm, [], 1), y_mm(iy0), z_mm(iz0));
zPSF_data = squeeze(abs(IQ(ix0, iy0, :))) ./ data_scale;
zPSF_fit = model_at(x_mm(ix0), y_mm(iy0), reshape(z_mm, [], 1));
xySlice_data = abs(IQ(:, :, iz0)) ./ data_scale;
xySlice_fit = model_at(reshape(x_mm, [], 1), reshape(y_mm, 1, []), z_mm(iz0));

figure('Name', 'RCA PSF fit'); tiledlayout(2, 3, 'TileSpacing', 'compact');
nexttile; plot(x_mm, xPSF_data, '.-'); hold on; plot(x_mm, xPSF_fit, 'LineWidth', 1.5); xline(x0 + halfwidth_lat_mm*[-1 1], ':'); hold off
xlim(x0 + 2*halfwidth_lat_mm*[-1 1]); xlabel('x [mm]'); title('x profile'); legend('Simulated', 'Fit', 'Fit window'); grid on

nexttile; semilogy(x_mm, xPSF_data + eps, '.-'); hold on; semilogy(x_mm, xPSF_fit + eps, 'LineWidth', 1.5); xline(x0 + halfwidth_lat_mm*[-1 1], ':'); hold off
ylim([1e-3 2]); xlabel('x [mm]'); title('x profile, log (arms extend past the fit window)'); grid on

nexttile; plot(z_mm, zPSF_data, '.-'); hold on; plot(z_mm, zPSF_fit, 'LineWidth', 1.5); xline(z0 + halfwidth_ax_mm*[-1 1], ':'); hold off
xlim(z0 + 2*halfwidth_ax_mm*[-1 1]); xlabel('z [mm]'); title('z profile'); grid on

slice_titles = {'xy slice: simulated', 'xy slice: fit', 'xy slice: simulated - fit'};
slice_data = {xySlice_data, xySlice_fit, xySlice_data - xySlice_fit};
slice_clims = {[0 1], [0 1], [-0.5 0.5]};
for k = 1:3
    ax = nexttile; imagesc(x_mm, y_mm, slice_data{k}'); axis image; clim(slice_clims{k}); title(slice_titles{k}); xlabel('x [mm]'); ylabel('y [mm]')
    rectangle(ax, 'Position', [x0 - halfwidth_lat_mm, y0 - halfwidth_lat_mm, 2*halfwidth_lat_mm*[1 1]], 'EdgeColor', 'w', 'LineStyle', ':') % Fit window
end
colorbar
