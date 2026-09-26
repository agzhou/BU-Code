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

load([datapath, '..\params.mat'])
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
% The real RCA PSF has long, low-amplitude cross-artifact arms (~10% of the peak, out to the edge of the volume)
% and an axial spread that grows with lateral offset. The model can follow neither, and over the whole volume the
% arms dominate the least-squares cost and can pull the fit to a spurious solution (sigma_wide of several mm, and
% the core amplitude underestimated ~5x). Fitting only near the scatterer avoids that.
fit_halfwidth_lat_mm = 1; % lateral (x and y) half-width of the fit window [mm] (~3.5*sigma_wide)
fit_halfwidth_ax_mm = 0.4; % axial (z) half-width of the fit window [mm] (~99% of the PSF energy)

ix = find(abs(x_mm - x0) <= fit_halfwidth_lat_mm);
iy = find(abs(y_mm - y0) <= fit_halfwidth_lat_mm);
iz = find(abs(z_mm - z0) <= fit_halfwidth_ax_mm);

% Everything is fit in mm on peak-normalized data, so all fit parameters are O(0.01 - 1). Fitting the raw data in
% meters (IQ ~ 1e10, sigmas ~ 1e-4) is badly scaled for finite-difference Jacobians.
D = abs(IQ(ix, iy, iz)); % Data in the fit window, dimensions are [x, y, z] like IQ
data_scale = max(D, [], 'all');
D = D ./ data_scale;

% Coordinate vectors oriented along their data dimension, so the model broadcasts to the [x, y, z] window
xr = reshape(x_mm(ix), [], 1);
yr = reshape(y_mm(iy), 1, []);
zr = reshape(z_mm(iz), 1, 1, []);

%% Define the model and the fit
PSF = @(x, y, z, x0, y0, z0, A, sigma_narrow, sigma_wide, sigma_axial) ...
        A.* ( exp(-(x - x0).^2 ./ (2*sigma_narrow^2) - (y - y0).^2 ./ (2*sigma_wide^2)) ...
        + exp(-(x - x0).^2 ./ (2*sigma_wide^2) - (y - y0).^2 ./ (2*sigma_narrow^2)) ) ...
        .* exp(-(z - z0).^2 ./ (2*sigma_axial^2));
% test = PSF(xr, yr, zr, x0, y0, z0, 1, 0.05, 0.25, 0.04); volumeViewer(test)

% Fit parameters: p = [A, sigma_narrow, sigma_wide/sigma_narrow, sigma_axial] (sigmas in mm)
%   The model is symmetric under sigma_narrow <-> sigma_wide, so fitting sigma_wide directly lets the optimizer
%   land on a mirrored "narrow > wide" solution. Fitting the ratio with a lower bound of 1 rules that out.
residual_fun = @(p) reshape(PSF(xr, yr, zr, x0, y0, z0, p(1), p(2), p(2)*p(3), p(4)) - D, [], 1);
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
p_all = zeros(size(p0_all));
cost_all = zeros(size(p0_all, 1), 1);
for k = 1:size(p0_all, 1)
    [p_all(k, :), cost_all(k)] = lsqnonlin(residual_fun, p0_all(k, :), lb, ub, opts);
end

[~, best] = min(cost_all);
p_fit = p_all(best, :);

fprintf('Fit from each start point (sigma in um):\n')
fprintf('  start %d: A = %.3f, sigma_narrow = %6.1f, sigma_wide = %7.1f, sigma_axial = %5.1f, cost = %.5f\n', ...
    [1:size(p0_all, 1); p_all(:, 1)'; p_all(:, 2)'*1e3; p_all(:, 2)'.*p_all(:, 3)'*1e3; p_all(:, 4)'*1e3; cost_all'])

%% Results
A_fit = p_fit(1) * data_scale; % Amplitude in the original IQ units
sigma_field_mm = [p_fit(2), p_fit(2)*p_fit(3), p_fit(4)]; % [narrow, wide, axial]
sigma_field_um = sigma_field_mm .* 1e3;
sigma_intensity_um = sigma_field_um ./ sqrt(2); % Intensity sigmas, same convention as calcPSFShape3D

model_at = @(x, y, z) PSF(x, y, z, x0, y0, z0, p_fit(1), sigma_field_mm(1), sigma_field_mm(2), sigma_field_mm(3));

D_fit = model_at(xr, yr, zr);
R2 = 1 - sum((D - D_fit).^2, 'all') / sum((D - mean(D, 'all')).^2, 'all'); % In the fit window
fprintf('Best fit (start %d): sigma_narrow / sigma_wide / sigma_axial = %.1f / %.1f / %.1f um (field), R^2 = %.3f in the fit window\n', ...
    best, sigma_field_um, R2)

%% Plot the fit against the data (through the scatterer, full extent of the volume)
[~, ix0] = min(abs(x_mm - x0)); [~, iy0] = min(abs(y_mm - y0)); [~, iz0] = min(abs(z_mm - z0)); % Voxel nearest the scatterer

xPSF_data = abs(IQ(:, iy0, iz0)) ./ data_scale;
xPSF_fit = model_at(reshape(x_mm, [], 1), y_mm(iy0), z_mm(iz0));
zPSF_data = squeeze(abs(IQ(ix0, iy0, :))) ./ data_scale;
zPSF_fit = model_at(x_mm(ix0), y_mm(iy0), reshape(z_mm, [], 1));
xySlice_data = abs(IQ(:, :, iz0)) ./ data_scale;
xySlice_fit = model_at(reshape(x_mm, [], 1), reshape(y_mm, 1, []), z_mm(iz0));

figure('Name', 'RCA PSF fit'); tiledlayout(2, 3, 'TileSpacing', 'compact');
nexttile; plot(x_mm, xPSF_data, '.-'); hold on; plot(x_mm, xPSF_fit, 'LineWidth', 1.5); xline(x0 + fit_halfwidth_lat_mm*[-1 1], ':'); hold off
xlim(x0 + 2*fit_halfwidth_lat_mm*[-1 1]); xlabel('x [mm]'); title('x profile'); legend('Simulated', 'Fit', 'Fit window'); grid on

nexttile; semilogy(x_mm, xPSF_data + eps, '.-'); hold on; semilogy(x_mm, xPSF_fit + eps, 'LineWidth', 1.5); xline(x0 + fit_halfwidth_lat_mm*[-1 1], ':'); hold off
ylim([1e-3 2]); xlabel('x [mm]'); title('x profile, log (arms extend past the fit window)'); grid on

nexttile; plot(z_mm, zPSF_data, '.-'); hold on; plot(z_mm, zPSF_fit, 'LineWidth', 1.5); hold off
xlim(z0 + fit_halfwidth_ax_mm*[-1 1]); xlabel('z [mm]'); title('z profile'); grid on

slice_titles = {'xy slice: simulated', 'xy slice: fit', 'xy slice: simulated - fit'};
slice_data = {xySlice_data, xySlice_fit, xySlice_data - xySlice_fit};
slice_clims = {[0 1], [0 1], [-0.5 0.5]};
for k = 1:3
    ax = nexttile; imagesc(x_mm, y_mm, slice_data{k}'); axis image; clim(slice_clims{k}); title(slice_titles{k}); xlabel('x [mm]'); ylabel('y [mm]')
    rectangle(ax, 'Position', [x0 - fit_halfwidth_lat_mm, y0 - fit_halfwidth_lat_mm, 2*fit_halfwidth_lat_mm*[1 1]], 'EdgeColor', 'w', 'LineStyle', ':') % Fit window
end
colorbar
