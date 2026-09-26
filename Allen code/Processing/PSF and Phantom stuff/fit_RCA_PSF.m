%% Description: semi-automatically calculate the ultrasound PSF shape in 3D, taking in simulation data
clearvars

%% Add the Processing folder to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
AllenProcessingCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing\");
addpath(genpath(AllenProcessingCodePath))

%% Load simulation data and parameters (assuming Verasonics)
% load('G:\My Drive\Data\PSF Simulations\RC15gV PSF sim - 11 angles from -5 to 5 deg\params.mat')
% load('G:\My Drive\Data\PSF Simulations\RC15gV PSF sim - 11 angles from -5 to 5 deg\PData.mat')
% load('G:\My Drive\Data\PSF Simulations\RC15gV PSF sim - 11 angles from -5 to 5 deg\PSF.mat')

datapath = uigetdir('G:\', 'Select the data path');
datapath = [datapath, '\'];

load([datapath, '..\params.mat'])
load([datapath, 'PData.mat'])
load([datapath, 'PSF.mat'])

%% Define coordinates
[x_mm, y_mm, z_mm] = getReconCoords3D(PData, P);

%% Plot the data
figure; imagesc(x_mm, z_mm, squeeze(max(abs(IQ), [], 1))'); title('xz MIP'); colormap gray; xlabel('x [mm]'); ylabel('z [mm]'); axis equal; axis tight
figure; imagesc(x_mm, y_mm, squeeze(max(abs(IQ), [], 3))'); title('xy MIP'); colormap gray; xlabel('x [mm]'); ylabel('y [mm]'); axis equal; axis tight

%% Get the position of the point scatterer - assuming there is only one scatterer
pos_wl = P.Media.MP(1, 1:3); % Position in wavelengths
pos_m = pos_wl .* P.wl; % [x, y, z] position in meters

%% Get the position of the point scatterer in terms of the reconstructed data coordinates

% This abs(PData.Origin) line is kind of a hack.. fix later if needed
pos_rc = round((pos_wl + abs(PData.Origin)) ./ PData.PDelta + 1); % Position in recon coords (add 1 because of Matlab indexing)


%% Fitting the RCA-specific PSF model to the data

x0 = pos_m(1); y0 = pos_m(2); z0 = pos_m(3);
sigma_narrow_test = 50e-6;
sigma_wide_test = 250e-6;
sigma_axial_test = 40e-6;
% x = linspace(x0 - sigma_wide*5, x0 + sigma_wide*5, 100);
% y = linspace(y0 - sigma_wide*5, y0 + sigma_wide*5, 100);
% z = linspace(z0 - sigma_wide*5, z0 + sigma_wide*5, 100);

[X, Y, Z] = meshgrid(x_mm./1e3, y_mm./1e3, z_mm./1e3);

PSF = @(x, y, z, x0, y0, z0, A, sigma_narrow, sigma_wide, sigma_axial) ...
        A.* ( exp(-(x - x0).^2 ./ (2*sigma_narrow^2) - (y - y0).^2 ./ (2*sigma_wide^2)) ...
        + exp(-(x - x0).^2 ./ (2*sigma_wide^2) - (y - y0).^2 ./ (2*sigma_narrow^2)) ) ...
        .* exp(-(z - z0).^2 ./ (2*sigma_axial^2));

% test = PSF(X, Y, Z, x0, y0, z0, 1, sigma_narrow, sigma_wide, sigma_axial);
% volumeViewer(test)
parametrized_fun = @(param_vec, X_design_matrix) PSF(X_design_matrix(:, 1), X_design_matrix(:, 2), X_design_matrix(:, 3), x0, y0, z0, param_vec(1), param_vec(2), param_vec(3), param_vec(4));

testfit = nlinfit([X(:), Y(:), Z(:)], abs(IQ(:)), parametrized_fun, [max(abs(IQ), [], 'all'), sigma_narrow_test, sigma_wide_test, sigma_axial_test]);

testPSF = reshape(parametrized_fun(testfit, [X(:), Y(:), Z(:)]), size(IQ));

sigma_um = testfit(2:end).*1e6

%% Test plots
figure; yyaxis left; plot(squeeze(abs(IQ(pos_rc(1), :, pos_rc(3))))); yyaxis right; plot(squeeze(abs(testPSF(pos_rc(1), :, pos_rc(3)))))

%%
% %% Define some parameters for the Gaussian PSF fit
% gfit_pixel_spacing = 0.01;
% fit_type = 'gauss2'; % Use a two-term Gaussian for the fit
% % fit_type = 'gauss1'; % Use a one-term Gaussian for the fit
% 
% %% Get x PSF and plot
% xPSF = squeeze(abs(IQ(pos_rc(1), :, pos_rc(3)))); % 1D x PSF
% xPSF_pixel_inds = 1:length(xPSF);
% [xPSF_GF, xPSF_GF_values, xPSF_pixel_inds_GF] = PSFGaussianFit(xPSF_pixel_inds, xPSF, gfit_pixel_spacing, fit_type); % Fit a Gaussian to the x PSF
% 
% % Plot
% figure; plot(xPSF_pixel_inds_GF, xPSF_GF_values, 'LineWidth', 2)
% hold on
% plot(1:length(xPSF), xPSF, ':', 'LineWidth', 2)
% hold off
% title('x PSF and Gaussian fit')
% legend('Gaussian fit', 'Simulated')
% 
% %% Get y PSF and plot
% yPSF = squeeze(abs(IQ(:, pos_rc(2), pos_rc(3)))); % 1D y PSF
% yPSF_pixel_inds = 1:length(yPSF);
% [yPSF_GF, yPSF_GF_values, yPSF_pixel_inds_GF] = PSFGaussianFit(yPSF_pixel_inds, yPSF, gfit_pixel_spacing, fit_type); % Fit a Gaussian to the x PSF
% 
% % Plot
% figure; plot(yPSF_pixel_inds_GF, yPSF_GF_values, 'LineWidth', 2)
% hold on
% plot(1:length(yPSF), yPSF, ':', 'LineWidth', 2)
% hold off
% title('y PSF and Gaussian fit')
% legend('Gaussian fit', 'Simulated')
% 
% %% Get z PSF and plot
% zPSF = squeeze(abs(IQ(pos_rc(1), pos_rc(2), :))); % 1D z PSF
% zPSF_pixel_inds = 1:length(zPSF);
% [zPSF_GF, zPSF_GF_values, zPSF_pixel_inds_GF] = PSFGaussianFit(zPSF_pixel_inds, zPSF, gfit_pixel_spacing, fit_type); % Fit a Gaussian to the x PSF
% 
% % Plot
% figure; plot(zPSF_pixel_inds_GF, zPSF_GF_values, 'LineWidth', 2)
% hold on
% plot(1:length(zPSF), zPSF, ':', 'LineWidth', 2)
% hold off
% title('z PSF and Gaussian fit')
% legend('Gaussian fit', 'Simulated')
% 
% %% Get 1/e and FWHM values for x, y, z PSFs by using the Gaussian fits
% % FWHM values
% [FWHM_GF_units(1)] = fwhm(xPSF_pixel_inds_GF, xPSF_GF_values);
% [FWHM_GF_units(2)] = fwhm(yPSF_pixel_inds_GF, yPSF_GF_values);
% [FWHM_GF_units(3)] = fwhm(zPSF_pixel_inds_GF, zPSF_GF_values);
% 
% % FWHM_wl = FWHM_GF_units .* PData.PDelta ./ gfit_pixel_spacing;
% FWHM_wl = FWHM_GF_units .* PData.PDelta;
% % FWHM_um = FWHM_wl .* P.wl ./ 1e6;
% FWHM_um = FWHM_wl .* P.wl .* 1e6;
% 
% % Sigma values (field)
% if strcmp(fit_type, 'gauss1')
%         % Get directly from the fit parameters if we use a single Gaussian
%         %   Matlab's Gauss1 model form: g(x) = a1*exp(-((x-b1)/c1)^2) --> we want sigma = c1 / sqrt(2)
%         [sigma_field_GF_units(1)] = xPSF_GF.c1 / sqrt(2);
%         [sigma_field_GF_units(2)] = yPSF_GF.c1 / sqrt(2);
%         [sigma_field_GF_units(3)] = zPSF_GF.c1 / sqrt(2);
% else % Otherwise, calculate it manually
%     [sigma_field_GF_units(1)] = fw_anymax(xPSF_pixel_inds_GF, xPSF_GF_values, 1/exp(1)) / 2 / sqrt(2); % fw_anymax is the full-width, so divide by 2, and then by sqrt(2) to get sigma
%     [sigma_field_GF_units(2)] = fw_anymax(yPSF_pixel_inds_GF, yPSF_GF_values, 1/exp(1)) / 2 / sqrt(2);
%     [sigma_field_GF_units(3)] = fw_anymax(zPSF_pixel_inds_GF, zPSF_GF_values, 1/exp(1)) / 2 / sqrt(2);
% 
% end
% 
% % sigma_field_wl = sigma_field_GF_units .* PData.PDelta ./ gfit_pixel_spacing;
% sigma_field_wl = sigma_field_GF_units .* PData.PDelta;
% sigma_field_um = sigma_field_wl .* P.wl .* 1e6;
% 
% % Sigma values (intensity)
% sigma_intensity_um = sigma_field_um ./ sqrt(2);