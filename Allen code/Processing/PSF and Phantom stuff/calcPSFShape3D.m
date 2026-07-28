%% Description: semi-automatically calculate the ultrasound PSF shape in 3D, taking in simulation data

%% Load simulation data and parameters (assuming Verasonics)
load('G:\My Drive\Data\RC15gV PSF sim - 11 angles from -5 to 5 deg\params.mat')
load('G:\My Drive\Data\RC15gV PSF sim - 11 angles from -5 to 5 deg\PData.mat')
load('G:\My Drive\Data\RC15gV PSF sim - 11 angles from -5 to 5 deg\PSF.mat')

%% Plot the data
figure; imagesc(squeeze(max(abs(IQ), [], 1))'); title('xz MIP'); colormap gray

%% Get the position of the point scatterer - assuming there is only one scatterer
pos_wl = P.Media.MP(1, 1:3); % Position in wavelengths
pos_m = pos_wl .* P.wl; % [x, y, z] position in meters

%% Get the position of the point scatterer in terms of the reconstructed data coordinates

% This abs(PData.Origin) line is kind of a hack.. fix later if needed
pos_rc = round((pos_wl + abs(PData.Origin)) ./ PData.PDelta + 1); % Position in recon coords (add 1 because of Matlab indexing)

%% Define some parameters for the Gaussian PSF fit
gfit_pixel_spacing = 0.01;


%% Get x PSF and plot
xPSF = squeeze(abs(IQ(pos_rc(1), :, pos_rc(3)))); % 1D x PSF
xPSF_pixel_inds = 1:length(xPSF);
[xPSF_GF, xPSF_GF_values, xPSF_pixel_inds_GF] = PSFGaussianFit(xPSF_pixel_inds, xPSF, gfit_pixel_spacing); % Fit a Gaussian to the x PSF

% Plot
figure; plot(xPSF_pixel_inds_GF, xPSF_GF_values, 'LineWidth', 2)
hold on
plot(1:length(xPSF), xPSF, ':', 'LineWidth', 2)
hold off
title('x PSF and Gaussian fit')
legend('Gaussian fit', 'Simulated')

%% Get y PSF and plot
yPSF = squeeze(abs(IQ(:, pos_rc(2), pos_rc(3)))); % 1D y PSF
yPSF_pixel_inds = 1:length(yPSF);
[yPSF_GF, yPSF_GF_values, yPSF_pixel_inds_GF] = PSFGaussianFit(yPSF_pixel_inds, yPSF, gfit_pixel_spacing); % Fit a Gaussian to the x PSF

% Plot
figure; plot(yPSF_pixel_inds_GF, yPSF_GF_values, 'LineWidth', 2)
hold on
plot(1:length(yPSF), yPSF, ':', 'LineWidth', 2)
hold off
title('y PSF and Gaussian fit')
legend('Gaussian fit', 'Simulated')

%% Get z PSF and plot
zPSF = squeeze(abs(IQ(pos_rc(1), pos_rc(2), :))); % 1D z PSF
zPSF_pixel_inds = 1:length(zPSF);
[zPSF_GF, zPSF_GF_values, zPSF_pixel_inds_GF] = PSFGaussianFit(zPSF_pixel_inds, zPSF, gfit_pixel_spacing); % Fit a Gaussian to the x PSF

% Plot
figure; plot(zPSF_pixel_inds_GF, zPSF_GF_values, 'LineWidth', 2)
hold on
plot(1:length(zPSF), zPSF, ':', 'LineWidth', 2)
hold off
title('z PSF and Gaussian fit')
legend('Gaussian fit', 'Simulated')

%% Get 1/e and FWHM values for x, y, z PSFs by using the Gaussian fits
% FWHM values
[FWHM_GF_units(1)] = fwhm(xPSF_pixel_inds_GF, xPSF_GF_values);
[FWHM_GF_units(2)] = fwhm(yPSF_pixel_inds_GF, yPSF_GF_values);
[FWHM_GF_units(3)] = fwhm(zPSF_pixel_inds_GF, zPSF_GF_values);

% FWHM_wl = FWHM_GF_units .* PData.PDelta ./ gfit_pixel_spacing;
FWHM_wl = FWHM_GF_units .* PData.PDelta;
% FWHM_um = FWHM_wl .* P.wl ./ 1e6;
FWHM_um = FWHM_wl .* P.wl .* 1e6;

% 1/e values
[OOE_GF_units(1)] = fw_anymax(xPSF_pixel_inds_GF, xPSF_GF_values, 1/exp(1));
[OOE_GF_units(2)] = fw_anymax(yPSF_pixel_inds_GF, yPSF_GF_values, 1/exp(1));
[OOE_GF_units(3)] = fw_anymax(zPSF_pixel_inds_GF, zPSF_GF_values, 1/exp(1));

% OOE_wl = OOE_GF_units .* PData.PDelta ./ gfit_pixel_spacing;
OOE_wl = OOE_GF_units .* PData.PDelta;
% OOE_um = OOE_wl .* P.wl ./ 1e6;
OOE_um = OOE_wl .* P.wl .* 1e6;