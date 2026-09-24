%% Description: semi-automatically calculate the ultrasound PSF shape in 3D, taking in simulation data
clearvars

%% Add the Processing folder to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
AllenProcessingCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing\");
addpath(genpath(AllenProcessingCodePath))

%% Load simulation data and parameters (assuming Verasonics)
% datapath = "G:\My Drive\Data\PSF Simulations\L22-14v PSF sim - 17 angles from -10 to 10 deg\";
datapath = uigetdir('G:\', 'Select the data path');
datapath = [datapath, '\'];

load([datapath, 'params.mat'])
load([datapath, 'PData.mat'])
load([datapath, 'PSF.mat'])

%% Define coordinates
[x_mm, z_mm] = getReconCoords2D(PData, P);

%% Plot the data
figure; imagesc(x_mm, z_mm, abs(IQ)); title('xz plane'); colormap gray; xlabel('x [mm]'); ylabel('z [mm]'); axis equal; axis tight

%% Get the position of the point scatterer - assuming there is only one scatterer
pos_wl = P.Media.MP(1, 1:3); % Position in wavelengths
pos_m = pos_wl .* P.wl; % [x, y, z] position in meters

%% Get the position of the point scatterer in terms of the reconstructed data coordinates

% This abs(PData.Origin) line is kind of a hack.. fix later if needed
pos_rc = round((pos_wl + abs(PData.Origin)) ./ PData.PDelta + 1); % Position in recon coords (add 1 because of Matlab indexing)

%% Define some parameters for the Gaussian PSF fit
gfit_pixel_spacing = 0.01;


%% Get x PSF and plot
xPSF = squeeze(abs(IQ(pos_rc(3), :))); % 1D x PSF
xPSF_pixel_inds = 1:length(xPSF);
[xPSF_GF, xPSF_GF_values, xPSF_pixel_inds_GF] = PSFGaussianFit(xPSF_pixel_inds, xPSF, gfit_pixel_spacing); % Fit a Gaussian to the x PSF

% Plot
figure; plot(xPSF_pixel_inds_GF, xPSF_GF_values, 'LineWidth', 2)
hold on
plot(1:length(xPSF), xPSF, ':', 'LineWidth', 2)
hold off
title('x PSF and Gaussian fit')
legend('Gaussian fit', 'Simulated')

%% Get z PSF and plot
zPSF = squeeze(abs(IQ(:, pos_rc(1)))); % 1D z PSF
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
% FWHM values: calculate from the fitted Gaussian's curve
[FWHM_GF_units(1)] = fwhm(xPSF_pixel_inds_GF, xPSF_GF_values);
[FWHM_GF_units(2)] = NaN;
[FWHM_GF_units(3)] = fwhm(zPSF_pixel_inds_GF, zPSF_GF_values);

% % FWHM values: calculate directly from the fit parameters
% [FWHM_GF_units(1)] = xPSF_GF.c1 * _;
% [FWHM_GF_units(2)] = NaN;
% [FWHM_GF_units(3)] = zPSF_GF.c1 * _;

% FWHM_wl = FWHM_GF_units .* PData.PDelta ./ gfit_pixel_spacing;
FWHM_wl = FWHM_GF_units .* PData.PDelta;
FWHM_um = FWHM_wl .* P.wl .* 1e6;

% sigma (field): calculate from the fitted Gaussian's curve 1/e values
% [sigma_field_GF_units(1)] = fw_anymax(xPSF_pixel_inds_GF, xPSF_GF_values, 1/exp(1)) / 2 / sqrt(2);
% [sigma_field_GF_units(2)] = NaN;
% [sigma_field_GF_units(3)] = fw_anymax(zPSF_pixel_inds_GF, zPSF_GF_values, 1/exp(1)) / 2 / sqrt(2);

% sigma (field): get directly from the fit parameters
%   Matlab's Gauss1 model form: g(x) = a1*exp(-((x-b1)/c1)^2) --> we want sigma = c1 / sqrt(2)
[sigma_field_GF_units(1)] = xPSF_GF.c1 / sqrt(2);
[sigma_field_GF_units(2)] = NaN;
[sigma_field_GF_units(3)] = zPSF_GF.c1 / sqrt(2);

% sigma_field_wl = sigma_field_GF_units .* PData.PDelta ./ gfit_pixel_spacing;
sigma_field_wl = sigma_field_GF_units .* PData.PDelta;
sigma_field_um = sigma_field_wl .* P.wl .* 1e6;

% sigma values (intensity)
sigma_intensity_um = sigma_field_um ./ sqrt(2);

