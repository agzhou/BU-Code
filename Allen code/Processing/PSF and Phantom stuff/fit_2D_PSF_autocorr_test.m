%% Load the PSF
datapath = uigetdir('G:\', 'Select the PSF data path');
load([datapath, '\PSF.mat'])
load([datapath, '\params.mat'])
load([datapath, '\PData.mat'])

% Define coordinates
[x_mm, z_mm] = getReconCoords2D(PData, P);

[zp, xp] = size(IQ);

%% Get the position of the point scatterer - assuming there is only one scatterer
pos_wl = P.Media.MP(1, 1:3); % Position in wavelengths
pos_m = pos_wl .* P.wl; % [x, y, z] position in meters

% Get the position of the point scatterer in terms of the reconstructed data coordinates
% This abs(PData.Origin) line is kind of a hack.. fix later if needed
pos_rc = round((pos_wl + abs(PData.Origin)) ./ PData.PDelta + 1); % Position in recon coords (add 1 because of Matlab indexing)
% ^ [x, y, z]

%% Use the PSF's spatial autocorrelation to calculate the g1,s(tau) - "single/uniform scatterer"

dzs = 0:10; % Vector of axial displacements [pixels] to calculate g1,s at
dxs = 0:10; % Vector of lateral displacements [pixels] to calculate g1,s at
dzs_m = dzs .* PData.PDelta(3)*P.wl; % Convert to m
dxs_m = dxs .* PData.PDelta(1)*P.wl;

zw = 20; % half width of the axial window to consider
xw = 20; % half width of the lateral window to consider

% check if the window will run past the data's total size
if pos_rc(3) + zw + max(dzs) > zp | pos_rc(3) - zw - max(dzs) < 1
    error('Window will run out of the valid range of the PSF data')
elseif pos_rc(1) + xw + max(dxs) > xp | pos_rc(1) - xw - max(dxs) < 1
    error('Window will run out of the valid range of the PSF data')
end

%% Calculate the autocorrelations
% Slice the centered window of the PSF
PSF_centered = IQ(pos_rc(3) - zw : pos_rc(3) + zw, pos_rc(1) - xw : pos_rc(1) + xw);

% First, only axial shifts
g1sd_z = zeros(length(dzs), 1); % g1,s as calculated from data, for axial shifts
for dzi = 1:length(dzs)
    dz = dzs(dzi);
    PSF_axial_shifted = IQ(dz + (pos_rc(3) - zw : pos_rc(3) + zw), pos_rc(1) - xw : pos_rc(1) + xw);
    g1sd_z(dzi) = sum(PSF_centered .* conj(PSF_axial_shifted), 'all') ./ sum(PSF_centered .* conj(PSF_centered), 'all');
end   
    
% Next, only lateral shifts
g1sd_x = zeros(length(dxs), 1); % g1,s as calculated from data, for lateral shifts
for dxi = 1:length(dxs)
    dx = dxs(dxi);
    PSF_lateral_shifted = IQ(pos_rc(3) - zw : pos_rc(3) + zw, dx + (pos_rc(1) - xw : pos_rc(1) + xw));
    g1sd_x(dxi) = sum(PSF_centered .* conj(PSF_lateral_shifted), 'all') ./ sum(PSF_centered .* conj(PSF_centered), 'all');
end   

%% Fit sigmas to the measured decorrelation curves
opts = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-14, 'StepTolerance', 1e-12);
lb = 1e-6; ub = 1000e-6;
% Fit axial
[sigma_z] = lsqnonlin(@(sigma_z) g1_model(sigma_z, dzs_m) - abs(g1sd_z).', 50e-6, lb, ub, opts);

% Fit lateral
[sigma_x] = lsqnonlin(@(sigma_x) g1_model(sigma_x, dxs_m) - abs(g1sd_x).', 50e-6, lb, ub, opts);

%% Visualize the fits
% Calculate g1 curves according to the model fit
g1sm_z = g1_model(sigma_z, dzs_m);
g1sm_x = g1_model(sigma_x, dxs_m);

% Plot
figure; plot(dzs, abs(g1sd_z), dzs, g1sm_z)
figure; plot(dxs, abs(g1sd_x), dxs, g1sm_x)

%% test
% test = g1_model(50e-6, dzs_m);
% figure; plot(test)
%% Helper functions

function [g1se] = g1_model(x, ds) % assumes dz = 0 when fitting laterally, or dx = 0 when fitting axially
    sigma = x(1);
    g1se = exp(-ds.^2./(4.*sigma.^2)); % Envelope only; the phase is fit later on in the vUS
end


