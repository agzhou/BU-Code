

clearvars
%% load params and stuff
IQpath = uigetdir('D:\Allen\Data\', 'Select the IQ data path');
IQpath = [IQpath, '\'];

% Load parameters
% if ~exist('P', 'var')
%     load([IQpath, '..\params.mat'])
% end
% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [IQpath, '..\params.mat']);
    load([params_pathname, params_filename])
end

% Load Verasonics reconstruction parameters: datapath\PData.mat
if ~exist('PData', 'var')
    load([IQpath, 'PData.mat'])
end

IQfilenameStructure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

%%
% addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Processing\Speckle tracking')
addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Simulations')
%% Load data
filenum = 1
load([IQpath, IQfilenameStructure, num2str(filenum)])

IQ = squeeze(IData + 1i .* QData);
clearvars IData QData

%%
%     volumeViewer(squeeze(abs(IQ(:, :, :, 1))))
% figure; imagesc(squeeze(max( abs(IQ(:, :, :, 1)), [], 1) )')
figure; imagesc(squeeze( abs(IQ(40, :, :, 1)) )')
%% Calculate g1T
num_g1_pts = 100;
g1 = sim_g1T(IQ, num_g1_pts);

%% Calculate time lags
taustep = 1/P.frameRate;
% tau = taustep:taustep:(P.numFramesPerBuffer * taustep);
tau = 0:taustep:((P.numFramesPerBuffer - 1) * taustep);
tau_ms = tau .* 1000; % Assuming even time spacing between frames

%%
% figure; plot(squeeze(abs(IQ(40, 42, 62, :))))
% figure; plot(tau_ms(1:num_g1_pts), squeeze(abs(g1(40, 40, 60, :)))); xlabel("Tau [ms]"); ylabel("|g_1|")
% figure; plot(squeeze(abs(g1(40, 45, 62, 1:50))))

%% Extract the vessel only - theoretically
PDI = sum(abs(IQ) .^ 2, 4) ./ size(IQ, 4);

% figure; histogram(PDI, 'Normalization', 'probability')
%%
frac = 0.5;
nonvessel_mask = PDI < frac * max(PDI, [], 'all'); % Mask of voxels "not in vessel"
vessel_mask = ~nonvessel_mask; % Mask of voxels "in vessel"
PDI_fracmax = PDI; PDI_fracmax(nonvessel_mask) = 0;

figure; imagesc(squeeze( max(PDI, [], 1) )')
figure; imagesc(squeeze( max(PDI_fracmax, [], 1) )')
% volumeViewer(PDI_fracmax)
%% Go through voxels "in vessel" and average the g1
g1_vessel_avg = calc_ROI_avg(g1, vessel_mask);
% figure; plot(squeeze(abs(g1_vessel_avg)))

% IQ_vessel_avg = calc_ROI_avg(IQ - mean(IQ, 4), vessel_mask);
%% Plot the average g1 "in vessel"
g1af = figure; plot(tau_ms(1:num_g1_pts), squeeze(abs(g1_vessel_avg)), 'LineWidth', 4); xlabel("Tau [ms]"); ylabel("|g_1|")
title("|g_1| average across vessel; flow speed = " + num2str(P.Mcr_SP.flow_v_mm_s) + " mm/s")
fontsize(g1af, 20, 'points')

%% Calculate CBV index and CBFspeed index
tau1_index_CBF = 2;
tau2_index_CBF = 3;
tau1_index_CBV = 2;

[CBFsi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
[CBFsi_va, CBVi_va] = g1_to_CBi(g1_vessel_avg, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)


% %%
% g1_tau1_cutoff = 0.25;
% % g1_tau1_cutoff = 0.1;
% 
% 
% [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);
% % CBFsi(~g1A_mask) = 0; % Remove noisy points from the CBFspeed index (in theory)

%% Plot CBV and CBFspeed indices
figure; imagesc(squeeze( max(CBVi, [], 1) )')
figure; imagesc(squeeze( max(CBFsi, [], 1) )')
% volumeViewer(CBVi)





%% Nonlinear Least Squares fitting of |g1T|

% Define a maximum tau to fit to (improve accuracy)
tau_max = 40 / 1e3; % [s]
tau_mask = tau < tau_max;
tau_range = tau(tau_mask);

% Define a function handle for |g1T|
% % Let x = [Ns, vx, vy, vz]
% Let x = [vx, vy, vz]

% Sigma is the 1/e * PSF max, so convert from FWHM
P.Mcr_SP.sigma = [300e-6, 300e-6, 150e-6] .* 1/(2*sqrt(2*log(2))); %%%% PSF testing %%%%
% P.Mcr_SP.sigma = [300e-6, 300e-6, 100e-6] .* 1/(2*sqrt(2*log(2))); %%%% PSF testing %%%%
% P.Mcr_SP.sigma = [300e-6, 300e-6, 10e-6]; %%%% PSF testing %%%%

Rs = P.Mcr_SP.scatterReflectivity;
% Rs = mean(voxel.data(:, 4)); % Average reflection coefficient of scatterer
% Re = SP.snr; % "Noise level of imaging system"
Re = 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re = 1e-8
% abs_g1T_model = @(x, tau) x(1)*Rs*(pi^1.5)*P.Mcr_SP.sigma(1)*P.Mcr_SP.sigma(2)*P.Mcr_SP.sigma(3) ...
%                           / ( x(1)*Rs*(pi^1.5)*P.Mcr_SP.sigma(1)*P.Mcr_SP.sigma(2)*P.Mcr_SP.sigma(3) + Re ) ...
%                           .* exp( -(x(2) .* tau).^2./(4.*P.Mcr_SP.sigma(1)^2) -(x(3) .* tau).^2./(4.*P.Mcr_SP.sigma(2)^2) -(x(4) .* tau).^2./(4.*P.Mcr_SP.sigma(3)^2) );
abs_g1T_model = @(x, tau) exp( -(x(1) .* tau).^2./(4.*P.Mcr_SP.sigma(1)^2) -(x(2) .* tau).^2./(4.*P.Mcr_SP.sigma(2)^2) -(x(3) .* tau).^2./(4.*P.Mcr_SP.sigma(3)^2) );


% Fit the unknown parameters (x) with the simulation and compare to the ground truth
% TESTING: INITIAL GUESS AND UPPER/LOWER BOUNDS
testx = [0, 0, P.Mcr_SP.flow_v_mm_s/1e3];
x0 = testx %%%%%%%% testing
% lb = [0, 0, 0, 0]
% ub = [Inf, 0, 0, 1000e-3]
lb = [0, 0, -1000]
ub = [0, 0, 1000e-3]
% x_fit = lsqcurvefit(abs_g1T_model, x0, tau, abs(voxel.g1))
% x_fit = lsqcurvefit(abs_g1T_model, x0, tau, abs(voxel.g1), lb, ub)
x_fit = lsqcurvefit(abs_g1T_model, x0, tau_range, abs(g1_vessel_avg(tau_mask)), lb, ub)

%% Input the ground truth parameters to see what the |g1T| model looks like

test_abs_g1T_model = abs_g1T_model(testx, tau_range);
figure; plot(tau_range*1e3, test_abs_g1T_model, 'LineWidth', 2); xlabel('tau [ms]'); ylabel("|g_1|")
hold on; plot(tau_range .* 1e3, abs(g1_vessel_avg(tau_mask)), 'LineWidth', 2);

% Plot the fit
abs_g1T_fit = abs_g1T_model(x_fit, tau_range);
plot(tau_range .* 1e3, abs_g1T_fit, 'LineWidth', 2);
% test2 = abs_g1T_model([7849, 0, 0, .0445], tau);
% test2 = abs_g1T_model([7849, 0, 0, .0845], tau);
% plot(tau .* 1e3, test2);
hold off
legend("|g_1| model with ground truth parameters", "|g_1| from ultrasound simulation", "Fit")