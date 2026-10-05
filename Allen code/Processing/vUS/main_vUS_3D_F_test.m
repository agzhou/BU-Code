%% Description:
%   Calculate flow velocity from volumetric fUS data using the vUS method
%   (Tang et al., 2020)

% Test whether fitting certain parameters results in a significantly better
% fit, with an F-test

%% Only use if needed: convert Jianbo/Bingxue's acquisition parameters to something I can use
% First: manually load an IQ file, like: load('E:\PROJ_tlfUS\IQdata\0806_2021_BL3_vUS_run1(good)\IQ-10-5-5000-1000-1-BL3-1.mat')
% if ~exist('P_old', 'var')
%     P_old = P; clearvars P
% end
% [P] = oldP2P(P_old);

%% Add the Processing folder to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
AllenProcessingCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing\");
addpath(genpath(AllenProcessingCodePath))
% ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
% addpath(genpath(ErrorFunctionCodePath))

%% Load data and things
if ~exist('PDpath', 'var')
    PDpath = uigetdir('D:\Allen\Data\', 'Select the Processed Data path');
    PDpath = [PDpath, '\'];
end

% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [PDpath, '..\params.mat']);
    load([params_pathname, params_filename])
end
clearvars params_filename params_pathname

% Choose and load the processing parameters: fUS_proc_params.mat file (from IQ2g1T_3D_loop.m)
[proc_params_filename, proc_params_pathname, ~] = uigetfile('*.mat', 'Select the processing params file', [PDpath, '\fUS_proc_params.mat']);
load([proc_params_pathname, proc_params_filename])
clearvars proc_params_filename proc_params_pathname

%% (Only if using continuous data) load some chunk of continuous superframes
% IQfilenameStructure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];
% 
% startFile = 1; endFile = 1;
% nsfpe = (endFile - startFile + 1)*P.numFramesPerBuffer;% Number of superframes per ensemble
% load([IQpath, IQfilenameStructure, num2str(startFile)])
% IQle = []; % long ensemble IQ
% 
% for fi = startFile:endFile
%     load([IQpath, IQfilenameStructure, num2str(fi)])
%     IQle = cat(3, IQle, IQ);
% end
% 
% IQ = IQle; clearvars IQle



%% Define some parameters

% Define the sigma values, from simulations
%   Single-Gaussian version
% sigma_SG = [121.8936, 121.8936, 44.3846].*1e-6; % Field-based sigma values (x, y, z) [m] for the RC15gV probe at 13.6 MHz and 11 x 2 angles from -5 to 5 deg (G:\My Drive\Data\RC15gV PSF sim - 11 angles from -5 to 5 deg)
% sigma_SG = [151.5410, 151.5410, 44.6252].*1e-6; % Field-based sigma values (x, y, z) [m] for the RC15gV probe at 13.6 MHz and 11 x 2 angles from -6 to 6 deg
sigma_SG = [130.1474, 130.1474, 44.7070].*1e-6; % Field-based sigma values (x, y, z) [m] for the RC15gV probe at 13.6 MHz and 5 x 2 angles from -6 to 6 deg

%   RCA-specific orthogonal Gaussians version
sigma_RCA = [57.6, 286.5, 50.1].*1e-6; % Field-based sigma values (narrow, wide, axial) [m] for the RC15gV probe at 13.6 MHz: 5 angles from -6 to 6 deg

xDim = 1; % Dimension of the data corresponding to x (lateral direction)
yDim = 2; % Dimension of the data corresponding to y (lateral direction)
zDim = 3; % Dimension of the data corresponding to z (axial direction)
fDim = 4; % Dimension of the data corresponding to frequency (or time)

if ~exist('ctp', 'var')
    ctp = 1:3; % Indices of which frequency Components To Process (typically [1, 2, 3]: negative, positive, all)
    ctp_labels = {"Down flows", "Up flows", "All flows"};
end

% Load Jianbo's colormaps
[VzCmap, VzCmapDn, VzCmapUp, pdiCmapUp, PhtmCmap] = Colormaps_fUS;

%% Create a struct for all the relevant processing parameters
dimensionality = 3; % 2D data
frameRate = P.frameRate;
wl = P.wl;
k0 = 2*pi/wl;
% PP = createStruct(xp, yp, zp, nf, nTau, xDim, yDim, zDim, fDim, dimensionality, faxis, freqMask, frameRate, wl, k0); % Processing Parameters ======> adjust as needed
PP = createStruct(xp, yp, zp, nf, nTau, xDim, yDim, zDim, fDim, dimensionality, frameRate, wl, k0); % Processing Parameters ======> adjust as needed

%% Load the g1 average across superframes and fit it

% Choose and load the g1 superframe average: g1_avg.mat file (from calc_g1_avg.m)
[g1_avg_filename, g1_avg_pathname, ~] = uigetfile('*.mat', 'Select the g1 average file', [PDpath, '\g1_avg.mat']);
load([g1_avg_pathname, g1_avg_filename], 'g1_avg')
g1 = g1_avg; clear g1_avg;

vs = size(g1{3}); vs = vs(1:end-1); % Volume size [voxels]
num_voxels = prod(vs);

t1i = 2; % Index for tau1 --> 2 for my code, because it calculates g1 starting at tau = 0


% for j = ctp
%     g1{j} = double(g1{j});
% end

% voxelTimeseriesGUI(g1{3}, tempPDI.^0.5, 'ComplexMode', 'abs')

% Create new variables for experimental g1, with spatial dimensions stacked
g1_exp = cell(size(g1)); % Cell array of experimental g1 data with spatial dimensions vectorized/stacked
num_voxels = size(g1{3}, xDim)*size(g1{3}, yDim)*size(g1{3}, zDim);
for j = ctp
    g1_exp{j} = reshape(g1{j}, num_voxels, nTau);
end

% ---- Loop through directional components and go through the fitting process ---- %
% for j = ctp
% for j = 1:2 % Fit only negative and positive frequencies (down and up flows)
for j = 3
    
    % Axial component of the blood flow's group velocity -- v_zgp
    [Vz0, tau_V] = findVzPhaseDiff(stackData(g1{j}, PP), PP); % v_zgp [m/s]
    % figure; imagesc(squeeze(max(unstackData(Vz0, PP), [], 1))'); colormap(VzCmap); axis equal; colorbar; clim([-30e-3, 30e-3])
    % volumeViewer(abs(unstackData(Vz0, PP)))

    % Mesh method for finding v_tgp0
    % [v_zgp0, v_xgp0, p0, DC0, F0, R20] = InitvUS2DParamsWithMesh(g1adj_stacked_j, Vz0, DCR0_j, FR_j, PP, sigma, tau);
    % [Vx0, R2_Vx0] = InitVx0WithMesh2D(stackData(g1{j}, PP), Vz0, DCR0_j, FR0_j, PP, sigma, tau);
    Vt0 = ones(size(Vz0)).* 5e-3; % TESTING: uniform initial v_xgp guess
    % figure; imagesc(unstackData(Vx0, PP)); axis equal; colorbar

    % **** TO DO: create a function that looks at the confidence in the
    % angle-based Vt0, and outputs an updated Vt0 if needed, plus searches
    % for the best p0 ****
    
    % Set parameters for findTauDecayed.m
    fp.absolute_tau_ss_cutoff_s = 10e-3;
    fp.too_fast_decay_s = 2e-3;

    % Adaptively find the tau range to fit over for each pixel, and another
    % quality mask
    [tau_decayed_ind, voxel_quality] = findTauDecayed(g1_exp{j}, tau, t1i, fp.absolute_tau_ss_cutoff_s, fp.too_fast_decay_s);
    % volumeViewer(unstackData(voxel_quality, PP))

    % Another mask
    temp_g1_tau1_mask = squeeze(abs(g1_exp{j}(:, t1i))) > 0.3; % Note: this does not account for the static component
    % volumeViewer(unstackData(temp_g1_tau1_mask, PP))

    temp_g1_tau2_mask = squeeze(abs(g1_exp{j}(:, t1i + 1))) > 0.2; % Note: this does not account for the static component
    % volumeViewer(unstackData(temp_g1_tau2_mask, PP))

    tempMask = and(temp_g1_tau1_mask, temp_g1_tau2_mask);

    % Struct for storing parameters of the different g1 models
    ifpv = zeros(num_voxels, 1); % initial fit parameter value matrix
    fpn_SG = {"v_tgp", "v_zgp", "F", "DC", "k", "a"}; % Fit parameter names (single Gaussian PSF model)
    fit_SG = initFitParamStruct(fpn_SG, ifpv);
    fpn_RCA = {"v_xgp", "v_ygp", "v_zgp", "F", "DC", "k", "a"}; % Fit parameter names (RCA-specific PSF model)
    fit_RCA = initFitParamStruct(fpn_RCA, ifpv);
    % fpn_C = {"C", "v_zgp", "F", "DC"}; % Fit parameter names (Combined parameter model)
    % fit_C = initFitParamStruct(fpn_C, ifpv);
    
    % Choose the mask to use to fit certain pixels or not
    % maskToUse = overall_mask_stacked; 
    % maskToUse = and(vesselAngleMask, stackData(temp_g1_tau1_mask, PP));
    % maskToUse = and(voxel_quality, tempMask);
    maskToUse = or(voxel_quality, tempMask);
    % maskToUse = voxel_quality;
    % maskToUse = vesselMask;
    % volumeViewer(unstackData(maskToUse, PP))

    maskToUseTrueInds = find(maskToUse).'; % indices where maskTouse is true
    
    % ---- Use batched solver ---- %
    [s, w] = gaussLegendre01(48);
    % Make vectors for bounds for fit parameters (and store in struct 'fb')
    fb.v_xgp = [0, 250e-3];
    fb.v_ygp = [0, 250e-3];
    fb.v_tgp = sqrt(fb.v_xgp.^2 + fb.v_ygp.^2);
    fb.v_zgp = [-50e-3; 50e-3];
    fb.F = [0, 1];
    fb.DC = [0, 1];
    fb.k = [2, 3];
    fb.a = [0, 1];
    % fb.C = [0, Inf];

    % Single-Gaussian PSF model
    tic
    fit_SG.lb = [fb.v_tgp(1), fb.v_zgp(1), fb.F(1), fb.DC(1), fb.k(1), fb.a(1)]; % Parameter lower bounds [v_tgp, v_zgp, F, DC, k, a]
    fit_SG.ub = [fb.v_tgp(2), fb.v_zgp(2), fb.F(2), fb.DC(2), fb.k(2), fb.a(2)]; % Parameter upper bounds
    
    fit_SG.xscale = [1e-2 1e-2 1 1 1 1];
    fit_SG.model = @(X, cols) vUS_3D_quad_batchModel(X, tau(cols), PP.k0, sigma_SG, s, w);
    % vi   = find(maskToUse);  
    nv = numel(maskToUseTrueInds);
    fit_SG.opts   = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', 2, 'jacobian', 'analytic', 'xscale', fit_SG.xscale);
    fit_SG.x0   = [Vt0(maskToUseTrueInds), Vz0(maskToUseTrueInds), ones(nv, 1), zeros(nv, 1), 2.5*ones(nv, 1), 0.7 .*ones(nv, 1)];
    fit_SG.opts = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', 2, 'jacobian', 'analytic', 'xscale', fit_SG.xscale);
    [fit_SG.x, fit_SG.cost, ~, fit_SG.conv] = fitBatchedLM(fit_SG.model, fit_SG.x0, g1_exp{j}(maskToUseTrueInds, :), fit_SG.lb, fit_SG.ub, fit_SG.opts);    
    % [fit_SG.x, fit_SG.cost] = vUS_3D_quad_fitBatched(g1_exp{j}(maskToUseTrueInds,:), tau_decayed_ind(maskToUseTrueInds), Vz0(maskToUseTrueInds), tau, PP.k0, sigma_SG, fit_SG.lb, fit_SG.ub);
    [fit_SG] = storeFitParams(fit_SG, fpn_SG, fit_SG.x, maskToUseTrueInds, PP); % Store/parse fitted parameters
    toc

    % RCA-specific PSF model
    tic
    fit_RCA.xscale = [1e-2 1e-2 1e-2 1 1 1 1];
    fit_RCA.lb = [fb.v_xgp(1), fb.v_ygp(1), fb.v_zgp(1), fb.F(1), fb.DC(1), fb.k(1), fb.a(1)]; % Parameter lower bounds [v_xgp, v_ygp, v_zgp, F, DC, k, a]
    fit_RCA.ub = [fb.v_xgp(2), fb.v_ygp(2), fb.v_zgp(2), fb.F(2), fb.DC(2), fb.k(2), fb.a(2)]; % Parameter upper bounds
    fit_RCA.model = @(X, cols) vUS_3D_quad_RCA_xy_batchModel(X, tau(cols), PP.k0, sigma_RCA, s, w);
    % vi   = find(maskToUse);  
    nv = numel(maskToUseTrueInds);
    fit_RCA.opts   = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', 2, 'jacobian', 'analytic', 'xscale', fit_RCA.xscale);
    fit_RCA.x0   = [Vt0(maskToUseTrueInds), Vt0(maskToUseTrueInds), Vz0(maskToUseTrueInds), ones(nv, 1), zeros(nv, 1), 2.5*ones(nv,1 ), 0.7 .*ones(nv, 1)];
    opts = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', 2, 'jacobian', 'analytic', 'xscale', fit_RCA.xscale);
    [fit_RCA.x, fit_RCA.cost, ~, fit_RCA.conv] = fitBatchedLM(fit_RCA.model, fit_RCA.x0, g1_exp{j}(maskToUseTrueInds, :), fit_RCA.lb, fit_RCA.ub, fit_RCA.opts);
    [fit_RCA] = storeFitParams(fit_RCA, fpn_RCA, fit_RCA.x, maskToUseTrueInds, PP); % Store/parse fitted parameters
    toc

    % Combined parameter model
    % fit_C.lb = [fb.C(1), fb.v_zgp(1), fb.F(1), fb.DC(1)];
    % fit_C.ub = [fb.C(2), fb.v_zgp(2), fb.F(2), fb.DC(2)];
    % fit_C.x0 = [500*ones(nv, 1), Vz0(maskToUseTrueInds), ones(nv, 1), zeros(nv, 1)];
    % fit_C.xscale = [500, 1e-2, 1, 1];  % C in units of 500 (its typical magnitude), v_zgp x10mm/s, F, DC
    % 
    % fit_C.model = @(X, cols) vUS_3D_combined_batchModel(X, tau(cols), PP.k0);
    % fit_C.opts = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', t1i, 'jacobian', 'analytic', 'xscale', fit_C.xscale);
    % tic
    % [fit_C.x, fit_C.cost, ~, fit_C.conv] = fitBatchedLM(fit_C.model, fit_C.x0, g1_exp{j}(maskToUseTrueInds, :), fit_C.lb, fit_C.ub, fit_C.opts);
    % [fit_C] = storeFitParams(fit_C, fpn_C, fit_C.x, maskToUseTrueInds, PP); % Store/parse fitted parameters
    % toc
   
end

% fit_sfa = createStruct(fit_SG, fit_RCA, fit_C, fpn_SG, fpn_RCA, fpn_C, maskToUse, maskToUseTrueInds, s, w, fb, nv, vs, num_voxels, t1i, tau_decayed_ind, voxel_quality);
fit_sfa = createStruct(fit_SG, fit_RCA, fpn_SG, fpn_RCA, maskToUse, maskToUseTrueInds, s, w, fb, nv, vs, num_voxels, t1i, tau_decayed_ind, voxel_quality);
% save([PDpath, 'fit_sfa.mat'], 'fit_sfa')
% clearvars fit_SG fit_RCA fit_C

%% Calculate the fitted g1 curves for each valid pixel (full model)
g1_model_SG_stacked = zeros(num_voxels, nTau);
g1_model_RCA_stacked = zeros(num_voxels, nTau);
% g1_model_C_stacked = zeros(num_voxels, nTau);

tic
parfor vi = 1:num_voxels % voxel index
    if maskToUse(vi) % If the voxel was fitted
        x_SG = [fit_SG.v_tgp_stacked(vi), fit_SG.v_zgp_stacked(vi), fit_SG.F_stacked(vi), fit_SG.DC_stacked(vi), fit_SG.k_stacked(vi), fit_SG.a_stacked(vi)];
        x_RCA = [fit_RCA.v_xgp_stacked(vi), fit_RCA.v_ygp_stacked(vi), fit_RCA.v_zgp_stacked(vi), fit_RCA.F_stacked(vi), fit_RCA.DC_stacked(vi), fit_RCA.k_stacked(vi), fit_RCA.a_stacked(vi)];
        % x_C = [fit_C.C_stacked(vi), fit_C.v_zgp_stacked(vi), fit_C.F_stacked(vi), fit_C.DC_stacked(vi)];

        % g1_model(vi, :) = vUS_3D_num_wrapper(x, tau, PP.k0, sigma);
        g1_model_SG_stacked(vi, :) = vUS_3D_quad_vec(x_SG, tau, PP.k0, sigma_SG, s, w);
        g1_model_RCA_stacked(vi, :) = vUS_3D_quad_RCA_xy_vec(x_RCA, tau, PP.k0, sigma_RCA, s, w);
        % g1_model_C_stacked(vi, :) = vUS_3D_combined_complex_Jac(x_C, tau, PP.k0);

    end
end

% Calculate fitting quality metrics
tau_mask = (fit_sfa.t1i:nTau);
fit_SG.R2 = calcR2(g1_model_SG_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2);
fit_RCA.R2 = calcR2(g1_model_RCA_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2);
% fit_C.R2 = calcR2(g1_model_C_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2);

fit_SG.R2_adj = calcR2_adj(g1_model_SG_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_SG));
fit_RCA.R2_adj = calcR2_adj(g1_model_RCA_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_RCA));
% fit_C.R2_adj = calcR2_adj(g1_model_C_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_C));

fit_SG.AIC = calcAIC(g1_model_SG_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_SG));
fit_RCA.AIC = calcAIC(g1_model_RCA_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_RCA));
% fit_C.AIC = calcAIC(g1_model_C_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_C));

% Save fit results
% save([PDpath, 'fit_results-', num2str(fi)], 'fit_SG', 'fit_RCA', 'fit_C', 'Vz0', 'Vt0')
% voxelTimeseriesGUI({g1{3}, unstackData(g1_model_RCA_stacked, PP)}, fit_RCA.v, 'ComplexMode', 'abs', 'Colormap', 'jet')
% clearvars fit_SG fit_RCA fit_C g1 g1_exp

toc

%% Fit the simplified version of the nested models
a_fixed = mean([median(fit_RCA.a(maskToUseTrueInds)), median(fit_SG.a(maskToUseTrueInds))])
k_fixed = mean([median(fit_RCA.k(maskToUseTrueInds)), median(fit_SG.k(maskToUseTrueInds))])
% ---- Loop through directional components and go through the fitting process ---- %
% for j = ctp
% for j = 1:2 % Fit only negative and positive frequencies (down and up flows)
for j = 3
    
    % Axial component of the blood flow's group velocity -- v_zgp
    [Vz0, tau_V] = findVzPhaseDiff(stackData(g1{j}, PP), PP); % v_zgp [m/s]
    % figure; imagesc(squeeze(max(unstackData(Vz0, PP), [], 1))'); colormap(VzCmap); axis equal; colorbar; clim([-30e-3, 30e-3])
    % volumeViewer(abs(unstackData(Vz0, PP)))

    % Mesh method for finding v_tgp0
    % [v_zgp0, v_xgp0, p0, DC0, F0, R20] = InitvUS2DParamsWithMesh(g1adj_stacked_j, Vz0, DCR0_j, FR_j, PP, sigma, tau);
    % [Vx0, R2_Vx0] = InitVx0WithMesh2D(stackData(g1{j}, PP), Vz0, DCR0_j, FR0_j, PP, sigma, tau);
    Vt0 = ones(size(Vz0)).* 5e-3; % TESTING: uniform initial v_xgp guess
    % figure; imagesc(unstackData(Vx0, PP)); axis equal; colorbar

    % **** TO DO: create a function that looks at the confidence in the
    % angle-based Vt0, and outputs an updated Vt0 if needed, plus searches
    % for the best p0 ****
    
    % Set parameters for findTauDecayed.m
    fp.absolute_tau_ss_cutoff_s = 10e-3;
    fp.too_fast_decay_s = 2e-3;

    % Adaptively find the tau range to fit over for each pixel, and another
    % quality mask
    [tau_decayed_ind, voxel_quality] = findTauDecayed(g1_exp{j}, tau, t1i, fp.absolute_tau_ss_cutoff_s, fp.too_fast_decay_s);
    % volumeViewer(unstackData(voxel_quality, PP))
    % 
    % % Another mask
    % temp_g1_tau1_mask = squeeze(abs(g1_exp{j}(:, t1i))) > 0.3; % Note: this does not account for the static component
    % % volumeViewer(unstackData(temp_g1_tau1_mask, PP))
    % 
    % temp_g1_tau2_mask = squeeze(abs(g1_exp{j}(:, t1i + 1))) > 0.2; % Note: this does not account for the static component
    % % volumeViewer(unstackData(temp_g1_tau2_mask, PP))
    % 
    % tempMask = and(temp_g1_tau1_mask, temp_g1_tau2_mask);

    % Struct for storing parameters of the different g1 models
    ifpv = zeros(num_voxels, 1); % initial fit parameter value matrix
    fpn_SG = {"v_tgp", "v_zgp", "F", "DC", "k", "a"}; % Fit parameter names (single Gaussian PSF model)
    fit_SG_AF = initFitParamStruct(fpn_SG, ifpv);
    fit_SG_KF = initFitParamStruct(fpn_SG, ifpv);
    fpn_RCA = {"v_xgp", "v_ygp", "v_zgp", "F", "DC", "k", "a"}; % Fit parameter names (RCA-specific PSF model)
    fit_RCA_AF = initFitParamStruct(fpn_RCA, ifpv);
    fit_RCA_KF = initFitParamStruct(fpn_RCA, ifpv);
    % fpn_C = {"C", "v_zgp", "F", "DC"}; % Fit parameter names (Combined parameter model)
    % fit_C = initFitParamStruct(fpn_C, ifpv);
    
    % Choose the mask to use to fit certain pixels or not
    % maskToUse = overall_mask_stacked; 
    % maskToUse = and(vesselAngleMask, stackData(temp_g1_tau1_mask, PP));
    % maskToUse = and(voxel_quality, tempMask);
    maskToUse = or(voxel_quality, tempMask);
    % maskToUse = voxel_quality;
    % maskToUse = vesselMask;
    % volumeViewer(unstackData(maskToUse, PP))

    maskToUseTrueInds = find(maskToUse).'; % indices where maskTouse is true
    
    % ---- Use batched solver ---- %
    % [s, w] = gaussLegendre01(48);
    % % Make vectors for bounds for fit parameters (and store in struct 'fb')
    % fb.v_xgp = [0, 250e-3];
    % fb.v_ygp = [0, 250e-3];
    % fb.v_tgp = sqrt(fb.v_xgp.^2 + fb.v_ygp.^2);
    % fb.v_zgp = [-50e-3; 50e-3];
    % fb.F = [0, 1];
    % fb.DC = [0, 1];
    % fb.k = [2, 3];
    % fb.a = [0, 1];
    % % fb.C = [0, Inf];

    % -------- a fixed -------- %
    % Single-Gaussian PSF model
    tic
    fit_SG_AF.lb = [fb.v_tgp(1), fb.v_zgp(1), fb.F(1), fb.DC(1), fb.k(1), a_fixed]; % Parameter lower bounds [v_tgp, v_zgp, F, DC, k, a]
    fit_SG_AF.ub = [fb.v_tgp(2), fb.v_zgp(2), fb.F(2), fb.DC(2), fb.k(2), a_fixed]; % Parameter upper bounds
    
    fit_SG_AF.xscale = [1e-2 1e-2 1 1 1 1];
    fit_SG_AF.model = @(X, cols) vUS_3D_quad_batchModel(X, tau(cols), PP.k0, sigma_SG, s, w);
    % vi   = find(maskToUse);  
    nv = numel(maskToUseTrueInds);
    fit_SG_AF.opts   = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', 2, 'jacobian', 'analytic', 'xscale', fit_SG_AF.xscale);
    fit_SG_AF.x0   = [Vt0(maskToUseTrueInds), Vz0(maskToUseTrueInds), ones(nv, 1), zeros(nv, 1), 2.5*ones(nv, 1), a_fixed .*ones(nv, 1)];
    fit_SG_AF.opts = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', 2, 'jacobian', 'analytic', 'xscale', fit_SG_AF.xscale);
    [fit_SG_AF.x, fit_SG_AF.cost, ~, fit_SG_AF.conv] = fitBatchedLM(fit_SG_AF.model, fit_SG_AF.x0, g1_exp{j}(maskToUseTrueInds, :), fit_SG_AF.lb, fit_SG_AF.ub, fit_SG_AF.opts);    
    % [fit_SG_AF.x, fit_SG_AF.cost] = vUS_3D_quad_fitBatched(g1_exp{j}(maskToUseTrueInds,:), tau_decayed_ind(maskToUseTrueInds), Vz0(maskToUseTrueInds), tau, PP.k0, sigma_SG, fit_SG_AF.lb, fit_SG_AF.ub);
    [fit_SG_AF] = storeFitParams(fit_SG_AF, fpn_SG, fit_SG_AF.x, maskToUseTrueInds, PP); % Store/parse fitted parameters
    toc

    % RCA-specific PSF model
    tic
    fit_RCA_AF.xscale = [1e-2 1e-2 1e-2 1 1 1 1];
    fit_RCA_AF.lb = [fb.v_xgp(1), fb.v_ygp(1), fb.v_zgp(1), fb.F(1), fb.DC(1), fb.k(1), a_fixed]; % Parameter lower bounds [v_xgp, v_ygp, v_zgp, F, DC, k, a]
    fit_RCA_AF.ub = [fb.v_xgp(2), fb.v_ygp(2), fb.v_zgp(2), fb.F(2), fb.DC(2), fb.k(2), a_fixed]; % Parameter upper bounds
    fit_RCA_AF.model = @(X, cols) vUS_3D_quad_RCA_xy_batchModel(X, tau(cols), PP.k0, sigma_RCA, s, w);
    % vi   = find(maskToUse);  
    nv = numel(maskToUseTrueInds);
    fit_RCA_AF.opts   = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', 2, 'jacobian', 'analytic', 'xscale', fit_RCA_AF.xscale);
    fit_RCA_AF.x0   = [Vt0(maskToUseTrueInds), Vt0(maskToUseTrueInds), Vz0(maskToUseTrueInds), ones(nv, 1), zeros(nv, 1), 2.5*ones(nv,1 ), a_fixed .*ones(nv, 1)];
    [fit_RCA_AF.x, fit_RCA_AF.cost, ~, fit_RCA_AF.conv] = fitBatchedLM(fit_RCA_AF.model, fit_RCA_AF.x0, g1_exp{j}(maskToUseTrueInds, :), fit_RCA_AF.lb, fit_RCA_AF.ub, fit_RCA_AF.opts);
    [fit_RCA_AF] = storeFitParams(fit_RCA_AF, fpn_RCA, fit_RCA_AF.x, maskToUseTrueInds, PP); % Store/parse fitted parameters
    toc

    % -------- k fixed -------- %
    % Single-Gaussian PSF model
    tic
    fit_SG_KF.lb = [fb.v_tgp(1), fb.v_zgp(1), fb.F(1), fb.DC(1), k_fixed, fb.a(1)]; % Parameter lower bounds [v_tgp, v_zgp, F, DC, k, a]
    fit_SG_KF.ub = [fb.v_tgp(2), fb.v_zgp(2), fb.F(2), fb.DC(2), k_fixed, fb.a(2)]; % Parameter upper bounds
    
    fit_SG_KF.xscale = [1e-2 1e-2 1 1 1 1];
    fit_SG_KF.model = @(X, cols) vUS_3D_quad_batchModel(X, tau(cols), PP.k0, sigma_SG, s, w);
    % vi   = find(maskToUse);  
    nv = numel(maskToUseTrueInds);
    fit_SG_KF.opts   = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', 2, 'jacobian', 'analytic', 'xscale', fit_SG_KF.xscale);
    fit_SG_KF.x0   = [Vt0(maskToUseTrueInds), Vz0(maskToUseTrueInds), ones(nv, 1), zeros(nv, 1), k_fixed*ones(nv, 1), 0.7 .*ones(nv, 1)];
    fit_SG_KF.opts = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', 2, 'jacobian', 'analytic', 'xscale', fit_SG_KF.xscale);
    [fit_SG_KF.x, fit_SG_KF.cost, ~, fit_SG_KF.conv] = fitBatchedLM(fit_SG_KF.model, fit_SG_KF.x0, g1_exp{j}(maskToUseTrueInds, :), fit_SG_KF.lb, fit_SG_KF.ub, fit_SG_KF.opts);    
    % [fit_SG_KF.x, fit_SG_KF.cost] = vUS_3D_quad_fitBatched(g1_exp{j}(maskToUseTrueInds,:), tau_decayed_ind(maskToUseTrueInds), Vz0(maskToUseTrueInds), tau, PP.k0, sigma_SG, fit_SG_KF.lb, fit_SG_KF.ub);
    [fit_SG_KF] = storeFitParams(fit_SG_KF, fpn_SG, fit_SG_KF.x, maskToUseTrueInds, PP); % Store/parse fitted parameters
    toc

    % RCA-specific PSF model
    tic
    fit_RCA_KF.xscale = [1e-2 1e-2 1e-2 1 1 1 1];
    fit_RCA_KF.lb = [fb.v_xgp(1), fb.v_ygp(1), fb.v_zgp(1), fb.F(1), fb.DC(1), k_fixed, fb.a(1)]; % Parameter lower bounds [v_xgp, v_ygp, v_zgp, F, DC, k, a]
    fit_RCA_KF.ub = [fb.v_xgp(2), fb.v_ygp(2), fb.v_zgp(2), fb.F(2), fb.DC(2), k_fixed, fb.a(2)]; % Parameter upper bounds
    fit_RCA_KF.model = @(X, cols) vUS_3D_quad_RCA_xy_batchModel(X, tau(cols), PP.k0, sigma_RCA, s, w);
    % vi   = find(maskToUse);  
    nv = numel(maskToUseTrueInds);
    fit_RCA_KF.opts   = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', 2, 'jacobian', 'analytic', 'xscale', fit_RCA_KF.xscale);
    fit_RCA_KF.x0   = [Vt0(maskToUseTrueInds), Vt0(maskToUseTrueInds), Vz0(maskToUseTrueInds), ones(nv, 1), zeros(nv, 1), k_fixed .*ones(nv, 1), 0.7 .*ones(nv, 1)];
    opts = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', 2, 'jacobian', 'analytic', 'xscale', fit_RCA_KF.xscale);
    [fit_RCA_KF.x, fit_RCA_KF.cost, ~, fit_RCA_KF.conv] = fitBatchedLM(fit_RCA_KF.model, fit_RCA_KF.x0, g1_exp{j}(maskToUseTrueInds, :), fit_RCA_KF.lb, fit_RCA_KF.ub, fit_RCA_KF.opts);
    [fit_RCA_KF] = storeFitParams(fit_RCA_KF, fpn_RCA, fit_RCA_KF.x, maskToUseTrueInds, PP); % Store/parse fitted parameters
    toc

    % Combined parameter model
    % fit_C.lb = [fb.C(1), fb.v_zgp(1), fb.F(1), fb.DC(1)];
    % fit_C.ub = [fb.C(2), fb.v_zgp(2), fb.F(2), fb.DC(2)];
    % fit_C.x0 = [500*ones(nv, 1), Vz0(maskToUseTrueInds), ones(nv, 1), zeros(nv, 1)];
    % fit_C.xscale = [500, 1e-2, 1, 1];  % C in units of 500 (its typical magnitude), v_zgp x10mm/s, F, DC
    % 
    % fit_C.model = @(X, cols) vUS_3D_combined_batchModel(X, tau(cols), PP.k0);
    % fit_C.opts = struct('window', tau_decayed_ind(maskToUseTrueInds), 't1i', t1i, 'jacobian', 'analytic', 'xscale', fit_C.xscale);
    % tic
    % [fit_C.x, fit_C.cost, ~, fit_C.conv] = fitBatchedLM(fit_C.model, fit_C.x0, g1_exp{j}(maskToUseTrueInds, :), fit_C.lb, fit_C.ub, fit_C.opts);
    % [fit_C] = storeFitParams(fit_C, fpn_C, fit_C.x, maskToUseTrueInds, PP); % Store/parse fitted parameters
    % toc
   
end

%% Calculate the fitted g1 curves for each valid pixel (full model)
g1_model_SG_AF_stacked = zeros(num_voxels, nTau);
g1_model_RCA_AF_stacked = zeros(num_voxels, nTau);
g1_model_SG_KF_stacked = zeros(num_voxels, nTau);
g1_model_RCA_KF_stacked = zeros(num_voxels, nTau);

parfor vi = 1:num_voxels % voxel index
    if maskToUse(vi) % If the voxel was fitted
        x_SG_AF = [fit_SG_AF.v_tgp_stacked(vi), fit_SG_AF.v_zgp_stacked(vi), fit_SG_AF.F_stacked(vi), fit_SG_AF.DC_stacked(vi), fit_SG_AF.k_stacked(vi), fit_SG_AF.a_stacked(vi)];
        x_RCA_AF = [fit_RCA_AF.v_xgp_stacked(vi), fit_RCA_AF.v_ygp_stacked(vi), fit_RCA_AF.v_zgp_stacked(vi), fit_RCA_AF.F_stacked(vi), fit_RCA_AF.DC_stacked(vi), fit_RCA_AF.k_stacked(vi), fit_RCA_AF.a_stacked(vi)];
        x_SG_KF = [fit_SG_KF.v_tgp_stacked(vi), fit_SG_KF.v_zgp_stacked(vi), fit_SG_KF.F_stacked(vi), fit_SG_KF.DC_stacked(vi), fit_SG_KF.k_stacked(vi), fit_SG_KF.a_stacked(vi)];
        x_RCA_KF = [fit_RCA_KF.v_xgp_stacked(vi), fit_RCA_KF.v_ygp_stacked(vi), fit_RCA_KF.v_zgp_stacked(vi), fit_RCA_KF.F_stacked(vi), fit_RCA_KF.DC_stacked(vi), fit_RCA_KF.k_stacked(vi), fit_RCA_KF.a_stacked(vi)];

        g1_model_SG_AF_stacked(vi, :) = vUS_3D_quad_vec(x_SG_AF, tau, PP.k0, sigma_SG, s, w);
        g1_model_RCA_AF_stacked(vi, :) = vUS_3D_quad_RCA_xy_vec(x_RCA_AF, tau, PP.k0, sigma_RCA, s, w);
        g1_model_SG_KF_stacked(vi, :) = vUS_3D_quad_vec(x_SG_KF, tau, PP.k0, sigma_SG, s, w);
        g1_model_RCA_KF_stacked(vi, :) = vUS_3D_quad_RCA_xy_vec(x_RCA_KF, tau, PP.k0, sigma_RCA, s, w);

    end
end

% Calculate fitting quality metrics
tau_mask = (fit_sfa.t1i:nTau);
fit_SG_AF.R2 = calcR2(g1_model_SG_AF_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2);
fit_RCA_AF.R2 = calcR2(g1_model_RCA_AF_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2);
fit_SG_KF.R2 = calcR2(g1_model_SG_KF_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2);
fit_RCA_KF.R2 = calcR2(g1_model_RCA_KF_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2);

fit_SG_AF.R2_adj = calcR2_adj(g1_model_SG_AF_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_SG));
fit_RCA_AF.R2_adj = calcR2_adj(g1_model_RCA_AF_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_RCA));
fit_SG_KF.R2_adj = calcR2_adj(g1_model_SG_KF_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_SG));
fit_RCA_KF.R2_adj = calcR2_adj(g1_model_RCA_KF_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_RCA));

fit_SG_AF.AIC = calcAIC(g1_model_SG_AF_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_SG));
fit_RCA_AF.AIC = calcAIC(g1_model_RCA_AF_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_RCA));
fit_SG_KF.AIC = calcAIC(g1_model_SG_KF_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_SG));
fit_RCA_KF.AIC = calcAIC(g1_model_RCA_KF_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, numel(fpn_RCA));

% Save fit results
% save([PDpath, 'fit_results-', num2str(fi)], 'fit_SG', 'fit_RCA', 'fit_C', 'Vz0', 'Vt0')
% voxelTimeseriesGUI({g1{3}, unstackData(g1_model_RCA_stacked, PP)}, fit_RCA.v, 'ComplexMode', 'abs', 'Colormap', 'jet')
% clearvars fit_SG fit_RCA fit_C g1 g1_exp

toc


%%
pS_SG = numel(fpn_SG) - 1;
pF_SG = numel(fpn_SG);
pS_RCA = numel(fpn_RCA) - 1;
pF_RCA = numel(fpn_RCA);
% F_a_SG = F_stat(g1_model_SG_AF_stacked(maskToUseTrueInds, tau_mask), g1_model_SG_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, pS_SG, pF_SG);
% F_a_RCA = F_stat(g1_model_RCA_AF_stacked(maskToUseTrueInds, tau_mask), g1_model_RCA_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, pS_RCA, pF_RCA);
% F_k_SG = F_stat(g1_model_SG_KF_stacked(maskToUseTrueInds, tau_mask), g1_model_SG_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, pS_SG, pF_SG);
% F_k_RCA = F_stat(g1_model_RCA_KF_stacked(maskToUseTrueInds, tau_mask), g1_model_RCA_stacked(maskToUseTrueInds, tau_mask), g1_exp{3}(maskToUseTrueInds, tau_mask), 2, pS_RCA, pF_RCA);
F_a_SG = F_stat_vUS(fit_SG_AF, fit_SG, pS_SG, pF_SG, n);
F_a_RCA = F_stat_vUS(fit_RCA_AF, fit_RCA, pS_RCA, pF_RCA, n);
F_k_SG = F_stat_vUS(fit_SG_KF, fit_SG, pS_SG, pF_SG, n);
F_k_RCA = F_stat_vUS(fit_RCA_KF, fit_RCA, pS_RCA, pF_RCA, n);

% p = 0.05;
n = tau_decayed_ind(maskToUseTrueInds) - t1i; % # of sample points for voxels that were fitted
n_fake = size(g1_model_SG_stacked(maskToUseTrueInds, tau_mask), 2); % # of sample points (total)
% Calculate p values
p_a_SG = fcdf(F_a_SG, pF_SG - pS_SG, n - pF_SG, "upper");
p_a_RCA = fcdf(F_a_RCA, pF_RCA - pS_RCA, n - pF_RCA, "upper");
p_k_SG = fcdf(F_k_SG, pF_SG - pS_SG, n - pF_SG, "upper");
p_k_RCA = fcdf(F_k_RCA, pF_RCA - pS_RCA, n - pF_RCA, "upper");




% Plot F-stat distributions against the theoretical one (note: this is
% actually wrong because each voxel's # of sample points (n) is different
figure; histogram(F_a_SG, 'Normalization', 'pdf', 'BinWidth', 0.1)
hold on
plot(0:0.01:100, fpdf(0:0.01:100, pF_SG - pS_SG, n_fake - pF_SG))
hold off

% Plot histograms of the F-stat for the a-fixed model
figure; hold on
fa = 0.4; ea = 0.1; nm = 'probability';
histogram(F_a_SG, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.1, 'Normalization', nm)
histogram(F_a_RCA, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.1, 'Normalization', nm)
hold off
xlabel("F-statistic: fixed a"); ylabel("Probability")
title("F-statistic: fixed a = " + num2str(a_fixed))
xline(0, 'r--', 'LineWidth', 2)
legend('Single Gaussian', 'RCA', 'Location', 'northwest')
% xlim([-2, 1])

% Plot histograms of the F-stat for the k-fixed model
figure; hold on
fa = 0.4; ea = 0.1; nm = 'probability';
histogram(F_k_SG, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.01, 'Normalization', nm)
histogram(F_k_RCA, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.01, 'Normalization', nm)
hold off
xlabel("F-statistic: fixed k"); ylabel("Probability")
title("F-statistic: fixed k = " + num2str(k_fixed))
xline(0, 'r--', 'LineWidth', 2)
legend('Single Gaussian', 'RCA', 'Location', 'northwest')
% xlim([-2, 1])

% Plot histograms of the p values for the a-fixed model
figure; hold on
fa = 0.4; ea = 0.1; nm = 'probability';
histogram(p_a_SG, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.01, 'Normalization', nm)
histogram(p_a_RCA, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.01, 'Normalization', nm)
hold off
xlabel("p-value: fixed a"); ylabel("Probability")
title("p-value: fixed a = " + num2str(a_fixed))
xline(0.05, 'r--', 'LineWidth', 2)
legend('Single Gaussian', 'RCA', 'Location', 'northwest')
% xlim([-2, 1])

% Plot histograms of the p values for the k-fixed model
figure; hold on
fa = 0.4; ea = 0.1; nm = 'probability';
histogram(p_k_SG, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.01, 'Normalization', nm)
histogram(p_k_RCA, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.01, 'Normalization', nm)
hold off
xlabel("p-value: fixed k"); ylabel("Probability")
title("p-value: fixed k = " + num2str(k_fixed))
xline(0.05, 'r--', 'LineWidth', 2)
legend('Single Gaussian', 'RCA', 'Location', 'northwest')
% xlim([-2, 1])

%% Visualize the experimental vs. fitted g1
voxelTimeseriesGUI({g1{3}, g1_model_SG}, fit_SG.v, 'DataNames', {'Data', 'Fit'}, 'ComplexMode', 'abs', 'Colormap', 'turbo')
voxelTimeseriesGUI({g1{3}, g1_model_RCA}, fit_RCA.v, 'DataNames', {'Data', 'Fit'}, 'ComplexMode', 'abs', 'Colormap', 'turbo')
voxelTimeseriesGUI({g1{3}, g1_model_C}, fit_C.C, 'DataNames', {'Data', 'Fit'}, 'ComplexMode', 'abs', 'Colormap', 'turbo')

%% Histograms of R^2 and adjusted R^2
figure; hold on
fa = 0.4; ea = 0.1; nm = 'probability';
histogram(fit_SG.R2, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.1', 'Normalization', nm)
histogram(fit_RCA.R2, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.1', 'Normalization', nm)
histogram(fit_C.R2, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.1', 'Normalization', nm)
hold off
legend('Single Gaussian', 'RCA', 'Combined', 'Location', 'northwest')
xlabel("R^2"); ylabel("Probability")
title("R^2")
xlim([-2, 1])

figure; hold on
fa = 0.4; ea = 0.1;
histogram(fit_SG.R2_adj, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.1', 'Normalization', nm)
histogram(fit_RCA.R2_adj, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.1', 'Normalization', nm)
histogram(fit_C.R2_adj, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', 0.1', 'Normalization', nm)
hold off
legend('Single Gaussian', 'RCA', 'Combined', 'Location', 'northwest')
xlabel("R^2"); ylabel("Probability")
title("Adjusted R^2")
xlim([-2, 1])

%% Histograms of AIC
figure; hold on
fa = 0.4; ea = 0.1;
bw = 1;
histogram(fit_SG.AIC, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw, 'Normalization', nm)
histogram(fit_RCA.AIC, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw, 'Normalization', nm)
histogram(fit_C.AIC, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw, 'Normalization', nm)
hold off
legend('Single Gaussian', 'RCA', 'Combined', 'Location', 'northeast')
xlabel("AIC"); ylabel("Probability")
title("AIC")

dAIC = fit_RCA.AIC - fit_C.AIC;
figure; histogram(dAIC, 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw, 'Normalization', nm)
xlabel("ΔAIC"); ylabel("Probability")
title("ΔAIC = AIC_{RCA} - AIC_{C}")
xline(0, 'r--', 'LineWidth', 2)
sum(dAIC < 0)/sum(maskToUse)

%% Visualize the histograms of fitted parameters
fs = 20;
% v_zgp
figure; hold on
bw_vzgp = 1e-3;
histogram(fit_SG.v_zgp_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_vzgp, 'Normalization', nm)
histogram(fit_RCA.v_zgp_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_vzgp, 'Normalization', nm)
histogram(fit_C.v_zgp_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_vzgp, 'Normalization', nm)
hold off
legend('Single Gaussian', 'RCA', 'Combined', 'Location', 'northwest')
xlabel("v_{zgp}"); ylabel("Probability")
title("v_{zgp}")
fontsize(fs, 'points')
% xlim([-2, 1])

% v_transverse_gp
figure; hold on
bw_vtgp = 1e-3;
histogram(fit_SG.v_tgp_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_vtgp, 'Normalization', nm)
histogram(sqrt(fit_RCA.v_xgp_stacked(maskToUseTrueInds).^2 + fit_RCA.v_ygp_stacked(maskToUseTrueInds).^2), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_vtgp, 'Normalization', nm)
hold off
legend('Single Gaussian', 'RCA', 'Location', 'northeast')
xlabel("v_{tgp}"); ylabel("Probability")
title("v_{tgp}")
fontsize(fs, 'points')
% xlim([-2, 1])

% C
figure; hold on
bw_C = 1e-3;
% histogram(fit_C.C_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_C, 'Normalization', nm)
histogram(fit_C.C_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_C, 'Normalization', nm)
hold off
legend('Combined', 'Location', 'northeast')
xlabel("C"); ylabel("Probability")
title("C")
fontsize(fs, 'points')
% xlim([-2, 1])

% a
figure; hold on
bw_a = 1e-3;
histogram(fit_SG.a_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_a, 'Normalization', nm)
histogram(fit_RCA.a_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_a, 'Normalization', nm)
hold off
legend('Single Gaussian', 'RCA', 'Location', 'northwest')
xlabel("a"); ylabel("Probability")
title("a")
fontsize(fs, 'points')

% k
figure; hold on
bw_k = 1e-2;
histogram(fit_SG.k_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_k, 'Normalization', nm)
histogram(fit_RCA.k_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_k, 'Normalization', nm)
hold off
legend('Single Gaussian', 'RCA', 'Location', 'northeast')
xlabel("k"); ylabel("Probability")
title("k")
fontsize(fs, 'points')

% F
figure; hold on
bw_F = 1e-2;
histogram(fit_SG.F_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_F, 'Normalization', nm)
histogram(fit_RCA.F_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_F, 'Normalization', nm)
histogram(fit_C.F_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_F, 'Normalization', nm)
hold off
legend('Single Gaussian', 'RCA', 'Combined', 'Location', 'northeast')
xlabel("F"); ylabel("Probability")
title("F")
xlim([0, 1])
fontsize(fs, 'points')

% DC
figure; hold on
bw_DC = 1e-2;
histogram(fit_SG.DC_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_DC, 'Normalization', nm)
histogram(fit_RCA.DC_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_DC, 'Normalization', nm)
histogram(fit_C.DC_stacked(maskToUseTrueInds), 'FaceAlpha', fa, 'EdgeAlpha', ea, 'BinWidth', bw_DC, 'Normalization', nm)
hold off
legend('Single Gaussian', 'RCA', 'Combined', 'Location', 'northeast')
xlabel("DC"); ylabel("Probability")
title("DC")
xlim([0, 1])
fontsize(fs, 'points')

%% Helper functions

% Plot the vUS fit at one point (z, x) against the experimental data. Do
% up, down, all flow directions separately
function plotg1pt(z, x, useF, useDC, tau, sigma, k0, g1exp, vUS, p, varargin)
    if nargin > 10
        F = varargin{1};
        if nargin > 11
            DC = varargin{2};
            X = paramsToX(z, x, useF, useDC, vUS, p, F, DC);
        else
            X = paramsToX(z, x, useF, useDC, vUS, p, F);
        end
    end
    
    testg1 = g1vUS2D_vec(X, tau(2:end), sigma, k0, useF, useDC);
    % figure; plot(squeeze(g1exp(z, x, :))); hold on; plot(testg1); hold off
    figure; plot(tau, squeeze(abs(g1exp(z, x, :))), '-x', 'LineWidth', 2); hold on; plot(tau(2:end), abs(testg1), '-o', 'LineWidth', 2); hold off; ylabel('|g1|'); xlabel('Time lag [s]'); legend('Data', 'Fit')
    % testspeed = sqrt(sum(X(1:2).^2))
    
end

% Convert fitted params into a vector X
%   Optional inputs: F, DC
function X = paramsToX(z, x, useF, useDC, vUS, p, varargin)
    if nargin > 6
        F = varargin{1};
        if nargin > 7
            DC = varargin{2};
        end
    end

    if useF
        if useDC
            % error('Have not added this in the code yet')
            X = [squeeze(vUS(z, x, :)); p(z, x); F(z, x); DC(z, x)];
        else
            X = [squeeze(vUS(z, x, :)); p(z, x); F(z, x)];
        end
    else
        X = [squeeze(vUS(z, x, :)); p(z, x)];
    end
end