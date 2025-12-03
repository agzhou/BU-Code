%% Description:
%       3D (tl-)fUS processing
%       Timing data should be processed with plotfUStiming.m first

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

savepath = uigetdir('D:\Allen\Data\', 'Select the save path');
savepath = [savepath, '\'];

addpath([cd, '\Speckle tracking']) % add path for the g1 calculation functions

% Load the timing data
[timingFilePathFN, timingFilePath] = uigetfile([IQpath, '..\Timing data\TD.mat'], 'Select the timing data');
timingFilePath = [timingFilePath, timingFilePathFN];
load(timingFilePath)
% load(timingFilePath, 'acqStart', 'airPuffOutput', 'daqStartTimetag', 'sfTimeTags', 'sfTimeTagsDAQStart', 'sfTimeTagsDAQStart_adj', 'sfWidth', 'sfWidth_adj', 'timeStamp')
%% Define some parameters

parameterPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound', 'Tau 1 index for CBFspeed', 'Tau 2 index for CBFspeed', 'Tau 1 index for CBV'};
parameterDefaults = {'1', '', '20', '', '2', '11', '2'};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});
numFiles = endFile - startFile + 1;
sv_threshold_lower = str2double(parameterUserInput{3});
sv_threshold_upper = str2double(parameterUserInput{4});
tau1_index_CBF = str2double(parameterUserInput{5});
tau2_index_CBF = str2double(parameterUserInput{6});
tau1_index_CBV = str2double(parameterUserInput{7});

clearvars parameterPrompt parameterDefaults parameterUserInput

taustep = 1/P.frameRate;
% tau = taustep:taustep:(P.numFramesPerBuffer * taustep);
tau = 0:taustep:((P.numFramesPerBuffer - 1) * taustep);
tau_ms = tau .* 1000; % Assuming even time spacing between frames

% tau1_index_CBF = 2;
% tau2_index_CBF = 6;
% tau1_index_CBV = 2;

%% Set up the High Pass Filter (parameters from the 2020 vUS paper)
HPF.fc = 25; % Cutoff frequency [Hz]
% 25 Hz corresponds to 1 mm/s

HPF.fs = P.frameRate; % Sampling frequency [Hz]
HPF.order = 4; % Butterworth filter order

[HPF.b, HPF.a] = butter(HPF.order, HPF.fc/(HPF.fs/2), 'high');

%% Define the mask manually for now

% load('E:\Allen BME-BOAS-27 Data Backup\AZ01 fUS\07-21-2025 awake RC15gV manual right whisker stim\coronal_mask_rep_07_24_2025.mat')
% load('I:\Ultrasound Data from 04-11-2025 to 05-08-2025\05-06-2025 AZ03 fUS pre-stroke\run 1 all frames stacked\coronal_mask_rep_07_31_2025.mat')
% load('J:\Ultrasound data from 7-21-2025\08-06-2025 AZ01 RCA fUS\coronal_mask_rep.mat')

%% Define other cropping
%     zstart = 40;
% %     zstart = 50;
%     zend = size(IQ, 3);
%     zstart = 15;
    zstart = 45;
    % zstart = 52;
%     zend = 105;
    % zend = 127;
    zend = 130;

%% Save proc params
numg1pts = 10; % Only calculate the first N points of the g1T
save([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'numg1pts', 'zstart', 'zend', 'HPF');

%% Main loop
% for filenum = startFile:endFile
% for filenum = [2:endFile]
% for filenum = 111:endFile
% for filenum = [endFile - 1:-1:startFile]
% for filenum = [116:endFile]
for filenum = 1

    % Load the IQ data
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    
    IQ = single(squeeze(IData + 1i .* QData));
    clearvars IData QData

    % figure; imagesc(squeeze(max(abs(IQ(:, :, :, 2)), [], 1))')
    
    % Crop the IQ first 
    IQm = IQ(:, :, zstart:zend, :);

%     figure; imagesc(squeeze(max(abs(IQm(:, :, :, 2)), [], 1))')

    %%%%%%%%%%%%%% IF USING THE PREDEFINED MASK %%%%%%%%%%%%
%     IQm = IQ;
%     IQm(coronal_mask_rep) = 0; % Apply the brain mask to the IQ: set the non-brain voxels equal to 0


    % SVD decluttering
%     [PP, EVs, V_sort] = getSVs2D(IQm);
    [xp, yp, zp, nf] = size(IQm);
    PP = reshape(IQm, [xp*yp*zp, nf]);
    tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
    [U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
    SVs = diag(S);
%     disp('Full SVD done')
    toc
    disp('SVs decomposed')

    % Plot one SVD subspace as an image
    subspace = 20;
    subspace_img = reshape(U(:, subspace) * SVs(subspace) * V(:, subspace)', [xp, yp, zp, nf]);
    figure; imagesc(squeeze(max(abs(subspace_img(:, :, :, 2)), [], 1))')
%     volumeViewer(abs(subspace_img(:, :, :, 2)))

    SSM = plotSSM(U, false);
%     SSM = plotSSM(U, true);
    [~, a_opt, b_opt] = fitSSM(SSM, false); % Get the optimal singular value thresholds
%     [~, a_opt, b_opt] = fitSSM(SSM, true); % Get the optimal singular value thresholds

    [IQf, noise] = applySVs2D(IQm, PP, SVs, V, sv_threshold_lower, sv_threshold_upper);
%     [IQf, noise] = applySVs2D(IQm, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')

    % Apply the HPF to the post-SVD clutter filtered data
    HPF.dim = length(size(IQf)); % Operate on the time dimension
    IQf_HPF = filter(HPF.b, HPF.a, IQf, [], HPF.dim);

%     volumeViewer(abs(IQf(:, :, :, 1)))
%     figure; imagesc(squeeze(abs(max(IQf(:, :, :, 1), [], 1)))')
    % clearvars IQ

    % Use the IQf with separated negative and positive frequency components
    [IQf_HPF_separated, IQf_HPF_FT_separated] = separatePosNegFreqs(IQf_HPF);

    % Test: plot the frequency spectra for each at some point
    for ind = 1:3
        figure; plot(squeeze(abs(IQf_HPF_FT_separated{ind}(40, 56, 9, :))))
    end

%     g1_n = g1T(IQf_separated{1}, numg1pts);
% %     [CBFsi_n, CBVi_n] = g1_to_CBi(g1_n, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%     g1_p = g1T(IQf_separated{2}, numg1pts);
%     [CBFsi_p, CBVi_p] = g1_to_CBi(g1_p, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)

    % Calculate g1T
    for ind = 1:3
        g1{ind} = g1T(IQf_HPF_separated{ind}, numg1pts);
    end

    % Calculate PDI
    [PDI] = calcPowerDoppler(IQf_HPF_separated, noise);
%     [CDI] = calcColorDoppler(IQf_FT_separated, P);

%     figure; imagesc(squeeze(max(PDI, [], 1))' .^ 0.5); colormap hot
    figure; imagesc(squeeze(max(PDI{1}, [], 1))' .^ 0.5); colormap hot
    figure; imagesc(squeeze(max(PDI{2}, [], 1))' .^ 0.5); colormap hot
    figure; imagesc(squeeze(max(PDI{3}, [], 1))' .^ 0.5); colormap hot
%     figure; imagesc(squeeze(max(PDI ./ noise, [], 1))' .^ 0.5); colormap hot
%     volumeViewer(PDI)

    save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'PDI', 'noise', '-v7.3');
%     save([savepath, 'g1-', num2str(filenum), '.mat'], 'g1', 'g1_n', 'g1_p', '-v7.3', '-nocompression');
%     save([savepath, 'g1-', num2str(filenum), '.mat'], 'g1', '-v7.3', '-nocompression');

    disp("fUS result for file " + num2str(filenum) + " saved" )

    toc
    
end

%% Convert g1 into CBV, CBFspeed, etc.

n_CBV = 2; % n for the new CBV index derivation

g1_tau1_cutoff = 0.2;
% g1_tau1_cutoff = 0.1;

% g1_tau1_cutoff = 0.0;
% tau_difference_cutoff = 0.2;

% for filenum = startFile:endFile
for filenum = 4:endFile
% for filenum = [endFile]
% for filenum = 110
%     load([savepath, 'g1-', num2str(filenum)], 'g1') % Load the saved g1 mat files
    load([savepath, 'fUSdata-', num2str(filenum)], 'g1') % Load the saved g1 mat files

    [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);
%     [g1A_mask] = createg1mask(g1Avg, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);

    [CBFsi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
    % [CBFsi, CBVi] = g1_to_CBi_NEW(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV, n_CBV);

%     CBFsi(~g1A_mask) = -Inf; % Remove noisy points from the CBFspeed index (in theory)
    CBFsi(~g1A_mask) = 0; % Remove noisy points from the CBFspeed index (in theory)

    % save([savepath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
    save([savepath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3');
%     save([savepath, 'tlfUSdatatest-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
    disp("tl-fUS result for file " + num2str(filenum) + " saved" )

end
% save([savepath, 'tlfUS_proc_params.mat'], 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff', 'g1A_mask');
save([savepath, 'tlfUS_proc_params.mat'], 'tau1_index_CBV', 'n_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff', 'g1A_mask');
% save([savepath, 'tlfUStest_proc_params.mat'], 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff');
figure; imagesc(squeeze(max(CBVi(:, :, :), [], 1) .^ 0.3)'); colormap hot
figure; imagesc(squeeze(max(CBVi(:, :, :), [], 3) .^ 0.3)'); colormap hot
vcmap = colormap_ULM;
figure; imagesc(squeeze(mean(CBFsi(:, :, :), 1))'); colormap(vcmap)

% generateTiffStack_multi({CBVi .^ 0.7}, [8.8, 8.8, 8], 'hot', 5)



%% Prepare template(s) for atlas registration
% Create templates for each hemodynamic parameter, averaging across superframes
CBVi_allSF_avg = mean(CBViallSF, 4);
CBFsi_allSF_avg = mean(CBFsiallSF, 4);
PDI_allSF_avg = mean(PDIallSF, 4);

voxel_size = PData.PDelta .* P.wl; % Voxel size (y, x, z) in meters
fUS_volume_dimensions_m = [P.Trans.numelements/2 * P.Trans.spacingMm / 1e3, P.Trans.numelements/2 * P.Trans.spacingMm / 1e3, (P.endDepthMM - P.startDepthMM)/1e3]; % Volume size in meters
fUS_volume_dimensions_voxels = PData.Size; % Volume size in voxels (from the recon PData)

% Adjust the sizes based on the pre-SVD/clutter filtering cropping
fUS_cropped_volume_dimensions_voxels = size(CBVi_allSF_avg);
fUS_cropped_volume_dimensions_m = fUS_cropped_volume_dimensions_voxels ./ fUS_volume_dimensions_voxels .* fUS_volume_dimensions_m;

targetVoxelSizePrompt = {'y Target Voxel Size [um]', 'x Target Voxel Size [um]', 'z Target Voxel Size [um]'};
% targetVoxelSizeDefaults = {'10', '10', '10'};
targetVoxelSizeDefaults = {'50', '50', '50'};
targetVoxelSizeUserInput = inputdlg(targetVoxelSizePrompt, 'Input Target Voxel Size', 1, targetVoxelSizeDefaults);

% Store target voxel size inputs and convert to meters
target_voxel_size(1) = str2double(targetVoxelSizeUserInput{1}) ./ 1e6;
target_voxel_size(2) = str2double(targetVoxelSizeUserInput{2}) ./ 1e6;
target_voxel_size(3) = str2double(targetVoxelSizeUserInput{3}) ./ 1e6;

prereg_interp_factor = voxel_size ./ target_voxel_size;

% Resample hemodynamic parameter template maps to the desired voxel size
CBVi_allSF_avg_rs = imresize3(CBVi_allSF_avg, 'Scale', prereg_interp_factor, 'Method', 'cubic');
CBFsi_allSF_avg_rs = imresize3(CBFsi_allSF_avg, 'Scale', prereg_interp_factor, 'Method', 'cubic');
PDI_allSF_avg_rs = imresize3(PDI_allSF_avg, 'Scale', prereg_interp_factor, 'Method', 'cubic');

% Store pre-registration parameters
prereg_params.orig_voxel_size = voxel_size;
prereg_params.fUS_volume_dimensions_m = fUS_volume_dimensions_m;
prereg_params.fUS_volume_dimensions_voxels = fUS_volume_dimensions_voxels;
prereg_params.fUS_cropped_volume_dimensions_voxels = fUS_cropped_volume_dimensions_voxels;
prereg_params.fUS_cropped_volume_dimensions_m = fUS_cropped_volume_dimensions_m;
prereg_params.target_voxel_size = target_voxel_size;
prereg_params.prereg_interp_factor = prereg_interp_factor;
prereg.params.orig_zstart = zstart;
prereg.params.orig_zend = zend;
% prereg_params. = 

save([savepath, 'fUS_avg_templates.mat'], 'CBVi_allSF_avg_rs', 'CBFsi_allSF_avg_rs', 'PDI_allSF_avg_rs', 'prereg_params')

% Resample activation maps to the desired voxel size
am_rCBV_rs = imresize3(am_rCBV, 'Scale', prereg_interp_factor, 'Method', 'cubic');
am_rCBFspeed_rs = imresize3(am_rCBFspeed, 'Scale', prereg_interp_factor, 'Method', 'cubic');
am_rPDI_rs = imresize3(am_rPDI, 'Scale', prereg_interp_factor, 'Method', 'cubic');

save([savepath, 'fUS_activationmaps.mat'], 'am_rCBV_rs', 'am_rCBFspeed_rs', 'am_rPDI_rs', 'prereg_params')


%% Helper functions

function [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF)

    g1A_T = {};
    
    g1A_T{1} = abs(g1(:, :, :, 2)) > g1_tau1_cutoff; % First treatment: tau1 is above some cutoff (make sure there is some actual blood signal there)
%     g1A_T{2} = abs(g1(:, :, :, tau1_index_CBF)) > abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > |g1(tau2)| (noise might have the g1 randomly increase with tau, but it should not happen with a voxel where there is a real blood signal)
%     g1A_T{3} = abs(g1(:, :, :, tau1_index_CBF)) > 2 .* abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
%     % g1A_T{4} = abs(g1(:, :, :, tau1_index_CBF)) - 1 .* abs(g1(:, :, :, tau2_index_CBF)) > tau_difference_cutoff; % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    
    g1A_mask = true(size(g1A_T{1})); % Mask of voxels to keep for the g1 treatments
    for i = 1:length(g1A_T)
        g1A_mask = and(g1A_mask, g1A_T{i});
    end

end
