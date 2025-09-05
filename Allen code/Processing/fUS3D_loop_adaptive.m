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

parameterPrompt = {'Start file number', 'End file number', 'Tau 1 index for CBFspeed', 'Tau 2 index for CBFspeed', 'Tau 1 index for CBV'};
parameterDefaults = {'1', '', '2', '11', '2'};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});
numFiles = endFile - startFile + 1;
tau1_index_CBF = str2double(parameterUserInput{3});
tau2_index_CBF = str2double(parameterUserInput{4});
tau1_index_CBV = str2double(parameterUserInput{5});

clearvars parameterPrompt parameterDefaults parameterUserInput

taustep = 1/P.frameRate;
% tau = taustep:taustep:(P.numFramesPerBuffer * taustep);
tau = 0:taustep:((P.numFramesPerBuffer - 1) * taustep);
tau_ms = tau .* 1000; % Assuming even time spacing between frames

% tau1_index_CBF = 2;
% tau2_index_CBF = 6;
% tau1_index_CBV = 2;

%% Define the mask manually for now

% load('E:\Allen BME-BOAS-27 Data Backup\AZ01 fUS\07-21-2025 awake RC15gV manual right whisker stim\coronal_mask_rep_07_24_2025.mat')
% load('I:\Ultrasound Data from 04-11-2025 to 05-08-2025\05-06-2025 AZ03 fUS pre-stroke\run 1 all frames stacked\coronal_mask_rep_07_31_2025.mat')
load('J:\Ultrasound data from 7-21-2025\08-06-2025 AZ01 RCA fUS\coronal_mask_rep.mat')
%% Set up the High Pass Filter
fc = 50; % Cutoff frequency [Hz]
fs = P.frameRate; % Sampling frequency [Hz]
HPF_order = 3; % Butterworth filter order

[HPF_b, HPF_a] = butter(HPF_order, fc/(fs/2), 'high');

%% Main loop with the Adaptive SVD Thresholding
for filenum = startFile:endFile
% for filenum = 2:endFile
% for filenum = [285:-1:189]
% for filenum = 1
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    IQ = single(squeeze(IData + 1i .* QData));
    clearvars IData QData

    figure; imagesc(squeeze(max(abs(IQ(:, :, :, 2)), [], 1))')

%     IQm = IQ(:, :, 40:end, :);
    IQm = IQ(:, :, 45:105, :);
%     IQm = IQ(:, :, 50:end, :);
%     IQm = IQ(:, :, 15:105, :);
%     figure; imagesc(squeeze(max(abs(IQm(:, :, :, 2)), [], 1))')

    %%%%%%%%%%%%%% IF USING THE MASK %%%%%%%%%%%%
%     IQm = IQ;
%     IQm(coronal_mask_rep) = 0; % Apply the brain mask to the IQ: set the non-brain voxels equal to 0
    
    % Apply the HPF
%     dim = length(size(IQm)); % Operate on the time dimension
%     IQm_HPF = filter(HPF_b, HPF_a, IQm, [], dim);

    IQm_HPF = IQm;
    % Determine the optimal SV thresholds with the spatial similarity matrix
    [xp, yp, zp, nf] = size(IQm_HPF);
    PP = reshape(IQm_HPF, [xp*yp*zp, nf]);
    tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
    [U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
    SVs = diag(S);
%     disp('Full SVD done')
    toc

    % Plot one SVD subspace as an image
    subspace = 50;
    subspace_img = reshape(U(:, subspace) * SVs(subspace) * V(:, subspace)', [xp, yp, zp, nf]);
    figure; imagesc(squeeze(max(abs(subspace_img(:, :, :, 2)), [], 1))')
    volumeViewer(abs(subspace_img(:, :, :, 2)))

    SSM = plotSSM(U, false);
%     SSM = plotSSM(U, true);
    [~, a_opt, b_opt] = fitSSM(SSM, false); % Get the optimal singular value thresholds
%     [~, a_opt, b_opt] = fitSSM(SSM, true); % Get the optimal singular value thresholds
    

%     [PP, EVs, V_sort] = getSVs2D(IQ);
%     disp('SVs decomposed')
    [IQf_HPF, noise] = applySVs2D(IQm_HPF, PP, SVs, V, a_opt, b_opt);
    [IQf_HPF, noise] = applySVs2D(IQm_HPF, PP, SVs, V, sv_threshold_lower, sv_threshold_upper);
%     disp('SVD filtered images put together')

%     volumeViewer(abs(IQf_HPF(:, :, :, 1)))
%     figure; imagesc(squeeze(abs(max(IQf_HPF(:, :, :, 1), [], 1)))')
    % clearvars IQ

    % Use the IQf with separated negative and positive frequency components
%     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
    
    numg1pts = 20; % Only calculate the first N points
%     g1_n = g1T(IQf_separated{1}, numg1pts);
% %     [CBFsi_n, CBVi_n] = g1_to_CBi(g1_n, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%     g1_p = g1T(IQf_separated{2}, numg1pts);
%     [CBFsi_p, CBVi_p] = g1_to_CBi(g1_p, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
    g1 = g1T(IQf_HPF, numg1pts);

%     [PDI] = calcPowerDoppler(IQf_separated);
    PDI = sum(abs(IQf_HPF) .^ 2, 4) ./ size(IQf_HPF, 4);
%     [CDI] = calcColorDoppler(IQf_FT_separated, P);

%     figure; imagesc(squeeze(max(PDI, [], 1))' .^ 0.5); colormap hot
%     figure; imagesc(squeeze(max(PDI ./ noise, [], 1))' .^ 0.5); colormap hot
%     volumeViewer(PDI)

%     save([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI', '-v7.3', '-nocompression');
%     disp("PDI and CDI for file " + num2str(filenum) + " saved" )
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'CBFsi', 'CBVi', 'PDI', 'CDI', '-v7.3', '-nocompression');
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'CBFsi', 'CBVi', 'PDI', 'CDI', 'g1_n', 'g1_p', 'CBFsi_n', 'CBVi_n', 'CBFsi_p', 'CBVi_p',  '-v7.3', '-nocompression');
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'g1_n', 'g1_p', 'PDI', 'CDI', '-v7.3', '-nocompression');
    save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'PDI', 'noise', 'SVs', 'SSM', 'a_opt', 'b_opt', '-v7.3', '-nocompression');
%     save([savepath, 'g1-', num2str(filenum), '.mat'], 'g1', 'g1_n', 'g1_p', '-v7.3', '-nocompression');
%     save([savepath, 'g1-', num2str(filenum), '.mat'], 'g1', '-v7.3', '-nocompression');

    disp("fUS result for file " + num2str(filenum) + " saved" )
%     disp("g1 result for file " + num2str(filenum) + " saved" )

    toc
    
end
% savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'tau1_index_CBF', 'tau2_index_CBF', 'tau1_index_CBV');
% savefast([savepath, 'fUS_proc_params.mat'], 'coronal_mask_rep', 'tau', 'tau_ms', 'numg1pts');
savefast([savepath, 'fUS_proc_params.mat'], 'tau', 'tau_ms', 'numg1pts', 'zstart', 'zend');
% savefast([savepath, 'PDI_CDI_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');

%% Convert g1 into CBV, CBFspeed, etc.

% g1_tau1_cutoff = 0.0;
g1_tau1_cutoff = 0.2;
% tau_difference_cutoff = 0.2;

for filenum = startFile:endFile
% for filenum = [1]
%     load([savepath, 'g1-', num2str(filenum)], 'g1') % Load the saved g1 mat files
    load([savepath, 'fUSdata-', num2str(filenum)], 'g1') % Load the saved g1 mat files

    [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);
%     [g1A_mask] = createg1mask(g1Avg, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);

    [CBFsi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)

%     CBFsi(~g1A_mask) = -Inf; % Remove noisy points from the CBFspeed index (in theory)
    CBFsi(~g1A_mask) = 0; % Remove noisy points from the CBFspeed index (in theory)

    save([savepath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
%     save([savepath, 'tlfUSdatatest-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
    disp("tl-fUS result for file " + num2str(filenum) + " saved" )

end
save([savepath, 'tlfUS_proc_params.mat'], 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff', 'g1A_mask');
% save([savepath, 'tlfUStest_proc_params.mat'], 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff');
figure; imagesc(squeeze(max(CBVi(:, :, :), [], 1) .^ 0.5)'); colormap hot
figure; imagesc(squeeze(max(CBVi(:, :, :), [], 3) .^ 0.5)'); colormap hot
vcmap = colormap_ULM;
figure; imagesc(squeeze(mean(CBFsi(:, :, :), 1))'); colormap(vcmap)

% generateTiffStack_multi({CBVi .^ 0.7}, [8.8, 8.8, 8], 'hot', 5)

%% Store all the updated CBVi and CBFsi across the experiment into one matrix
load([savepath, 'tlfUSdata-', num2str(1), '.mat'], 'CBFsi', 'CBVi')
CBViallSF = zeros([size(CBVi), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
CBViallSF(:, :, :, 1) = CBVi;

CBFsiallSF = zeros([size(CBFsi), endFile - startFile + 1]); % Matrix with the CBFsi for every superframe
CBFsiallSF(:, :, :, 1) = CBFsi;

for filenum = startFile + 1:endFile
    load([savepath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi')
    CBViallSF(:, :, :, filenum) = CBVi;
    CBFsiallSF(:, :, :, filenum) = CBFsi;
end

%% Store all the PDI across the experiment into one matrix (with separated frequencies)
% % load([savepath, 'PDI_CDI-', num2str(1), '.mat'], 'PDI', 'CDI')
% load([savepath, 'fUSdata-', num2str(1), '.mat'], 'PDI', 'CDI')
% % PDIallSF = cell([length(PDI), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
% PDIallSF = cell([size(PDI)]); % Matrix with the CBVi for every superframe
% % PDIallSF(:,  1) = PDI;
% CDIallSF = cell([size(CDI)]); % Matrix with the CBVi for every superframe
% % CDIallSF(:,  1) = CDI;
% 
% % for filenum = startFile + 1:endFile
% for filenum = startFile:endFile
% %     load([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI')
%     load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'CDI')
% %     PDI = load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI')
% %     CDI = load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'CDI')
% 
%     for i = 1:3
% %         PDIallSF1 = cat(4, PDIallSF1, PDI{1});
% %         CDIallSF1 = cat(4, CDIallSF1, CDI{1});
% %         PDIallSF2 = cat(4, PDIallSF2, PDI{2});
% %         CDIallSF2 = cat(4, CDIallSF2, CDI{2});
% %         PDIallSF3 = cat(4, PDIallSF3, PDI{3});
% %         CDIallSF3 = cat(4, CDIallSF3, CDI{3});
%         PDIallSF{i} = cat(4, PDIallSF{i}, PDI{i});
%         CDIallSF{i} = cat(4, CDIallSF{i}, CDI{i});
%     end
% end

%% Store all the PDI across the experiment into one matrix
% load([savepath, 'PDI_CDI-', num2str(1), '.mat'], 'PDI', 'CDI')
% load([savepath, 'fUSdata-', num2str(1), '.mat'], 'PDI', 'CDI')
load([savepath, 'fUSdata-', num2str(1), '.mat'], 'PDI')
PDIallSF = zeros([size(PDI), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
PDIallSF(:, :, :, 1) = PDI;

% for filenum = startFile + 1:endFile
for filenum = startFile:endFile
%     load([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI')
    load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI')

    PDIallSF(:, :, :, filenum) = PDI;
end

%% Visualize the PDI and CDI across the experiment
% mwr = 30:50; % MIP window range
% mdim = 1; % MIP dimension
% newsize = size(CBViallSF);
% newsize(mdim) = 1; % Set the size of the new variable to 1
% CBViMIPStack = zeros(newsize);
PDIMIPStack = squeeze(max(PDIallSF{1}(30:50, :, :, :), [], 1));
CDIMIPStack = squeeze(max(CDIallSF{1}(30:50, :, :, :), [], 1));

%
generateTiffStack_acrossframes(PDIallSF{3} .^ 0.7, [8.8, 8.8, 8], 'hot', 1:80)

%% Visualize the singular thresholds across the experiment
a_allSF = zeros(numFiles, 1);
b_allSF = zeros(numFiles, 1);

for filenum = startFile:endFile
    load([savepath, 'fUSdata-', num2str(filenum)], 'a_opt', 'b_opt') % Load the saved SV thresholds
    a_allSF(filenum) = a_opt;
    b_allSF(filenum) = b_opt;
end

% Plot the calculated optimal singular thresholds across the experiment
figure;
yyaxis left
plot(a_allSF)
yyaxis right
plot(b_allSF)
title('Optimal lower and upper singular value thresholds')
xlabel('Superframe number')
ylabel('Singular value index')
legend('Lower', 'Upper')

%% Visualize the average CBVi across the volume, for each measurement in the experiment
avg_CBVi_involume_allSF = zeros(numFiles, 1);

for filenum = startFile:endFile
    avg_CBVi_involume_allSF(filenum) = mean(CBViallSF(:, :, :, filenum), 'all', 'omitnan');
end

figure; plot(avg_CBVi_involume_allSF)
figure; yyaxis left; plot(a_allSF); yyaxis right; plot(avg_CBVi_involume_allSF)
figure; scatter(a_allSF, avg_CBVi_involume_allSF); xlabel('Lower SV threshold'); ylabel('Average CBV index across the volume')
%% Visualize the CBVi across the experiment
% mwr = 30:50; % MIP window range
% mdim = 1; % MIP dimension
% newsize = size(CBViallSF);
% newsize(mdim) = 1; % Set the size of the new variable to 1
% CBViMIPStack = zeros(newsize);
% CBViMIPStack = squeeze(max(CBViallSF(30:50, :, :, :), [], 1));

generateTiffStack_acrossframes(CBViallSF .^ 0.7, [8.8, 8.8, 8], 'hot', 1:80)

%% Visualize the CBFsi across the experiment
vcmap = colormap_ULM;
generateTiffStack_acrossframes_MeanIPs(CBFsiallSF .^ 0.1, [8.8, 8.8, 8], vcmap, 1:80)
%% Check different MIPs across superframes
yr = 20:40;
generateTiffStack_acrossframes(CBViallSF .^ 0.7, [8.8, 8.8, 8], 'hot', yr)
% generateTiffStack_acrossframes(CBViallSF .^ 1, [8.8, 8.8, 8], 'hot', yr)


%% Separate each trial
ah = 3; % Approximate a cutoff value for analog high

ind_above_ah = find(TD.airPuffOutput > ah); % Get indices of the air puff output above analog high
ind_shift_below_ah = find(TD.airPuffOutput(ind_above_ah - 1) < ah); % See which indices above analog high have an analog low when shifted by -1 (rising edge)
ind_rising_edge = ind_above_ah(ind_shift_below_ah); % Store the original indices for the rising edges
% hold on
% plot(ind_rising_edge, ones(size(ind_rising_edge)) .* 5, 'o')
% hold off

stim_starts_gap = (P.Mcr_fcp.apis.seq_length_s - P.Mcr_fcp.apis.stim_length_s) * P.daqrate; % How long we expect the stim gap to be between the end of one stim to the start of the next
stim_prestart_baseline = (P.Mcr_fcp.apis.delay_time_ms / 1e3) * P.daqrate; % The duration between the baseline period and the corresponding stim start
stim_starts = ind_rising_edge([true; diff(ind_rising_edge) > stim_starts_gap]); % Add a 1/true at the beginning index for the first stim

% Plot the air puff signal and the calculated start points of each stim period
figure; plot(TD.airPuffOutput)
hold on
plot(stim_starts, ones(size(stim_starts)) .* 5, 'o')
hold off

clearvars ind_above_ah ind_shift_below_ah ind_rising_edge
% figure; plot(TD.sfTimeTagsDAQStart_adj) % plot the time tags for each superframe, adjusted to match the DAQ sampling rate

trial_windows = cell(size(stim_starts)); % Cell array of size (# trials, 1). Each cell contains the time points (according to the DAQ rate) that correspond to that trial.
trial_sf = cell(size(trial_windows));    % Cell array of size (# trials, 1). Each cell contains the superframe indices that started within that trial.

sfStarts = (TD.sfTimeTagsDAQStart_adj - TD.sfWidth_adj); % Adjust the superframe time tags so each index is at the start of the superframe acquisition

% Go through each trial within the run and assign the trial timepoints and the corresponding superframe indices
for trial = 1:length(trial_windows)
    trial_windows{trial} = stim_starts(trial) - stim_prestart_baseline : stim_starts(trial) + stim_starts_gap;

    trial_sf{trial} = find(sfStarts >= trial_windows{trial}(1) & sfStarts <= trial_windows{trial}(end));
end
clearvars trial

%% Resample the trials for the hemodynamic parameters
interp_factor = 100;
[trial_CBVi_usi] = resampleTrials(CBViallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);
[trial_CBFsi_usi] = resampleTrials(CBFsiallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);
[trial_PDI_usi] = resampleTrials(PDIallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);

%% Calculate the relative hemodynamic changes for each trial

[trial_CBVi_usi_baseline, trial_rCBV_usi] = fUS_calc_rHP(trial_CBVi_usi, P, interp_factor);
[trial_CBFsi_usi_baseline, trial_rCBFspeed_usi] = fUS_calc_rHP(trial_CBFsi_usi, P, interp_factor);
[trial_PDI_usi_baseline, trial_rPDI_usi] = fUS_calc_rHP(trial_PDI_usi, P, interp_factor);

%% Trial average the relative hemodynamic changes

rCBV_TA = fUS_trialAverage(trial_rCBV_usi);
rCBFspeed_TA = fUS_trialAverage(trial_rCBFspeed_usi);
rPDI_TA = fUS_trialAverage(trial_rPDI_usi);

%% Correlation on the trial averaged rCBV

% Resample the stim pattern/predicted HRF
trial_stim_pattern = zeros(P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor, 1);
trial_stim_pattern(P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor : ...
    P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor + ...
    P.Mcr_fcp.apis.stim_length_s * P.daqrate / interp_factor) = 1;
figure; plot(trial_stim_pattern); title('Trial stim pattern')

zt = 2;
[r_rCBV, z_rCBV, am_rCBV] = activationMap3D(rCBV_TA, trial_stim_pattern, zt);

% volumeViewer(r_rCBV)
% volumeViewer(z_rCBV)
% volumeViewer(am_rCBV)
figure; imagesc(squeeze(max(r_rCBV(:, :, :), [], 1))'); colorbar; colormap jet; title('Correlation map coronal MIP'); clim([0, 1]) %clim([-1, 1])]
figure; imagesc(squeeze(max(z_rCBV(:, :, :), [], 1))'); colorbar; colormap jet; title('z-score map coronal MIP');
% figure; imagesc(squeeze(mean(z_rCBV(:, :, :), 1))'); colormap jet; clim([0, 1]) % clim([-1, 1])
% figure; imagesc(am_rCBV); colormap jet; title("Activation Map (rCBV) with z threshold = " + num2str(zt))
figure; imagesc(squeeze(max(am_rCBV(:, :, :), [], 1))'); colorbar; colormap jet; title("Activation Map (rCBV) coronal MIP with z threshold = " + num2str(zt))

% generateTiffStack_multi({r_rCBV}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({z_rCBV}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({am_rCBV}, [8.8, 8.8, 8], 'jet', 5)

%% Correlation on the trial averaged rCBFspeed

% Resample the stim pattern/predicted HRF
trial_stim_pattern = zeros(P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor, 1);
trial_stim_pattern(P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor : ...
    P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor + ...
    P.Mcr_fcp.apis.stim_length_s * P.daqrate / interp_factor) = 1;
figure; plot(trial_stim_pattern); title('Trial stim pattern')

zt = 2;
[r_rCBFspeed, z_rCBFspeed, am_rCBFspeed] = activationMap3D_boxfilt(rCBFspeed_TA, trial_stim_pattern, zt);

% volumeViewer(r_rCBFspeed)
% volumeViewer(z_rCBFspeed)
% volumeViewer(am_rCBFspeed)
figure; imagesc(squeeze(max(r_rCBFspeed(:, :, :), [], 1))'); colorbar; colormap jet; title('Correlation map coronal MIP'); clim([0, 1]) %clim([-1, 1])]
figure; imagesc(squeeze(max(z_rCBFspeed(:, :, :), [], 1))'); colorbar; colormap jet; title('z-score map coronal MIP');
% figure; imagesc(squeeze(mean(z_rCBFspeed(:, :, :), 1))'); colormap jet; clim([0, 1]) % clim([-1, 1])
% figure; imagesc(am_rCBFspeed); colormap jet; title("Activation Map (rCBV) with z threshold = " + num2str(zt))
figure; imagesc(squeeze(max(am_rCBFspeed(:, :, :), [], 1))'); colorbar; colormap jet; title("Activation Map (rCBFspeed) coronal MIP with z threshold = " + num2str(zt))

% generateTiffStack_multi({r_rCBFspeed}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({z_rCBFspeed}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({am_rCBFspeed}, [8.8, 8.8, 8], 'jet', 5)

%% Correlation on the trial averaged rPDI

% Resample the stim pattern/predicted HRF
trial_stim_pattern = zeros(P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor, 1);
trial_stim_pattern(P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor : ...
    P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor + ...
    P.Mcr_fcp.apis.stim_length_s * P.daqrate / interp_factor) = 1;
figure; plot(trial_stim_pattern); title('Trial stim pattern')

zt = 2;
[r_rPDI, z_rPDI, am_rPDI] = activationMap3D(rPDI_TA, trial_stim_pattern, zt);

% volumeViewer(r_rPDI)
% volumeViewer(z_rPDI)
% volumeViewer(am_rPDI)
figure; imagesc(squeeze(max(r_rPDI(:, :, :), [], 1))'); colorbar; colormap jet; title('Correlation map coronal MIP'); clim([0, 1]) %clim([-1, 1])]
figure; imagesc(squeeze(max(z_rPDI(:, :, :), [], 1))'); colorbar; colormap jet; title('z-score map coronal MIP');
% figure; imagesc(squeeze(mean(z_rPDI(:, :, :), 1))'); colormap jet; clim([0, 1]) % clim([-1, 1])
% figure; imagesc(am_rPDI); colormap jet; title("Activation Map (rCBV) with z threshold = " + num2str(zt))
figure; imagesc(squeeze(max(am_rPDI(:, :, :), [], 1))'); colorbar; colormap jet; title("Activation Map (rPDI) coronal MIP with z threshold = " + num2str(zt))

% generateTiffStack_multi({r_rPDI}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({z_rPDI}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({am_rPDI}, [8.8, 8.8, 8], 'jet', 5)
generateTiffStack_multi({squeeze(mean(PDIallSF, 4)) .^ 0.5, am_rPDI}, [8.8, 8.8, 8], 'jet', 5)
%% Plot activation at each slice
for slice = 1:10
    my_inds = (slice-1)*5:slice*5;
    my_inds = my_inds+1;
    figure; imagesc(squeeze(mean(test(my_inds, :, :), 1))')
    title(num2str(mean(my_inds)))
end

kernel = ones(3, 3, 3);
kernel(2, 2, 2) = 3;%sum(kernel, 'all');

test_r_CBVi_relative_change_conv = convn(test_r_CBVi_relative_change, kernel, 'same');

%% Helper functions

function [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF)

    g1A_T = {};
    
    g1A_T{1} = abs(g1(:, :, :, 2)) > g1_tau1_cutoff; % First treatment: tau1 is above some cutoff (make sure there is some actual blood signal there)
    g1A_T{2} = abs(g1(:, :, :, tau1_index_CBF)) > abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > |g1(tau2)| (noise might have the g1 randomly increase with tau, but it should not happen with a voxel where there is a real blood signal)
    g1A_T{3} = abs(g1(:, :, :, tau1_index_CBF)) > 2 .* abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    % g1A_T{4} = abs(g1(:, :, :, tau1_index_CBF)) - 1 .* abs(g1(:, :, :, tau2_index_CBF)) > tau_difference_cutoff; % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    
    g1A_mask = true(size(g1A_T{1})); % Mask of voxels to keep for the g1 treatments
    for i = 1:length(g1A_T)
        g1A_mask = and(g1A_mask, g1A_T{i});
    end

end

