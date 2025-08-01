%% Description:
%       2D (tl-)fUS processing
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

%% Define some parameters (add this to a prompt later)

parameterPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound', 'Tau 1 index for CBFspeed', 'Tau 2 index for CBFspeed', 'Tau 1 index for CBV'};
parameterDefaults = {'1', '', '20', '500', '2', '10', '2'};
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

%% Create a brain mask
load([IQpath, IQfilenameStructure, num2str(5)])
IQ = squeeze( IData + 1i .* QData );
IQfs = squeeze(mean(IQ, 3));
figure; imagesc(abs(squeeze(IQfs)))
% figure; imagesc(abs(squeeze(max(IQfs, [], 2))'))

% Draw a ROI on the coronal MIP
coronal_roi = images.roi.Freehand;
coronal_roi.draw;
coronal_mask = createMask(coronal_roi);
% coronal_mask = drawfreehand;
%% Mask the IQ
IQm = IQ; % IQ masked
coronal_mask_rep = repmat(~coronal_mask, 1, 1, size(IQ, 3));
IQm(coronal_mask_rep) = 0;
figure; imagesc(abs(squeeze(IQm(:, :, 1))))

%% Main loop with the masking
for filenum = startFile:endFile
% for filenum = 2:endFile
% for filenum = [2]
% for filenum = 1
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    
    IQ = squeeze(IData + 1i .* QData);
    clearvars IData QData

    % Apply the mask
    IQm = IQ; % IQ masked
    IQm(coronal_mask_rep) = 0;
    
    % SVD decluttering
    
    % Determine the optimal SV thresholds with the spatial similarity matrix
    [zp, xp, nf] = size(IQm);
    % PP = reshape(IQm, [zp*xp, nf]);
    % tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
    % [U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
    % disp('Full SVD done')
    % toc
    % 
    % SSM = plotSSM(U, false);
    % [~, a_opt, b_opt] = fitSSM(SSM, false); % Get the optimal singular value thresholds
    % SVs = diag(S);

    % Get the filtered IQ
    % [IQf] = applySVs1D(IQm, PP, SVs, V, sv_threshold_lower, sv_threshold_upper); % with fixed SV thresholds
    % [IQf_opt] = applySVs1D(IQ, PP, SVs, V, a_opt, b_opt); % with optimal SV thresholds
    % disp('SVD filtered images put together')

    [PP, EVs, V_sort] = getSVs1D(IQm);
    % disp('SVs decomposed')
    [IQf] = applySVs1D(IQm, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    % disp('SVD filtered images put together')

%     figure; imagesc(squeeze(abs(IQf(:, :, 1))) .^ 0.5)

    % clearvars IQ

    % Use the IQf with separated negative and positive frequency components
%     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
%     
    numg1pts = 20; % Only calculate the first N points
%     g1_n = g1T(IQf_separated{1}, numg1pts);
% %     [CBFsi_n, CBVi_n] = g1_to_CBi(g1_n, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%     g1_p = g1T(IQf_separated{2}, numg1pts);
%     [CBFsi_p, CBVi_p] = g1_to_CBi(g1_p, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
    g1 = g1T(IQf, numg1pts);
%     g1 = g1T(IQf);
%     [CBFsi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
% %     savefast([savepath, 'fUSdata-', num2str(filenum), '.mat'], g1, CBFi, CBVi);
% 
%     [PDI] = calcPowerDoppler(IQf_separated);
%     [CDI] = calcColorDoppler(IQf_FT_separated, P);

    PDI = sum(abs(IQf) .^ 2, 3) ./ size(IQf, 3);
%     figure; imagesc(squeeze(PDI_test .^ 0.5)); colormap hot

%     save([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI', '-v7.3', '-nocompression');
%     disp("PDI and CDI for file " + num2str(filenum) + " saved" )
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'CBFsi', 'CBVi', 'PDI', 'CDI', '-v7.3', '-nocompression');
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'CBFsi', 'CBVi', 'PDI', 'CDI', 'g1_n', 'g1_p', 'CBFsi_n', 'CBVi_n', 'CBFsi_p', 'CBVi_p',  '-v7.3', '-nocompression');
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'g1_n', 'g1_p', 'PDI', 'CDI', '-v7.3', '-nocompression');
    save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'PDI', '-v7.3', '-nocompression');
%     save([savepath, 'g1-', num2str(filenum), '.mat'], 'g1', 'g1_n', 'g1_p', '-v7.3', '-nocompression');
%     save([savepath, 'g1-', num2str(filenum), '.mat'], 'g1', '-v7.3', '-nocompression');

    disp("fUS result for file " + num2str(filenum) + " saved" )
%     disp("g1 result for file " + num2str(filenum) + " saved" )

    toc
    
end
% savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'tau1_index_CBF', 'tau2_index_CBF', 'tau1_index_CBV');
savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'numg1pts');
% savefast([savepath, 'PDI_CDI_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');



%% Convert g1 into CBV, CBFspeed, etc.

g1_tau1_cutoff = 0.0;
% tau_difference_cutoff = 0.2;

for filenum = startFile:endFile
% for filenum = [288]
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
save([savepath, 'tlfUS_proc_params.mat'], 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff');
% save([savepath, 'tlfUStest_proc_params.mat'], 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff');
figure; imagesc(squeeze(CBVi .^ 0.5)); colormap hot
% figure; imagesc(squeeze(max(CBVi(:, :), [], 2) .^ 0.5)'); colormap hot
vcmap = colormap_ULM;
figure; imagesc(squeeze(CBFsi)); colormap(vcmap)

% generateTiffStack_multi({CBVi .^ 0.7}, [8.8, 8.8, 8], 'hot', 5)

%% Store all the updated CBVi and CBFsi across the experiment into one matrix
load([savepath, 'tlfUSdata-', num2str(1), '.mat'], 'CBFsi', 'CBVi')
CBViallSF = zeros([size(CBVi), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
CBViallSF(:, :, 1) = CBVi;

CBFsiallSF = zeros([size(CBFsi), endFile - startFile + 1]); % Matrix with the CBFsi for every superframe
CBFsiallSF(:, :, 1) = CBFsi;

for filenum = startFile + 1:endFile
    load([savepath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi')
    CBViallSF(:, :, filenum) = CBVi;
    CBFsiallSF(:, :, filenum) = CBFsi;
end

% generateTiffStack_acrossframes(CBViallSF .^ 0.5, [8.8, 8.8, 8], 'hot')
%% Store all the PDI across the experiment into one matrix
% load([savepath, 'PDI_CDI-', num2str(1), '.mat'], 'PDI', 'CDI')
% load([savepath, 'fUSdata-', num2str(1), '.mat'], 'PDI', 'CDI')
load([savepath, 'fUSdata-', num2str(2), '.mat'], 'PDI')
% PDIallSF = cell([length(PDI), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
% PDIallSF = cell([size(PDI)]); 
PDIallSF = zeros([size(PDI), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
% PDIallSF(:,  1) = PDI;
% PDIallSF(:, :, 1) = PDI;
% CDIallSF = cell([size(CDI)]); % Matrix with the CBVi for every superframe
% CDIallSF(:,  1) = CDI;

% for filenum = startFile + 1:endFile
for filenum = startFile:endFile
%     load([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI')
%     load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'CDI')
    load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI')
%     PDI = load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI')
%     CDI = load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'CDI')

%     for i = 1:3
%         PDIallSF{i} = cat(3, PDIallSF{i}, PDI{i});
%         CDIallSF{i} = cat(3, CDIallSF{i}, CDI{i});
%     end

    if iscell(PDI)
        PDIallSF(:, :, filenum) = PDI{3};
    else
        PDIallSF(:, :, filenum) = PDI;
    end
end

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

hold on
plot(sfStarts, ones(length(sfStarts), 1), 'x') % Plot the start of each superframe
hold off

clearvars trial

%% NEW TEST OF UPSAMPLING AND INTERPOLATING EACH TRIAL %% (07/15/2025)

% trial_CBVi_us = cell(size(trial_sf)); % Store each resampled trial individually
% % zeros([size(CBViallSF(:, :, 1)), P.daqrate * P.Mcr_fcp.apis.seq_length_s]);
% 
% % Add the CBVi timepoints we do have to the corresponding time point in the
% % daqrate sampling space
% for trial = 1:length(trial_windows)
% % for trial = 1
%     disp("Resampling trial " + num2str(trial))
%     trial_CBVi_us{trial} = NaN([size(CBViallSF(:, :, 1)), P.daqrate * P.Mcr_fcp.apis.seq_length_s]);
%     temp_indices = sfStarts(trial_sf{trial});
%     temp_indices_shifted = temp_indices - trial_windows{trial}(1) + 1; % Shift the indices so they correspond to a trial start at 1
%     trial_CBVi_us{trial}(:, :, temp_indices_shifted) = CBViallSF(:, :, trial_sf{trial});
% end
% 
% figure; plot(squeeze(trial_CBVi_us{1}(50, 50, :)), 'o-')
% 
% %%%% Resample and interpolate %%%%
% trial_CBVi_usi = cell(size(trial_sf)); % Store each resampled trial individually
% testfactor = 100;
% 
% interp_times = 1:testfactor:P.daqrate * P.Mcr_fcp.apis.seq_length_s; % Time points at which we calculate an interpolated value
% for trial = 1:length(trial_windows)
% % for trial = 1
%     disp("Resampling trial " + num2str(trial))
% %     trial_CBVi_usi{trial} = NaN([size(CBViallSF(:, :, 1)), P.daqrate * P.Mcr_fcp.apis.seq_length_s]);
%     temp_indices = sfStarts(trial_sf{trial});
%     temp_indices_shifted = temp_indices - trial_windows{trial}(1) + 1; % Shift the indices so they correspond to a trial start at 1
%     trial_CBVi_usi{trial} = spline(temp_indices_shifted, CBViallSF(:, :, trial_sf{trial}), interp_times);
% end
% 
% % figure; plot(squeeze(trial_CBVi_usi{1}(50, 50, :)), 'o-')
% 
% % Inspect the interpolation
% figure; plot(squeeze(trial_CBVi_us{1}(50, 50, :)), 'o-')
% hold on
% plot(interp_times, squeeze(trial_CBVi_usi{1}(50, 50, :)), '--')
% hold off

%% Resample the trials for the hemodynamic parameters
interp_factor = 100;
[trial_CBVi_usi] = resampleTrials(CBViallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);
[trial_CBFsi_usi] = resampleTrials(CBFsiallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);
[trial_PDI_usi] = resampleTrials(PDIallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);

% Inspect the interpolation
% figure; plot(interp_times, squeeze(trial_PDI_usi{1}(50, 50, :)), '--')

%% Calculate the relative hemodynamic changes for each trial

% trial_CBVi_usi_baseline = cell(size(trial_sf));
% trial_rCBV_usi = cell(size(trial_sf));
% 
% for trial = 1:length(trial_windows)
%     trial_CBVi_usi_baseline{trial} = mean(trial_CBVi_usi{trial}(:, :, 1 : P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor), 3);
%     trial_rCBV_usi{trial} = (trial_CBVi_usi{trial} - trial_CBVi_usi_baseline{trial}) ./ trial_CBVi_usi_baseline{trial};
% end

[trial_CBVi_usi_baseline, trial_rCBV_usi] = fUS_calc_rHP(trial_CBVi_usi, P, interp_factor);
[trial_CBFsi_usi_baseline, trial_rCBFs_usi] = fUS_calc_rHP(trial_CBFsi_usi, P, interp_factor);
[trial_PDI_usi_baseline, trial_rPDI_usi] = fUS_calc_rHP(trial_PDI_usi, P, interp_factor);

%% Trial average the relative hemodynamic changes

rCBV_TA = fUS_trialAverage(trial_rCBV_usi);
rCBFs_TA = fUS_trialAverage(trial_rCBFs_usi);
rPDI_TA = fUS_trialAverage(trial_rPDI_usi);

%% Correlation on the trial averaged rCBV

% Resample the stim pattern/predicted HRF
trial_stim_pattern = zeros(P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor, 1);
trial_stim_pattern(P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor : ...
    P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor + ...
    P.Mcr_fcp.apis.stim_length_s * P.daqrate / interp_factor) = 1;
figure; plot(trial_stim_pattern); title('Trial stim pattern')

zt = 12;
[r_rCBV, z_rCBV, am_rCBV] = activationMap2D(rCBV_TA, trial_stim_pattern, zt);

figure; imagesc(r_rCBV); colormap jet; colorbar; clim([0, 1])
figure; imagesc(z_rCBV); colormap jet; colorbar
figure; imagesc(am_rCBV); colormap jet; title("Activation Map (rCBV) with z threshold = " + num2str(zt)); colorbar

%% Correlation on the trial averaged rCBFspeed

% Resample the stim pattern/predicted HRF
trial_stim_pattern = zeros(P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor, 1);
trial_stim_pattern(P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor : ...
    P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor + ...
    P.Mcr_fcp.apis.stim_length_s * P.daqrate / interp_factor) = 1;
% figure; plot(trial_stim_pattern); title('Trial stim pattern')

zt = 12;
[r_rCBFs, z_rCBFs, am_rCBFs] = activationMap2D(rCBFs_TA, trial_stim_pattern, zt);

figure; imagesc(r_rCBFs); colormap jet; clim([0, 1]); colorbar
figure; imagesc(z_rCBFs); colormap jet; colorbar
figure; imagesc(am_rCBFs); colormap jet; title("Activation Map (rCBFs) with z threshold = " + num2str(zt)); colorbar

%% Correlation on the trial averaged rPDI

% Resample the stim pattern/predicted HRF
trial_stim_pattern = zeros(P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor, 1);
trial_stim_pattern(P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor : ...
    P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor + ...
    P.Mcr_fcp.apis.stim_length_s * P.daqrate / interp_factor) = 1;
% figure; plot(trial_stim_pattern); title('Trial stim pattern')

zt = 12;
[r_rPDI, z_rPDI, am_rPDI] = activationMap2D(rPDI_TA, trial_stim_pattern, zt);

figure; imagesc(r_rPDI); colormap jet; clim([0, 1]); colorbar
figure; imagesc(z_rPDI); colormap jet; colorbar
figure; imagesc(am_rPDI); colormap jet; title("Activation Map (rPDI) with z threshold = " + num2str(zt)); colorbar

%% Remove points outside of the brain region (manually selected)
figure; imagesc(trial_CBVi_usi_baseline{1} .^ 0.5); % colormap hot % CBVi map
brain_mask = roipoly; % manually define the ROI
figure; imagesc(brain_mask)

am_rCBV_inbrain = am_rCBV;
am_rCBV_inbrain(~brain_mask) = 0;
figure; imagesc(am_rCBV_inbrain); colormap jet; title("Activation Map (rCBV) masked to the brain with z threshold = " + num2str(zt))

%% Look at the timecourse from a ROI
figure; imagesc(am_rCBV_inbrain); colormap jet; title("Activation Map (rCBV) masked to the brain with z threshold = " + num2str(zt))
roi_mask = roipoly; % manually define the ROI
figure; imagesc(roi_mask)

numPtsUSI = P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor; % # of time points per trial for the upsampling
% Calculate the timecourse from the average within that ROI
roi_rCBV_TA = zeros(size(rCBV_TA, 3), 1);
% repmat(roi_mask, [1, 1, stim_pattern.trial_duration])
for ti = 1:numPtsUSI
% for ti = 1
     temp_rCBV_TA = rCBV_TA(:, :, ti);
     temp_roi_rCBV_avg = mean(temp_rCBV_TA(roi_mask));
     roi_rCBV_TA(ti) = temp_roi_rCBV_avg;
end

% Plot the average timecourse in the ROI
% figure; plot(roi_rCBV_TA)
figure; plot(smoothdata(roi_rCBV_TA, 'movmean', 30))

%% Look at the timecourse from a random ROI
figure; imagesc(am_rCBV_inbrain); colormap jet; title("Activation Map (rCBV) masked to the brain with z threshold = " + num2str(zt))
random_roi_mask = roipoly; % manually define the ROI
figure; imagesc(random_roi_mask)

numPtsUSI = P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor; % # of time points per trial for the upsampling
% Calculate the timecourse from the average within that ROI
random_roi_rCBV_TA = zeros(size(rCBV_TA, 3), 1);
% repmat(roi_mask, [1, 1, stim_pattern.trial_duration])
for ti = 1:numPtsUSI
% for ti = 1
     temp_random_rCBV_TA = rCBV_TA(:, :, ti);
     temp_random_roi_rCBV_avg = mean(temp_random_rCBV_TA(random_roi_mask));
     random_roi_rCBV_TA(ti) = temp_random_roi_rCBV_avg;
end

% Plot the average timecourse in the ROI
figure; plot(random_roi_rCBV_TA)
figure; plot(smoothdata(random_roi_rCBV_TA, 'movmean', 30))

%% Helper functions

function [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF)

    g1A_T = {};
    
    g1A_T{1} = abs(g1(:, :, 2)) > g1_tau1_cutoff; % First treatment: tau1 is above some cutoff (make sure there is some actual blood signal there)
    g1A_T{2} = abs(g1(:, :, tau1_index_CBF)) > abs(g1(:, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > |g1(tau2)| (noise might have the g1 randomly increase with tau, but it should not happen with a voxel where there is a real blood signal)
%     g1A_T{3} = abs(g1(:, :, tau1_index_CBF)) > 2 .* abs(g1(:, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    % g1A_T{4} = abs(g1(:, :, :, tau1_index_CBF)) - 1 .* abs(g1(:, :, :, tau2_index_CBF)) > tau_difference_cutoff; % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    
    g1A_mask = true(size(g1A_T{1})); % Mask of voxels to keep for the g1 treatments
    for i = 1:length(g1A_T)
        g1A_mask = and(g1A_mask, g1A_T{i});
    end

end

% Resample trials and interpolate between the hemodynamic data
% function [data_resampled] = resampleTrials(data, trial_sf, trial_windows, sfStarts, P, interp_factor)
% 
%     % Resample and interpolate
%     data_resampled = cell(size(trial_sf)); % Store each resampled trial individually
% %     interp_factor = 100; % Factor by which to "decimate" the daq rate 
% %     for interpolation timepoints
%     
%     interp_times = 1:interp_factor:P.daqrate * P.Mcr_fcp.apis.seq_length_s; % Time points at which we calculate an interpolated value
%     for trial = 1:length(trial_windows)
%         disp("Resampling trial " + num2str(trial))
%         temp_indices = sfStarts(trial_sf{trial});
%         temp_indices_shifted = temp_indices - trial_windows{trial}(1) + 1; % Shift the indices so they correspond to a trial start at 1
%         data_resampled{trial} = spline(temp_indices_shifted, data(:, :, trial_sf{trial}), interp_times);
%     end
% end

% Calculate r(Hemodynamic parameter) -- relative change
% function [data_baseline, data_relative_change] = calculateRelativeChange(data, P, interp_factor)
%     data_baseline = cell(size(data));
%     data_relative_change = cell(size(data));
%     
%     for trial = 1:length(data)
%         data_baseline{trial} = mean(data{trial}(:, :, 1 : P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor), 3);
%         data_relative_change{trial} = (data{trial} - data_baseline{trial}) ./ data_baseline{trial};
%     end
% end
% 
% % Trial average [the relative change in] a hemodynamic parameter (assumed
% % to be a cell array with each cell a separate trial with the same # of sample points)
% function [data_trial_average] = trialAverage(data)
%     data_trial_average = data{1};
%     if length(data) > 1
%         for trial = 2:length(data)
%             data_trial_average = data_trial_average + data{trial};
%         end
%     end
%     data_trial_average = data_trial_average ./ length(data);
% end