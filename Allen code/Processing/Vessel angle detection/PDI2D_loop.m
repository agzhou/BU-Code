%% Description:
%       2D (tl-)fUS processing
%       Timing data should be processed with plotfUStiming.m first

clearvars

%% 1. Load acquisition parameters, timing info, beamforming parameters, etc.
IQpath = uigetdir('G:\', 'Select the IQ data path');
IQpath = [IQpath, '\'];

codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
AllenProcessingCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing");
addpath(AllenProcessingCodePath)

% Load parameters
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
% Create vectors defining x and z coordinates in actual length units
x_mm = (1:PData.Size(2)) .* PData.PDelta(1) .* P.wl .* 1e3; % x [mm]
z_mm = (1:PData.Size(1)) .* PData.PDelta(3) .* P.wl .* 1e3; % z [mm]

IQfilenameStructure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

savepath = uigetdir([IQpath, '\..'], 'Select the save path');
savepath = [savepath, '\'];

% addpath([cd, '\Speckle tracking']) % add path for the g1 calculation functions

% Load or define the timing data
% TDchoices = {'Yes', 'No - Manually Define'};
% TDopts.Default = TDchoices{1}; TDopts.Interpreter = 'none';
% TDanswer = questdlg('Use automatically-acquired stim timing data?', 'Timing data source', TDchoices{1}, TDchoices{2}, TDopts);
% switch TDanswer
%     case TDchoices{1}
%         [timingFilePathFN, timingFilePath] = uigetfile([IQpath, '..\Timing data\TD.mat'], 'Select the timing data');
%         timingFilePath = [timingFilePath, timingFilePathFN];
%         load(timingFilePath)
%     case TDchoices{2}
%         manualTimingPrompt = {'Baseline duration [s]', 'Stim duration [s]', 'Rest duration [s]', 'Offset duration [s]', 'Number of trials'};
%         manualTimingDefaults = {'5', '5', '20', '0', '10'};
%         manualTimingUserInput = inputdlg(manualTimingPrompt, 'Input Parameters', 1, manualTimingDefaults);
% 
%         % Store the user inputs for stim timing parameters into the corresponding variables
%         baseline_duration = str2double(manualTimingUserInput{1});
%         stim_duration = str2double(manualTimingUserInput{2});
%         rest_duration = str2double(manualTimingUserInput{3});
%         offset_duration = str2double(manualTimingUserInput{4});
%         num_trials = str2double(manualTimingUserInput{5});
%         stim_sample_rate = 1000; % [Hz]
% 
%         % trial_stim_pattern = 
% 
% end

% load(timingFilePath, 'acqStart', 'airPuffOutput', 'daqStartTimetag', 'sfTimeTags', 'sfTimeTagsDAQStart', 'sfTimeTagsDAQStart_adj', 'sfWidth', 'sfWidth_adj', 'timeStamp')

%% 2. Define some processing parameters

procParamsPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound'};
procParamsDefaults = {'1', '', '20', num2str(P.numFramesPerBuffer)};
procParamsUserInput = inputdlg(procParamsPrompt, 'Input Parameters', 1, procParamsDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(procParamsUserInput{1});
endFile = str2double(procParamsUserInput{2});
numFiles = endFile - startFile + 1;
sv_threshold_lower = str2double(procParamsUserInput{3});
sv_threshold_upper = str2double(procParamsUserInput{4});

clearvars procParamsPrompt procParamsDefaults procParamsUserInput

%% 3. Loop through files: Clutter Filtering and Power Doppler calculation
% for filenum = startFile:endFile
% for filenum = 2:endFile
% for filenum = startFile
% for filenum = 1
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    
    % IQ = squeeze(IData + 1i .* QData);
    % clearvars IData QData
    
    % SVD decluttering
%     [xp, yp, zp, nf] = size(IQ);
    
    [PP, EVs, V_sort] = getSVs1D(IQ);
    disp('SVs decomposed')
    [IQf, noise] = applySVs1D(IQ, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')

%     figure; imagesc(squeeze(abs(IQf(:, :, 1))) .^ 0.5)

    % clearvars IQ

    % Use the IQf with separated negative and positive frequency components
%     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
%     [PDI] = calcPowerDoppler(IQf_separated);
%     [CDI] = calcColorDoppler(IQf_FT_separated, P);

    PDI = sum(abs(IQf) .^ 2, 3) ./ size(IQf, 3);
    % PDI = sum(abs(IQf) .^ 2, 3) ./ size(IQf, 3) ./ noise;
    % figure; imagesc(x_mm, z_mm, squeeze(PDI .^ 0.5)); colormap hot; colorbar; title('Power Doppler'); xlabel('x [mm]'); ylabel('z [mm]')
    % figure; imagesc(x_mm, z_mm, squeeze(abs(IQ(:, :, 1)))); colorbar; title('IQ'); xlabel('x [mm]'); ylabel('z [mm]')

%     save([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI', '-v7.3', '-nocompression');
%     disp("PDI and CDI for file " + num2str(filenum) + " saved" )
    save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', '-v7.3', '-nocompression');

    disp("fUS result for file " + num2str(filenum) + " saved" )
%     disp("g1 result for file " + num2str(filenum) + " saved" )

    toc
    
end
% savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'tau1_index_CBF', 'tau2_index_CBF', 'tau1_index_CBV');
save([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');
% savefast([savepath, 'PDI_CDI_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');

%% 4. Store all the PDI across the experiment into one matrix
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


%% 5. Separate each trial
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

%% 6. Resample the trials for the hemodynamic parameters
interp_factor = 100;
% interp_factor = 1000;

% [trial_CBVi_usi] = resampleTrials(CBViallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);
% [trial_CBFsi_usi] = resampleTrials(CBFsiallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);
[trial_PDI_usi] = resampleTrials(PDIallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);

% Inspect the interpolation
% figure; plot(interp_times, squeeze(trial_PDI_usi{1}(50, 50, :)), '--')

%% 7. Calculate the relative hemodynamic changes for each trial
[trial_PDI_usi_baseline, trial_rPDI_usi] = fUS_calc_rHP(trial_PDI_usi, P, interp_factor);

%% 8. Trial average the relative hemodynamic changes
rPDI_TA = fUS_trialAverage(trial_rPDI_usi);

%% 9. Correlation on the trial averaged rPDI

% Resample the stim pattern/predicted HRF
trial_stim_pattern = zeros(P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor, 1);
trial_stim_pattern(P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor : ...
    P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor + ...
    P.Mcr_fcp.apis.stim_length_s * P.daqrate / interp_factor) = 1;
% figure; plot(trial_stim_pattern); title('Trial stim pattern')

zt = 2;
[r_rPDI, z_rPDI, am_rPDI] = activationMap2D(rPDI_TA, trial_stim_pattern, zt);
% figure; imagesc(x_mm, z_mm, squeeze(PDI .^ 0.5)); colormap hot; colorbar; title('Power Doppler'); xlabel('x [mm]'); ylabel('z [mm]')
figure; imagesc(x_mm, z_mm, r_rPDI); colormap jet; clim([0, 1]); colorbar; title('Correlation map'); xlabel('x [mm]'); ylabel('z [mm]')
% figure; imagesc(x_mm, z_mm, z_rPDI); colormap jet; colorbar; title('z-score map'); xlabel('x [mm]'); ylabel('z [mm]')
% figure; imagesc(x_mm, z_mm, am_rPDI); colormap jet; title("Activation Map (rPDI) with z threshold = " + num2str(zt)); colorbar; xlabel('x [mm]'); ylabel('z [mm]')

%% Remove points outside of the brain region (manually selected)
% figure; imagesc(trial_CBVi_usi_baseline{1} .^ 0.5); % colormap hot % CBVi map
% brain_mask = roipoly; % manually define the ROI
% figure; imagesc(brain_mask)
% 
% am_rCBV_inbrain = am_rCBV;
% am_rCBV_inbrain(~brain_mask) = 0;
% figure; imagesc(am_rCBV_inbrain); colormap jet; title("Activation Map (rCBV) masked to the brain with z threshold = " + num2str(zt))
% 
% %% Look at the timecourse from a ROI (rCBV)
% figure; imagesc(am_rCBV_inbrain); colormap jet; title("Activation Map (rCBV) masked to the brain with z threshold = " + num2str(zt))
% roi_mask = roipoly; % manually define the ROI
% figure; imagesc(roi_mask)
% 
% numPtsUSI = P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor; % # of time points per trial for the upsampling
% % Calculate the timecourse from the average within that ROI
% roi_rCBV_TA = zeros(size(rCBV_TA, 3), 1);
% % repmat(roi_mask, [1, 1, stim_pattern.trial_duration])
% for ti = 1:numPtsUSI
% % for ti = 1
%      temp_rCBV_TA = rCBV_TA(:, :, ti);
%      temp_roi_rCBV_avg = mean(temp_rCBV_TA(roi_mask));
%      roi_rCBV_TA(ti) = temp_roi_rCBV_avg;
% end
% 
% % Plot the average timecourse in the ROI
% figure; plot((1:length(roi_rCBV_TA)) .* interp_factor ./ P.daqrate, roi_rCBV_TA)
% figure; plot((1:length(roi_rCBV_TA)) .* interp_factor ./ P.daqrate, smoothdata(roi_rCBV_TA, 'movmean', 30))
% 
% %% Look at the timecourse from a ROI (rCBFspeed)
% figure; imagesc(am_rCBFs); colormap jet; title("Activation Map (rCBFspeed) with z threshold = " + num2str(zt))
% CBFs_roi_mask = roipoly; % manually define the ROI
% figure; imagesc(CBFs_roi_mask)
% 
% numPtsUSI = P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor; % # of time points per trial for the upsampling
% % Calculate the timecourse from the average within that ROI
% roi_rCBFs_TA = zeros(size(rCBFs_TA, 3), 1);
% % repmat(roi_mask, [1, 1, stim_pattern.trial_duration])
% for ti = 1:numPtsUSI
% % for ti = 1
%      temp_rCBFs_TA = rCBFs_TA(:, :, ti);
%      temp_roi_rCBFs_avg = mean(temp_rCBFs_TA(CBFs_roi_mask));
%      roi_rCBFs_TA(ti) = temp_roi_rCBFs_avg;
% end
% 
% % Plot the average timecourse in the ROI
% figure; plot((1:length(roi_rCBFs_TA)) .* interp_factor ./ P.daqrate, roi_rCBFs_TA)
% figure; plot((1:length(roi_rCBFs_TA)) .* interp_factor ./ P.daqrate, smoothdata(roi_rCBFs_TA, 'movmean', 30))


%% 10. Look at the timecourse from a ROI (rPDI)
figure; imagesc(am_rPDI); colormap jet; title("Activation Map (rPDI) masked to the brain with z threshold = " + num2str(zt))
roi_mask = roipoly; % manually define the ROI
figure; imagesc(roi_mask)

numPtsUSI = P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor; % # of time points per trial for the upsampling
% Calculate the timecourse from the average within that ROI
roi_rPDI_TA = zeros(size(rPDI_TA, 3), 1);
% repmat(roi_mask, [1, 1, stim_pattern.trial_duration])
for ti = 1:numPtsUSI
% for ti = 1
     temp_rPDI_TA = rPDI_TA(:, :, ti);
     temp_roi_rPDI_avg = mean(temp_rPDI_TA(roi_mask));
     roi_rPDI_TA(ti) = temp_roi_rPDI_avg;
end

% Plot the average timecourse in the ROI
figure; plot((1:length(roi_rPDI_TA)) .* interp_factor ./ P.daqrate, roi_rPDI_TA, 'LineWidth', 2); xlabel('Time [s]'); ylabel('rPDI'); title("rPDI ROI timecourse")
mmws = 30; % Movmean window size (in units of the trial interpolation rate)
figure; plot((1:length(roi_rPDI_TA)) .* interp_factor ./ P.daqrate, smoothdata(roi_rPDI_TA, 'movmean', mmws), 'LineWidth', 2); xlabel('Time [s]'); ylabel('rPDI'); title("rPDI ROI timecourse, moving mean over " + num2str(mmws/length(roi_rPDI_TA) * P.Mcr_fcp.apis.seq_length_s) + "s")

%% Look at the timecourse from a random ROI
% figure; imagesc(am_rCBV_inbrain); colormap jet; title("Activation Map (rCBV) masked to the brain with z threshold = " + num2str(zt))
% random_roi_mask = roipoly; % manually define the ROI
% figure; imagesc(random_roi_mask)
% 
% numPtsUSI = P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor; % # of time points per trial for the upsampling
% % Calculate the timecourse from the average within that ROI
% random_roi_rCBV_TA = zeros(size(rCBV_TA, 3), 1);
% % repmat(roi_mask, [1, 1, stim_pattern.trial_duration])
% for ti = 1:numPtsUSI
% % for ti = 1
%      temp_random_rCBV_TA = rCBV_TA(:, :, ti);
%      temp_random_roi_rCBV_avg = mean(temp_random_rCBV_TA(random_roi_mask));
%      random_roi_rCBV_TA(ti) = temp_random_roi_rCBV_avg;
% end
% 
% % Plot the average timecourse in the ROI
% figure; plot(random_roi_rCBV_TA)
% figure; plot(smoothdata(random_roi_rCBV_TA, 'movmean', 30))