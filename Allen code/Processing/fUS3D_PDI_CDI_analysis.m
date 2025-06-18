%% Assign the superframe trial binning to PDI and CDI
% PDIallSFadj = smoothdata(PDIallSF{3}, 4, "sgolay", 9); % SMOOTH THE PDI
PDIallSFadj = PDIallSF{3};
CDIallSFadj = CDIallSF{3};
% PDIallSFadj = smoothdata(PDIallSF{3}, 4, "movmean", 3); % SMOOTH THE PDI

trial_PDI = cell(size(trial_sf)); % use the all frequency PDI
trial_CDI = cell(size(trial_sf));

minNumPts = Inf;
for trial = 1:length(trial_sf)
    trial_PDI{trial} = PDIallSFadj(:, :, :, trial_sf{trial});
    trial_CDI{trial} = CDIallSFadj(:, :, :, trial_sf{trial});
    minNumPts = min(minNumPts, length(trial_sf{trial})); % Get the minimum number of measurement points across all trials
end

%% Get the mean or max CBVi or PDI etc. within each trial's stimulation period

trial_sf_stimon = cell(size(trial_windows));    % Cell array of size (# trials, 1). Each cell contains the superframe indices that correspond to the stimulus period within that trial.
trial_sf_baseline = cell(size(trial_windows));    % Cell array of size (# trials, 1). Each cell contains the superframe indices that correspond to the baseline period within that trial.
stim_length = P.Mcr_fcp.apis.stim_length_s * P.daqrate; % Stim length, adjusted for the DAQ rate

% Get the superframe indices corresponding to the baseline and stim periods
% within each trial
for trial = 1:length(trial_windows)
    trial_sf_baseline{trial} = find(sfStarts >= trial_windows{trial}(1) & sfStarts <= (trial_windows{trial}(1) + stim_prestart_baseline));
    trial_sf_stimon{trial} = find(sfStarts >= (trial_windows{trial}(1) + stim_prestart_baseline) & sfStarts <= (trial_windows{trial}(1) + stim_prestart_baseline + stim_length));
end
clearvars trial

% Get a square wave approximation of when the stim period is, in the
% superframe timing
trial_stim_pattern = cell(size(trial_windows)); % Cell array of size (# trials, 1). Each cell contains a timeseries of the whole trial, with a square wave approximation of the stimulus within that trial.
for trial = 1:length(trial_windows)
    trial_stim_pattern{trial} = zeros(size(trial_sf{trial}));
%     trial_stim_pattern{trial}(stim_starts(trial) : stim_starts(trial) + stim_length) = 1;
    trial_stim_pattern{trial}(find(sfStarts >= (trial_windows{trial}(1) + stim_prestart_baseline) & sfStarts <= (trial_windows{trial}(1) + stim_prestart_baseline + stim_length)) - trial_sf{trial}(1) + 1) = 1;
end

% % Test - plot the square wave stim approximation for every trial
% figure; hold on
% for trial = 1:length(trial_windows)
%     plot(trial_stim_pattern{trial})
% end
% hold off
    
% Store the actual CBVi or PDI etc. values within the baseline and stim periods
avg_PDI_stimon = cell(size(trial_sf_stimon));
avg_PDI_baseline = cell(size(trial_sf_baseline));
avg_CDI_stimon = cell(size(trial_sf_stimon));
avg_CDI_baseline = cell(size(trial_sf_baseline));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PDI_relative_change = cell(size(trial_sf)); % Relative change of PDI, per trial, compared to the mean at baseline of that trial
CDI_relative_change = cell(size(trial_sf)); % Relative change of CDI, per trial, compared to the mean at baseline of that trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for trial = 1:length(trial_windows)
    avg_PDI_baseline{trial} = mean(PDIallSFadj(:, :, :, trial_sf_baseline{trial}), 4);
    avg_PDI_stimon{trial} = mean(PDIallSFadj(:, :, :, trial_sf_stimon{trial}), 4);

    avg_CDI_baseline{trial} = mean(CDIallSFadj(:, :, :, trial_sf_baseline{trial}), 4);
    avg_CDI_stimon{trial} = mean(CDIallSFadj(:, :, :, trial_sf_stimon{trial}), 4);
end

% % Test - plot the avg PDI for baseline and stim periods in trial 1
% figure; imagesc(squeeze(max(avg_PDI_baseline{1}(30:50, :, :), [], 1))' .^ 0.5); colormap hot
% figure; imagesc(squeeze(max(avg_PDI_stimon{1}(30:50, :, :), [], 1))' .^ 0.5); colormap hot

% Percent change of PDI and CDI for each trial, compared to the mean at baseline
for trial = 1:length(trial_windows)
    temp_avg_CDI_baseline_trial = avg_CDI_baseline{trial};
    temp_avg_CDI_baseline_trial(temp_avg_CDI_baseline_trial == 0) = 1;
    PDI_relative_change{trial} = (trial_PDI{trial} - avg_PDI_baseline{trial}) ./ avg_PDI_baseline{trial} .* 100;
    CDI_relative_change{trial} = (trial_CDI{trial} - avg_CDI_baseline{trial}) ./ temp_avg_CDI_baseline_trial .* 100;
end

%% Smoothed percent change of CBVi and CBFsi for each trial, compared to the mean at baseline
PDI_relative_change_smoothed = cell(size(trial_sf)); % Relative change of PDI, per trial, compared to the mean at baseline of that trial
CDI_relative_change_smoothed = cell(size(trial_sf)); % Relative change of CDI, per trial, compared to the mean at baseline of that trial
rCBparam_smoothing_window = 5; % smoothing window size for rCBV and rCBFspeed

for trial = 1:length(trial_windows)
    PDI_relative_change_smoothed{trial} = smoothdata(PDI_relative_change{trial}, 4, "movmean", rCBparam_smoothing_window);
    CDI_relative_change_smoothed{trial} = smoothdata(CDI_relative_change{trial}, 4, "movmean", rCBparam_smoothing_window);
end

%%
figure; imagesc(squeeze(max(trial_PDI{1}(40, :, :, 1), [], 1) .^ 0.5)'); colormap hot
testinvessel = squeeze(PDI_relative_change_smoothed{1}(40, 47, 23:34, :));
testinvessel_avg = mean(testinvessel, 1);
figure; plot(testinvessel_avg)

% tempwholebrain = ;
% testwholebrain = squeeze(PDI_relative_change_smoothed{1}(:, 47, 18:38, :));

% test1 = squeeze(CBVi_relative_change{1}(40, 47, 28, :));
test1 = squeeze(PDI_relative_change_smoothed{1}(40, 47, 28, :));
figure; plot(test1)

test2 = squeeze(PDI_relative_change_smoothed{1}(40, 47, 27, :));
figure; yyaxis left; plot(test2); % ylim([0, 20])
yyaxis right
plot(trial_stim_pattern{1}); % ylim([0, 1])

test3 = squeeze(PDI_relative_change_smoothed{1}(40, 47, 29, :));
figure; plot(test3)


test4 = squeeze(PDI_relative_change_smoothed{1}(40, 36, 52, :));
figure; plot(test4)

for trial = 1:length(trial_windows)
    test = PDI_relative_change_smoothed{trial};
    test(test > 50) = 0; % remove insanely large changes which are probably from artifacts
%     figure; imagesc(squeeze(max(max(CBVi_relative_change{trial}(1:60, :, :, :), [], 4), [], 1) .^ 0.7)'); colormap hot
    % figure; imagesc(squeeze(max(max(PDI_relative_change_smoothed{trial}(1:60, :, :, :), [], 4), [], 1) .^ 0.7)'); colormap hot
    figure; imagesc(squeeze(max(max(test(1:60, :, :, :), [], 4), [], 1) .^ 1)'); colormap hot
%     figure; imagesc(squeeze(max(mean(CBVi_relative_change{trial}, 4), [], 1) .^ 0.7)'); colormap hot
%     figure; imagesc(squeeze(max(max(CBFsi_relative_change{trial}, [], 4), [], 1) .^ 1)'); colormap hot
%     figure; imagesc(squeeze(mean(max(CBFsi_relative_change{trial}, [], 4), 1) .^ 2)'); colormap hot
end

%% Do the correlation stuff
r_PDI_relative_change = [];
z_PDI_relative_change = [];

r_CDI_relative_change_smoothed = [];
z_CDI_relative_change_smoothed = [];

activationMap = [];
zt = 1;

for trial = 1:length(trial_windows)
%     [r_CBVi_relative_change(:, :, :, trial), z_CBVi_relative_change(:, :, :, trial)] = corrCoef3D(CBVi_relative_change{trial}, trial_stim_pattern{trial});
%     [r_CBFsi_relative_change(:, :, :, trial), z_CBFsi_relative_change(:, :, :, trial)] = corrCoef3D(CBFsi_relative_change{trial}, trial_stim_pattern{trial});
    % [r_PDI_relative_change_smoothed(:, :, :, trial), z_PDI_relative_change_smoothed(:, :, :, trial)] = corrCoef3D(PDI_relative_change_smoothed{trial}, trial_stim_pattern{trial});
    % [r_CDI_relative_change_smoothed(:, :, :, trial), z_CDI_relative_change_smoothed(:, :, :, trial)] = corrCoef3D(CDI_relative_change_smoothed{trial}, trial_stim_pattern{trial});
    [r_PDI_relative_change_smoothed(:, :, :, trial), z_PDI_relative_change_smoothed(:, :, :, trial), activationMap(:, :, :, trial)] = activationMap3D(PDI_relative_change_smoothed{trial}, trial_stim_pattern{trial}, zt);
end

%% Plot some slices of the correlation coefficient
% figure; imagesc(squeeze(max(PDIallSF{3}(:, :, :, 1), [], 3)) .^ 0.5); colormap hot
% Trial n
for test_trialnum = 1:length(trial_windows)
    % figure; imagesc(squeeze(max(r_PDI_relative_change_smoothed(30:50, :, :, test_trialnum), [], 1))'); colormap hot; clim([0, 1])
    figure; imagesc(squeeze(max(z_PDI_relative_change_smoothed(30:50, :, :, test_trialnum), [], 1))')
end

%% Correlation of the trial average and stim
% need to make some way for the # of points to be the same..
% [r_PDI_relative_change_trialavg, z_PDI_relative_change_trialavg] = corrCoef3D(PDI_relative_change_smoothed_, trial_stim_pattern{trial});

%% Trial average the correlation stuff
r_PDI_relative_change_trialavg = mean(r_PDI_relative_change_smoothed, 4);
z_PDI_relative_change_trialavg = mean(z_PDI_relative_change_smoothed, 4);

r_CDI_relative_change_trialavg = mean(r_CDI_relative_change_smoothed, 4);
z_CDI_relative_change_trialavg = mean(z_CDI_relative_change_smoothed, 4);

% volumeViewer(r_PDI_relative_change_trialavg)
volumeViewer(z_PDI_relative_change_trialavg)

%%
zscore_mask = z_PDI_relative_change_trialavg < 2.5;
r_CBVi_relative_change_trialavg_thresholded = r_PDI_relative_change_trialavg;
r_CBVi_relative_change_trialavg_thresholded(zscore_mask) = 0;

volumeViewer(r_CBVi_relative_change_trialavg_thresholded)
