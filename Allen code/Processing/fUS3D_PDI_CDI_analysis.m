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

tempwholebrain = 'alive';
testwholebrain = squeeze(PDI_relative_change_smoothed{1}(:, 47, 18:38, :));

% test1 = squeeze(CBVi_relative_change{1}(40, 47, 28, :));
test1 = squeeze(PDI_relative_change_smoothed{1}(40, 47, 28, :));
figure; plot(test1)
test2 = squeeze(PDI_relative_change_smoothed{1}(40, 47, 27, :));
figure; plot(test2)
test3 = squeeze(PDI_relative_change_smoothed{1}(40, 47, 29, :));
figure; plot(test3)


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
%%
% Do the correlation stuff
r_PDI_relative_change = [];
z_PDI_relative_change = [];

r_CDI_relative_change_smoothed = [];
z_CDI_relative_change_smoothed = [];

for trial = 1:length(trial_windows)
%     [r_CBVi_relative_change(:, :, :, trial), z_CBVi_relative_change(:, :, :, trial)] = corrCoef3D(CBVi_relative_change{trial}, trial_stim_pattern{trial});
%     [r_CBFsi_relative_change(:, :, :, trial), z_CBFsi_relative_change(:, :, :, trial)] = corrCoef3D(CBFsi_relative_change{trial}, trial_stim_pattern{trial});
    [r_PDI_relative_change_smoothed(:, :, :, trial), z_PDI_relative_change_smoothed(:, :, :, trial)] = corrCoef3D(PDI_relative_change_smoothed{trial}, trial_stim_pattern{trial});
    [r_CDI_relative_change_smoothed(:, :, :, trial), z_CDI_relative_change_smoothed(:, :, :, trial)] = corrCoef3D(CDI_relative_change_smoothed{trial}, trial_stim_pattern{trial});

end
%%
r_PDI_relative_change_trialavg = mean(r_PDI_relative_change_smoothed, 4);
z_PDI_relative_change_trialavg = mean(z_PDI_relative_change_smoothed, 4);

r_CDI_relative_change_trialavg = mean(r_CDI_relative_change_smoothed, 4);
z_CDI_relative_change_trialavg = mean(z_CDI_relative_change_smoothed, 4);

volumeViewer(r_PDI_relative_change_trialavg)
%%
zscore_mask = z_PDI_relative_change_trialavg < 1;
r_CBVi_relative_change_trialavg_thresholded = r_PDI_relative_change_trialavg;
r_CBVi_relative_change_trialavg_thresholded(zscore_mask) = 0;