%% Assign the superframe trial binning to PDI
% CBViallSFadj = smoothdata(CBViallSF, 4, "sgolay", 9); % SMOOTH THE CBVi
% CBViallSFadj = smoothdata(CBViallSF, 4, "movmean", 3); % SMOOTH THE CBVi
PDIallSFadj = PDIallSF{3}; % all frequencies

trial_PDI = cell(size(trial_sf)); % use the all frequency PDI
minNumPts = Inf;

for trial = 1:length(trial_sf)
    trial_PDI{trial} = PDIallSF{3}(:, :, :, trial_sf{trial});
    minNumPts = min(minNumPts, length(trial_sf{trial})); % Get the minimum number of measurement points across all trials
end

%% Get the mean or max CBVi or PDI etc. within each trial's stimulation period

% Store the actual CBVi or PDI etc. values within the baseline and stim periods
% max_CBVi_stimon = cell(size(trial_sf_stimon));
avg_PDI_stimon = cell(size(trial_sf_stimon));
avg_PDI_baseline = cell(size(trial_sf_baseline));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PDI_relative_change = cell(size(trial_sf)); % Relative change of PDI, per trial, compared to the mean at baseline of that trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for trial = 1:length(trial_windows)
    avg_PDI_baseline{trial} = mean(PDIallSFadj(:, :, :, trial_sf_baseline{trial}), 4);
%     max_CBVi_stimon{trial} = max(CBViallSFsmoothed(:, :, :, trial_sf_stimon{trial}), [], 4);
    avg_PDI_stimon{trial} = mean(PDIallSFadj(:, :, :, trial_sf_stimon{trial}), 4);
end

% Percent change of CBVi and CBFsi for each trial, compared to the mean at baseline
for trial = 1:length(trial_windows)
    PDI_relative_change{trial} = (trial_PDI{trial} - avg_PDI_baseline{trial}) ./ avg_PDI_baseline{trial} .* 100;
end

%% Smoothed percent change of CBVi and CBFsi for each trial, compared to the mean at baseline
PDI_relative_change_smoothed = cell(size(trial_sf)); % Relative change of CBVi, per trial, compared to the mean at baseline of that trial
PDIparam_smoothing_window = 5; % smoothing window size for rCBV and rCBFspeed

for trial = 1:length(trial_windows)

    PDI_relative_change_smoothed{trial} = smoothdata(PDI_relative_change{trial}, 4, "movmean", PDIparam_smoothing_window);
end

%%
% testinvessel = squeeze(CBVi_relative_change{1}(40, 47, 18:38, :));
testinvessel = squeeze(PDI_relative_change_smoothed{1}(40, 47, 18:38, :));
testinvessel_avg = mean(testinvessel, 1);
figure; plot(testinvessel_avg)

% test1 = squeeze(CBVi_relative_change{1}(40, 47, 28, :));
test1 = squeeze(PDI_relative_change_smoothed{1}(40, 47, 28, :));
figure; plot(test1)
test2 = squeeze(PDI_relative_change_smoothed{1}(40, 47, 27, :));
figure; plot(test2)
test3 = squeeze(PDI_relative_change_smoothed{1}(40, 47, 29, :));
figure; plot(test3)

test4 = squeeze(PDI_relative_change_smoothed{1}(40, 21, 28, :));
figure; plot(test4)
test5 = squeeze(PDI_relative_change_smoothed{1}(40, 21, 27, :));
figure; plot(test5)

test6 = squeeze(PDI_relative_change_smoothed{1}(15, 34, 25, :)); % MCA?
figure; plot(test6)

for trial = 1:length(trial_windows)
%     figure; imagesc(squeeze(max(max(CBVi_relative_change{trial}(1:60, :, :, :), [], 4), [], 1) .^ 0.7)'); colormap hot
    figure; imagesc(squeeze(max(max(PDI_relative_change_smoothed{trial}(1:60, :, :, :), [], 4), [], 1) .^ 0.7)'); colormap hot
%     figure; imagesc(squeeze(max(mean(CBVi_relative_change{trial}, 4), [], 1) .^ 0.7)'); colormap hot
%     figure; imagesc(squeeze(max(max(CBFsi_relative_change{trial}, [], 4), [], 1) .^ 1)'); colormap hot
%     figure; imagesc(squeeze(mean(max(CBFsi_relative_change{trial}, [], 4), 1) .^ 2)'); colormap hot
end

%% Do the correlation stuff on PDI
r_PDI_relative_change = [];
z_PDI_relative_change = [];

r_PDI_relative_change_smoothed = [];
z_PDI_relative_change_smoothed = [];

for trial = 1:length(trial_windows)
    [r_PDI_relative_change_smoothed(:, :, :, trial), z_PDI_relative_change_smoothed(:, :, :, trial)] = corrCoef3D(PDI_relative_change_smoothed{trial}, trial_stim_pattern{trial});

end

r_PDI_relative_change_trialavg = mean(r_PDI_relative_change_smoothed, 4);
z_PDI_relative_change_trialavg = mean(r_PDI_relative_change_smoothed, 4);

zscore_mask = z_PDI_relative_change_trialavg < 1;
r_PDI_relative_change_trialavg_thresholded = r_PDI_relative_change_trialavg;
r_PDI_relative_change_trialavg_thresholded(zscore_mask) = 0;