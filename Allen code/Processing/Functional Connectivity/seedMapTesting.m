
data = PDIallBlocks_reg - mean(PDIallBlocks_reg, 4);

refVol = mean(PDIallBlocks_reg, 4); % Reference volume: mean over all frames

% figure; imagesc(squeeze(max(refVol, [], 2))' .^ 0.5); colormap gray; axis image
figure; imagesc(squeeze(refVol(100, :, :))' .^ 0.5); colormap gray; axis image

%% Set up the Band Pass Filter
fc = [0.01, 0.1]; % Cutoff frequencies [Hz]
fs = mean(diff(TD.sfTimeTags))/(1-bo); % Sampling frequency [Hz]
BPF_order = 3; % Butterworth filter order

[BPF_b, BPF_a] = butter(BPF_order, fc./(fs/2), 'bandpass');

%% Try Band-Pass filtering
dim = length(size(data)); % Operate on the time dimension

% BPF on the original data
data_BPF = filter(BPF_b, BPF_a, data, [], dim);
% BPF on the "denoised" data
data_denoised_BPF = filter(BPF_b, BPF_a, data_denoised, [], dim);


%% Seed stuff on the original data
% seedCoords = [100, 158, 36];
seedCoords = [100, 149, 41];
seed = squeeze(data(seedCoords(1), seedCoords(2), seedCoords(3), :)); % Seed timecourse
figure; plot(seed)

[r, z] = seedCorrMap(seed, data);

% Plot the seed correlation map
figure; imagesc(squeeze(max(r, [], 1))')

%% Seed stuff on the bandpass-filtered data
seedCoords = [100, 158, 36];
% seedCoords = [100, 149, 41];
seed_BPF = squeeze(data_BPF(seedCoords(1), seedCoords(2), seedCoords(3), :)); % Seed timecourse
figure; plot(seed_BPF)

[r, z] = seedCorrMap(seed_BPF, data_BPF);

% Plot the seed correlation map
figure; imagesc(squeeze(max(r, [], 1))')