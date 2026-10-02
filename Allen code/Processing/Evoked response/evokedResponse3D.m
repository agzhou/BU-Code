%% Description: processing for 3D fUS data, in an evoked response experimental paradigm

%% Add the Processing folder (and subfolders) to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
AllenProcessingCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing\");
addpath(genpath(AllenProcessingCodePath))

%% Choose processed data directory, params, and RF timetags data file
if ~exist('PDpath', 'var')
    PDpath = uigetdir('D:\Allen\Data\', 'Select the Processed data path');
    PDpath = [PDpath, '\'];
end

% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [PDpath, '..\params.mat']);
    load([params_pathname, params_filename])
end

% Load RF time tags: params.mat
if ~exist('RFtimeTags', 'var')
    % Choose and load the RFTimeTagsPerChannel.mat file (from the acquisition and running readRFTimeTags.m)
    [RFTT_filename, RFTT_pathname, ~] = uigetfile('*.mat', 'Select the RF timetags file', [PDpath, '..\RFTimeTagsPerChannel.mat']);
    load([RFTT_pathname, RFTT_filename])
end

% Load trigger (stimulus) data: triggerData.mat
% if ~exist('inScanData', 'var')
    % Choose and load the RFTimeTagsPerChannel.mat file (from the acquisition and running readRFTimeTags.m)
    [triggerData_filename, triggerData_pathname, ~] = uigetfile('*.mat', 'Select the trigger data file', [PDpath, '..\triggerData.mat']);
    load([triggerData_pathname, triggerData_filename])
% end

clearvars params_filename params_pathname RFTT_filename RFTT_pathname triggerData_filename triggerData_pathname

%% Adjust RF timetags
% Create RF timetag-per-superframe vector and subtract so it starts at 0.
%   RFtimeTags is a (# frames per superframe, # channels, # superframes matrix of RF time tags, in seconds, of the start of each frame's acquisition.
RFTT_start = min(RFtimeTags, [], 'all'); % Get the first timetag value to subtract from the rest, so it starts from 0
channelDim = 2;
RFTT = squeeze(mean(RFtimeTags - RFTT_start, channelDim)); % Also collapse the channel dimension, since time tags should be the same across channels
sfStarts = RFTT(1, :).'; % Superframe start times [s]

% Crop the stimulus (air puff) signal so it's aligned with the start of data
% acquisition
analogHigh = 2; % Trigger value past which we treat the signal as 'on'
indFirstOn = find(inScanData > analogHigh, 1, 'first') - 1; % Index (in DAQ samples) at which the first stim 'on' starts
indSFStart_DAQunits = indFirstOn - P.apis.delay_time_ms/1e3 * P.daqrate; % Index (in DAQ samples) at which the data acquisition started
stim = inScanData(indSFStart_DAQunits:end); % Crop the stim/trigger data to start when the data acquisition started
stimTimestamps = timeStamp(indSFStart_DAQunits:end); stimTimestamps = stimTimestamps - stimTimestamps(1); % Cropped timestamps for the stim/trigger data
stimOnsetTimestamps = ((0:P.numTrials-1).' .* P.apis.seq_length_s) + P.apis.delay_time_ms/1e3; % Time [s] at the stim onset for each trial (assumes repeated trials of the same structure)
stimEndTimestamps = stimOnsetTimestamps + P.apis.stim_length_s; % Time [s] at the stim ends for each trial (assumes repeated trials of the same structure)
clearvars inScanData timeStamp

% Visualize stim and superframe times
figure
plot(stimTimestamps, stim)
maxValue = max(stim); % max value of the stim, for shading
minValue = min(stim);
yshade = [maxValue, minValue, minValue, maxValue]; % Go top right and clockwise for the shading patch vertices
sfWidth = P.numFramesPerBuffer/P.frameRate; % Duration of a superframe [s]
hold on
for sftt = sfStarts.'
    xshade = [sftt + sfWidth, sftt + sfWidth, sftt, sftt];
    patch(xshade, yshade, 'g', 'FaceAlpha', .3) % Plot the shaded region
end
clearvars xshade yshade maxValue minValue

% Determine which superframes are present while the stimulus is on

%% Create a struct for all the relevant timing parameters and save
TD = createStruct(RFTT, sfStarts, stim, stimTimestamps, stimOnsetTimestamps, stimEndTimestamps); % Processing Parameters ======> adjust as needed

if ~exist('TDsavepath', 'var')
    TDsavepath = uigetdir([PDpath, '..\'], 'Select the path to save the timing data in');
    TDsavepath = [TDsavepath, '\'];
end
save([TDsavepath, 'TD.mat'], 'TD')

%% Calculate vUS (g1 fits) for each pre-computed g1 file
% See main_vUS_3D.m

%% Go through pre-computed data files (for all superframes in the experiment) and store

% For now, go through only PDI and CDI
% trialsToUse = [1:10];
% for trialInd = trialsToUse
% 
% end
PDIallSF = cell(3, 1);
CDIallSF = cell(3, 1);
for fi = 1:RFcount

    disp(fi)
    load([PDpath, 'fUSdata-', num2str(fi)], 'PDI', 'CDI')
    for j = 1:3
        PDIallSF{j} = cat(4, PDIallSF{j}, PDI{j});
        CDIallSF{j} = cat(4, CDIallSF{j}, CDI{j});
    end
end
% Average PDI and CDI across superframes
PDIA = cell(3, 1);
CDIA = cell(3, 1);
for j = 1:3
    PDIA{j} = mean(PDIallSF{j}, 4);
    CDIA{j} = mean(CDIallSF{j}, 4);
end
% save([PDpath, 'PDIA_CDIA.mat'], "PDIA", "CDIA")

% Load v from vUS fits
VallSF_SG = cell(3, 1);
VallSF_RCA = cell(3, 1);
CallSF_C = cell(3, 1);
for fi = 1:RFcount

    disp(fi)
    load([PDpath, 'fit_results-', num2str(fi)], 'fit_SG', 'fit_RCA', 'fit_C')
    % for j = 1:3
    for j = 3
        VallSF_SG{j} = cat(4, VallSF_SG{j}, fit_SG.v);
        VallSF_RCA{j} = cat(4, VallSF_RCA{j}, fit_RCA.v);
        CallSF_C{j} = cat(4, CallSF_C{j}, fit_C.C);
    end
end

%%
[evoked, tGrid, info] = periStimulusAverage(PDIallSF, sfStarts, stimOnsetTimestamps);
testStim = zeros(size(tGrid));
testStim(tGrid > 5 & tGrid < 10) = 1;

%% Align the timing of the stim and ultrasound acquisition
% cleanStim = zeros(P.apis.seq_length_s * P.numTrials * P.daqrate, 1);
cleanStim = zeros(size(TD.stimTimestamps)); % Vector with 0s when stim is off, and 1 when stim is on
for ti = 1:P.numTrials
    cleanStim(TD.stimTimestamps > TD.stimOnsetTimestamps(ti) & TD.stimTimestamps < TD.stimEndTimestamps(ti)) = 1;
end

%% Evoked response analysis: peri-stimulus averaging and GLM
% Both functions take the same data types and timing, all in the same time base (t = 0 is the start of acquisition):
%   - superframe start times (TD.sfStarts), and the superframe duration (sfWidth)
%   - the stim: onset times for the peri-stimulus average, and the full-rate stim waveform (cleanStim) for the GLM
% Each data type is a cell (one per j) of (x, y, z, superframe) arrays; empty cells (e.g. vUS for j = 1, 2) are skipped.
sfWidth = P.numFramesPerBuffer/P.frameRate; % Duration of a superframe [s]
delay_s = P.apis.delay_time_ms/1e3; % Time before the first stim onset [s]

dataTypes = struct();
dataTypes.PDI = PDIallSF;
dataTypes.CDI = CDIallSF;
% dataTypes.vUS_SG = VallSF_SG;
% dataTypes.vUS_RCA = VallSF_RCA;
assert(size(PDIallSF{3}, 4) == numel(TD.sfStarts), 'Number of superframes in the data (%d) doesn''t match the RF time tags (%d)', size(PDIallSF{3}, 4), numel(TD.sfStarts))

% Peri-stimulus average of each trial, interpolated onto a common time grid (percent change for PDI, difference from baseline for the signed types)
preTime = min(5, delay_s); % Baseline duration before the stim onset [s] (limited by the time available before the first trial)
postTime = P.apis.seq_length_s - delay_s; % Time after the stim onset [s]; through the end of the trial
[evoked, tGrid, evokedInfo] = periStimulusAverage(dataTypes, TD.sfStarts, TD.stimOnsetTimestamps, ...
    'SFWidth', sfWidth, 'Pre', preTime, 'Post', postTime, ...
    'Normalize', struct('PDI', 'percent', 'default', 'diff'));
% evoked.PDI{3} is (x, y, z, # peri-stimulus time points); evokedInfo.sem has the standard error across trials

% GLM: HRF-convolved stim at the DAQ rate, window-averaged over each superframe. Betas in percent change for PDI
% The acquisition can run past the end of the stim recording (the protocol is over, so the stim is off): pad with 'off' so those superframes are used
nPad = ceil((TD.sfStarts(end) + sfWidth) * P.daqrate) - numel(cleanStim);
cleanStimPadded = [cleanStim; zeros(max(nPad, 0), 1)];
hrfPeak = 1; % HRF peak time [s]. The air puff response is fast: in 09-22-2026 AZ07, activation (z-scores) was strongest for peaks of 0.5-1 s, and much weaker for 3 s. Worth sweeping per dataset
[glmRes, glmInfo] = glmActivationMap(dataTypes, TD.sfStarts, cleanStimPadded, ...
    'Fs', P.daqrate, 'SFWidth', sfWidth, 'HRFPeak', hrfPeak, ...
    'Normalize', struct('PDI', 'percent', 'default', 'none'));
% glmRes.PDI{3}.t, .beta, .z, .p, .q are (x, y, z) maps

% Visualize PDI (j = 3): active voxels from the GLM, and their average peri-stimulus time course
jShow = 3;
activeMask = glmRes.PDI{jShow}.q < 0.05 & glmRes.PDI{jShow}.beta > 0; % FDR-corrected, positive activation
figure
subplot(1, 2, 1)
imagesc(squeeze(max(glmRes.PDI{jShow}.z, [], 1)).'); colorbar; axis image
title('GLM z-score (max projection along dim 1)')
subplot(1, 2, 2)
if any(activeMask, 'all')
    evokedPDI = reshape(evoked.PDI{jShow}, [], numel(tGrid)); % (voxels x time)
    plot(tGrid, mean(evokedPDI(activeMask(:), :), 1, 'omitnan'), 'k', 'LineWidth', 1.5) % Mean across active voxels of the trial-averaged response
    hold on
    xline(0, 'g'); xline(P.apis.stim_length_s, 'r') % Stim onset and offset
    xlabel('Time from stim onset [s]'); ylabel('\DeltaPDI [%]')
    title(sprintf('Active voxels (n = %d)', nnz(activeMask)))
    clearvars evokedPDI
else
    title('No active voxels at q < 0.05')
end

% save([PDpath, 'evokedResponse.mat'], 'evoked', 'tGrid', 'evokedInfo', 'glmRes', 'glmInfo', '-v7.3')

%% TESTING: get vessel angle from superframe-averaged PDI (or CDI?) maps
% Load the averaged PDI and CDI maps
[PDIA_CDIA_filename, PDIA_CDIA_pathname, ~] = uigetfile('*.mat', 'Select the averaged PDI and CDI file');
load([PDIA_CDIA_pathname, PDIA_CDIA_filename])

j = 3;

% USF = 5; % Upsampling factor
% spacing = [1, 1];
% sigmas = [1*USF:1:5*USF];
% tau = 0.5;
% brightondark = true;
% vesselnessThreshold = 0.05;
% minBranchLengthPix = 2*USF;
% minSegLengthPix = 2*USF;
% gamma = 0.5;
% % figure; imagesc(PDIA{j} .^ gamma)
% PDI_US_j = imresize(PDIA{j}, USF, "bilinear");
% % figure; imagesc(PDI_US_j .^ gamma)
% % PDIN = PDI ./ max(PDIN, [], 'all'); % Normalized PDI [0, 1]
% % [vessels, angleMap, vesselness, vesselMask, segLabel, skel] = vesselAngle2D(abs(unstackData(Vz03, PP)) .* 1, sigmas, spacing, tau, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix);
% [vessels, angleMap, vesselness, vesselMask, segLabel, skel] = vesselAngle2D(abs(PDI_US_j .^ gamma), sigmas, spacing, tau, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix);
% figure; imagesc(vesselness); axis square
% figure; imagesc(vesselMask); axis square
% figure; imagesc(skel); axis square
% figure; h = imagesc(angleMap); colormap hsv; colorbar; axis square; xlabel('x [mm]'); ylabel('z [mm]'); title('Vessel angle'); set(h, 'AlphaData', ~isnan(angleMap)) % make pixels transparent if the angle = NaN

spacing = [1, 1, 0.5];
sigmas = [1:0.5:3];
tau_va = 0.5;
brightondark = true;
vesselnessThreshold = 0.05;
minBranchLengthPix = 5;
minSegLengthPix = 3;
% PDIN = PDI ./ max(PDIN, [], 'all'); % Normalized PDI [0, 1]
% [vessels, angleMap, vesselness, vesselMask, segLabel, skel] = vesselAngle2D(abs(unstackData(Vz03, PP)) .* 1, sigmas, spacing, tau_va, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix);
% [vessels, angleMap, vesselness, vesselMask, segLabel, skel] = vesselAngle2D(abs(PDIA{j}(zrange, xrange) .^ gamma), sigmas, spacing, tau_va, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix);
% [vessels, angleMap, vesselness, vesselMask, segLabel, skel] = vesselAngle2D(abs(CDIA{j}(zrange, xrange)), sigmas, spacing, tau_va, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix);
[vessels, dir1, dir2, dir3, vesselness, vesselMask, segLabel, skel] = vesselAngle3D(abs(CDIA{j}(:, :, :)), sigmas, spacing, tau_va, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix);
% test = CDIA{3}; test(test >0 ) = 0; test = abs(test);
% [vessels, angleMap, vesselness, vesselMask, segLabel, skel] = vesselAngle2D(test(zrange, xrange), sigmas, spacing, tau_va, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix);
% figure; imagesc(squeeze(max(vesselness, [], 1))')
% figure; imagesc(squeeze(max(vesselMask, [], 1))')
% figure; imagesc(skel)
% figure; h = imagesc(angleMap); colormap hsv; colorbar; axis equal; xlabel('x'); ylabel('z'); title('Vessel angle'); set(h, 'AlphaData', ~isnan(angleMap)) % make pixels transparent if the angle = NaN
% % vesselAngles = deg2rad(angleMap); % Vessel angles [rad]
% vesselAngles = angleMap; % Vessel angles [deg]
% 
% vesselAngleMask = ~isnan(stackData(vesselAngles, PP)); % Mask to avoid NaN voxels in the vessel angle mask
% 
% figure; imagesc(CDIA{3}); colormap(VzCmap); axis equal; colorbar

save([PDpath, 'vesselMask.mat'], 'vesselMask')
