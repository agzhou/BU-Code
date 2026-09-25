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

%% Go through pre-computed data files (for all superframes in the experiment) and store

% Need to make a g1-to-vUS loop function before passing in things to
% here...

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
save([PDpath, 'PDIA_CDIA.mat'], "PDIA", "CDIA")

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
