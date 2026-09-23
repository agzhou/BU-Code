

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

% Get corresponding times of "stim on"

% Determine which superframes are present while the stimulus is on

