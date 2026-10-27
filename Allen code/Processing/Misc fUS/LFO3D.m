%% Description: Look for low frequency oscillations in 3D fUS data

%% Add the Processing folder (and subfolders) to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
AllenProcessingCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing\");
addpath(genpath(AllenProcessingCodePath))

%% Choose processed data directory, params, and RF timetags data file
if ~exist('PDpath', 'var')
    PDpath = uigetdir('D:\Allen\Data\', 'Select the Processed data (PDI, CDI) path');
    PDpath = [PDpath, '\'];
end

if ~exist('FRpath', 'var')
    FRpath = uigetdir(PDpath, 'Select the vUS fit results data path');
    FRpath = [FRpath, '\'];
end

% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [PDpath, '..\params.mat']);
    load([params_pathname, params_filename])
end

% Load RF time tags: RFtimeTags.mat
if ~exist('RFtimeTags', 'var')
    % Choose and load the RFTimeTagsPerChannel.mat file (from the acquisition and running readRFTimeTags.m)
    [RFTT_filename, RFTT_pathname, ~] = uigetfile('*.mat', 'Select the RF timetags file', [PDpath, '..\RFTimeTagsPerChannel.mat']);
    load([RFTT_pathname, RFTT_filename])
end

% % Load trigger (stimulus) data: triggerData.mat
% % if ~exist('inScanData', 'var')
%     % Choose and load the RFTimeTagsPerChannel.mat file (from the acquisition and running readRFTimeTags.m)
%     [triggerData_filename, triggerData_pathname, ~] = uigetfile('*.mat', 'Select the trigger data file', [PDpath, '..\triggerData.mat']);
%     load([triggerData_pathname, triggerData_filename])
% % end

clearvars params_filename params_pathname RFTT_filename RFTT_pathname triggerData_filename triggerData_pathname

if ~exist('PData', 'var')
    % Choose and load the RFTimeTagsPerChannel.mat file (from the acquisition and running readRFTimeTags.m)
    [PData_filename, PData_pathname, ~] = uigetfile('*.mat', 'Select the PData file', [PDpath, '..\PData.mat']);
    load([PData_pathname, PData_filename])
end

load([PDpath, 'fUS_proc_params.mat'])
[x_mm, y_mm, z_mm] = getReconCoords3D(PData, P);
if ~exist('voxelRange', 'var')
    voxelRange = Mcr_voxelRange;
end
x_mm = x_mm(voxelRange{1});
y_mm = y_mm(voxelRange{2});
z_mm = z_mm(voxelRange{3});

%% Adjust RF timetags
% Create RF timetag-per-superframe vector and subtract so it starts at 0.
%   RFtimeTags is a (# frames per superframe, # channels, # superframes matrix of RF time tags, in seconds, of the start of each frame's acquisition.
RFTT_start = min(RFtimeTags, [], 'all'); % Get the first timetag value to subtract from the rest, so it starts from 0
channelDim = 2;
RFTT = squeeze(mean(RFtimeTags - RFTT_start, channelDim)); % Also collapse the channel dimension, since time tags should be the same across channels
sfStarts = RFTT(1, :).'; % Superframe start times [s]



% Visualize stim and superframe times
figure
% plot(stimTimestamps, stim)
maxValue = 1; % max value of the stim, for shading
minValue = 0;
yshade = [maxValue, minValue, minValue, maxValue]; % Go top right and clockwise for the shading patch vertices
frameRate = checkFrameRate(P);
sfWidth = P.numFramesPerBuffer/frameRate; % Duration of a superframe [s]
hold on
for sftt = sfStarts.'
    xshade = [sftt + sfWidth, sftt + sfWidth, sftt, sftt];
    patch(xshade, yshade, 'g', 'FaceAlpha', .3) % Plot the shaded region
end
clearvars xshade yshade maxValue minValue

% Get superframe center timestamps
sfCenters = sfStarts + sfWidth / 2;

% Determine which superframes are present while the stimulus is on

%% Create a struct for all the relevant timing parameters and save
TD = createStruct(RFTT, sfStarts); % Processing Parameters ======> adjust as needed

if ~exist('TDsavepath', 'var')
    TDsavepath = uigetdir([PDpath, '..\'], 'Select the path to save the timing data in');
    TDsavepath = [TDsavepath, '\'];
end
save([TDsavepath, 'TD.mat'], 'TD')

%% Calculate vUS (g1 fits) for each pre-computed g1 file
% See main_vUS_3D.m

%% Go through pre-computed data files (for all superframes in the experiment) and store
fDim = 4;
% For now, go through only PDI and CDI
% trialsToUse = [1:10];
% for trialInd = trialsToUse
% 
% end
PDIallSF = cell(3, 1);
% CDIallSF = cell(3, 1);
for fi = 1:RFcount

    disp(fi)
    load([PDpath, 'fUSdata-', num2str(fi)], 'PDI', 'CDI')
    for j = 1:3
        PDIallSF{j} = cat(fDim, PDIallSF{j}, PDI{j});
        % CDIallSF{j} = cat(fDim, CDIallSF{j}, CDI{j});
    end
end
% Average PDI and CDI across superframes
PDIA = cell(3, 1);
% CDIA = cell(3, 1);
for j = 1:3
    PDIA{j} = mean(PDIallSF{j}, fDim);
    % CDIA{j} = mean(CDIallSF{j}, fDim);
end
% save([PDpath, 'PDIA_CDIA.mat'], "PDIA", "CDIA")
voxelTimeseriesGUI(PDIallSF{3}, PDIA{3}.^0.5, 'Colormap', 'turbo')

% Load v from vUS fits
VallSF_SG = cell(3, 1);
VallSF_RCA = cell(3, 1);
CallSF_C = cell(3, 1);
for fi = 1:RFcount

    disp(fi)
    load([FRpath, 'fit_results-', num2str(fi)], 'fit_SG', 'fit_RCA', 'fit_C')
    % for j = 1:3
    for j = 3
        VallSF_SG{j} = cat(fDim, VallSF_SG{j}, fit_SG.v);
        VallSF_RCA{j} = cat(fDim, VallSF_RCA{j}, fit_RCA.v);
        CallSF_C{j} = cat(fDim, CallSF_C{j}, fit_C.C);
    end
end

voxelTimeseriesGUI(VallSF_RCA{3}, mean(VallSF_RCA{3}, fDim), 'Colormap', 'turbo')

% Load proxies from tl-fUS fits
CBViallSF = cell(3, 1);
% CBFSiallSF = cell(3, 1);
for fi = 1:RFcount

    disp(fi)
    load([PDpath, 'tlfUSdata-', num2str(fi)], 'CBVi', 'CBFSi')
    for j = 1:3
        CBViallSF{j} = cat(fDim, CBViallSF{j}, CBVi{j});
        % CBFSiallSF{j} = cat(fDim, CBFSiallSF{j}, CBFSi{j});
    end
end
voxelTimeseriesGUI({(PDIallSF{3} - mean(PDIallSF{3}, fDim)) ./ max(PDIallSF{3} - mean(PDIallSF{3}, fDim), [], fDim), CBViallSF{3}}, PDIA{3}.^0.5, 'Colormap', 'turbo')

% Average CBVi across superframes
CBViA = cell(3, 1);
% CDIA = cell(3, 1);
for j = 1:3
    CBViA{j} = mean(CBViallSF{j}, fDim);
end

%% Smoothing test
PDIallSF_orig = PDIallSF;
% gk = [2, 2, 4]; % 3D Gaussian smoothing kernel
gk = [1, 1, 2];
for sfi = 1:size(PDIallSF{3}, fDim)
    PDIallSF{3}(:, :, :, sfi) = imgaussfilt3(PDIallSF_orig{3}(:, :, :, sfi), gk);

end

%% Resample the data onto a uniform time grid, so we can FFT
% rr = 10; % Resampling rate [Hz]
% 
% % PDIallSF(isnan(PDIallSF)) = eps;
% PDIrs = permute(interp1(sfCenters(1:size(PDIallSF{3}, 4)), permute(PDIallSF{3}, [4, 1, 2, 3]), stimTimestamps(1):1/rr:stimTimestamps(end), "linear", "extrap"), [2, 3, 4, 1]);
% % CDIrs = permute(interp1(sfCenters, permute(abs(CDIallSF{3}), [4, 1, 2, 3]), stimTimestamps(1):1/rr:stimTimestamps(end), "linear", "extrap"), [2, 3, 4, 1]);

%% Make a spectrogram of PDI
vs = [length(voxelRange{1}), length(voxelRange{2}), length(voxelRange{3})]; % Volume size [voxels]
% FOR NOW: assume the superframe rate is constant
sfRate = 1/median(diff(sfStarts)); % Superframe rate [Hz]

ws = 20; % Window size (number of samples to include) for the sliding FFT
nw = RFcount - ws + 1; % Number of windows to compute the FFT over (max possible)

fftLen = ws; % FFT length (adjust if needed)

faxis = (-floor(fftLen/2):ceil(fftLen/2)-1) * sfRate / fftLen; % fftshifted FFT bin frequencies [Hz] (zero at index floor(fftLen/2)+1)

s = zeros([vs, nw, fftLen], 'single');

for wi = 1:nw % Go through windows and calculate the spectrogram
    s(:, :, :, wi, :) = fftshift(fft(PDIallSF{3}(:, :, :, wi:wi+ws-1), fftLen, fDim), fDim); % shift only the frequency dim, not the spatial dims
end

figure; imagesc(sfStarts(1:nw), faxis, abs(squeeze(s(40, 40, 40, :, :)))'); axis xy; xlabel('Time [s]'); ylabel('Frequency [Hz]')

%% Inspect the spectrogram per voxel
voxelSpectrogramGUI(s, PDIA{3}.^0.5, 'Time', sfStarts(1:nw), 'Freq', faxis, 'Colormap', 'turbo')