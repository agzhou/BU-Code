%% Description: writing my own code to call rcabeam from Matlab.
% Dependencies: Parallel Computing Toolbox, Matlab R2026b (for pycupyarray)

%% Add the rcabeam folder (and subfolders) to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
rcabeamCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "GitHub"))), '\') + "\rcabeam\");
addpath(genpath(rcabeamCodePath))

%% Load RF and params
RFpath = uigetdir('G:\Allen\Data\', 'Select the RF data path');
RFpath = [RFpath, '\'];
if ~exist('P', 'var')
    load([RFpath, 'params.mat']) % load acquisition parameters
end

IQpath = uigetdir(RFpath, 'Select the folder to save reconstructed data to');
IQpath = [IQpath, '\'];

% Prompt for parameter user input
parameterPrompt = {'Start file number', 'End file number'};
parameterDefaults = {'1', ''};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

startFile = str2double(parameterUserInput{1}); % File to start reconstructing from
endFile = str2double(parameterUserInput{2});   % File to stop reconstructing on

clear parameterPrompt parameterDefaults parameterUserInput

%% Get parameters
nf = P.numFramesPerBuffer; % # of frames per superframe
na = P.na; % # of angles (for RC or CR)
nspa = P.Resource.RcvBuffer(1).rowsPerFrame / nf / (2*na); % of samples per acquisition
ntc = P.Resource.Parameters.numRcvChannels; % # of total channels
nc = length(P.Trans.Connector); % # of probe channels (row + col)
wl = P.wl;
Trans = P.Trans;
numElements = P.Trans.numelements/2; % # of rows or cols
angles = P.angles;
c = P.Resource.Parameters.speedOfSound; % Speed of sound [m/s]

% Sampling
sampleMode = P.Receive(1).sampleMode; % Verasonics ADC sample mode
spw = P.Receive(1).samplesPerWave; % Samples per wavelength
decimSampleRate = P.Receive(1).decimSampleRate * 1e6; % [Hz]
demodFreq = P.Receive(1).demodFrequency * 1e6; % [Hz]

%% Load RF file
filenum = 1; % Testing
RFData = loadBinaryRF([RFpath, P.fsFileName, sprintf('%05d', filenum), '.rf'], P);
RFData = unstackFrames(RFData, P);
RFData = reshape(RFData, [nspa, 2*na, ntc, nf]); % Reshape to dimensions (samples, angles, channels, frames)
RFData = RFData(:, :, P.Trans.Connector, :); % Map the receive channels properly

%% If RF is undersampled, do the necessary stuff

sampleDim = 1;
numRFSamples = size(RFData, 1);
sampleVec = (0:numRFSamples-1).';
tSample = sampleVec ./ decimSampleRate; % [Micro]-Time corresponding to each sample

% Demodulate (remove carrier) in the time domain
RF_demod = 2.* single(RFData) .* exp(1i.*2.*pi.* -demodFreq .* tSample); % Multiply by 2 because the original signal is real{envelope * exp(2j*pi*demodFreq*t)} = 1/2*envelope * exp(2j*pi*demodFreq*t) + 1/2*envelope * exp(-2j*pi*demodFreq*t) --> property of complex numbers. Then, demodulating gives you 1/2*envelope + 1/2*envelope*exp(-4j*pi*demodFreq*t), and you keep the first term

% FFT and lowpass the demodulated RF
RFFD = fft(RF_demod, numRFSamples, sampleDim); % RF in the frequency domain
BW = abs(diff(P.Trans.Bandwidth)) * 1e6; % Bandwidth of the probe [Hz]
faxis = linspace(-decimSampleRate/2, decimSampleRate/2, size(RFFD, sampleDim));
LPF_mask = abs(faxis) <= min(BW/2, decimSampleRate/4);
% figure; yyaxis left; plot(faxis, abs(RFFD(:, 1, 120, 1))); yyaxis right; plot(faxis, LPF_mask)
RFFD_LPF = fftshift(RFFD, sampleDim); RFFD_LPF(~LPF_mask, :, :, :) = 0; % Apply the low pass filter
% figure; yyaxis left; plot(faxis, abs(RFFD_LPF(:, 1, 120, 1))); yyaxis right; plot(faxis, LPF_mask)

RF_LPF = ifft(ifftshift(RFFD_LPF, sampleDim), [], sampleDim);

% Get the start sample time of each acquisition and shift the RF for each
% acquisition accordingly (start sample time wrt each acquisition
% corresponds to the first element that fired). Phase depends on the
% carrier frequency, demodFreq.
elementx = (wl .* Trans.ElementPos(1:numElements, 1))';
elementy = (wl .* Trans.ElementPos(numElements + 1 : numElements*2, 2))';
x_end = max(elementx); % End element position [m]
y_end = max(elementy);
d_correction_1way = Trans.lensCorrection * wl; % TX and RX each pass through the lens
pathConst = wl * P.TW.peak + 2 * d_correction_1way - 2 * P.startDepth * wl;
t0 = -(pathConst + x_end * abs(sin(angles))) / c;

RF = RF_LPF .* repmat(exp(-1i.*2.*pi .* demodFreq .* t0), 1, 2);


%% Pass the demodulated but nonbeamformed data to the GPU
RFG = gpuArray(RF);
RFG_cupy = pycupyarray(RFG);

%% Run the rcabeam beamforming function
BF = py.opw(RFG_cupy);

%% Convert the output back to a Matlab array (gpu or cpu)
IQ = gather(BF);
