%% Description: writing my own code to call rcabeam from Matlab.
% Dependencies: Parallel Computing Toolbox, Matlab R2026b (for pycupyarray)

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

%% Load RF file
filenum = 1; % Testing
RFData = loadBinaryRF([RFpath, P.fsFileName, sprintf('%05d', filenum), '.rf'], P);
RFData = unstackFrames(RFData, P);

%% If RF is undersampled, do the necessary stuff
sampleMode = P.Receive(1).sampleMode; % Verasonics ADC sample mode
spw = P.Receive(1).samplesPerWave; % Samples per wavelength
decimSampleRate = P.Receive(1).decimSampleRate * 1e6; % [Hz]
demodFreq = P.Receive(1).demodFrequency * 1e6; % [Hz]

sampleDim = 1; numRFSamples = size(RFData, 1);
RFFD = fft(single(RFData), nextpow2(numRFSamples), sampleDim); % RF in the frequency domain
RFFD = RFFD(1:numRFSamples, :, :); % Crop the number of samples

% Up-mix the frequency domain RF to the original (pre-undersampling) baseband
% figure; plot(squeeze(abs(RFFD(:, 111, 3) )))
sampleVec = (0:numRFSamples-1).';
tSample = sampleVec ./ decimSampleRate;
RFFD_upmix = RFFD .* exp(1i.*2.*pi.* demodFreq .* tSample);

RF_upmix = ifft(RFFD_upmix, numRFSamples, sampleDim);

clear RFFD RFFD_upmix

%% I/Q Demodulation


%% Pass the demodulated but nonbeamformed data to the GPU
RFG = gpuArray(RF_upmix);
RFG_cupy = pycupyarray(RFG);

%% Run the rcabeam beamforming function
BF = py.opw(RFG_cupy);

%% Convert the output back to a Matlab array (gpu or cpu)
IQ = gather(BF);
