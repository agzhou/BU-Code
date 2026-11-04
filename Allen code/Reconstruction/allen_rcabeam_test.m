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
sampleMode = P.Receive.sampleMode;

%% I/Q Demodulation


%% Pass the demodulated but nonbeamformed data to the GPU
% (gpuArray(..))
% (pycupyarray(..))


%% Run the rcabeam beamforming function


%% Convert the output back to a Matlab array (gpu or cpu)

