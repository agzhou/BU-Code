%% Description: main script to run for processing RCA simulation data for KK beamforming

%% Add paths

%% Load data
datapath = 'U:\Projects\Ultrasound\Datasets\Allen Data\RCA Verasonics simulations';
% ...

%% Reshape Verasonics RF data (CODE DOES NOT ACCOUNT FOR MULTIPLE FRAMES YET)
% Originally: [# samples * # angles, # elements]
% Reshape to: [# samples, # elements, # angles]
% numFrames = P.numFramesPerBuffer;
numChannels = P.Resource.Parameters.numRcvChannels;
ntaTX = P.na*2; % R-C and C-R total # of transmit plane waves
anglesTXList = P.angles; anglesTXList = anglesTXList(:); % Make into a column vector
% anglesTX = [[anglesTXList; zeros(size(anglesTXList))], [zeros(size(anglesTXList)); anglesTXList]]; % List of transmit angles. Dimensions [ntaTX, 2 (x and y angle)]
anglesTX = listToAnglesRCA(anglesTXList, 'TX');
maTX = max(anglesTXList); % max TX angle
daTX = mean(diff(anglesTXList));
% if all(diff(anglesTXList) == daTX)
%     error('TX angle increment is not uniform')
% end
numSamples = P.Resource.RcvBuffer.rowsPerFrame/ntaTX; % # samples per element per acquisition (plane wave)
RFData_scrambled = permute(reshape(RcvData', [numChannels, numSamples, ntaTX]), [2, 1, 3]);
% figure; imagesc(squeeze(RFData(:, :, 1)))
RFData_allElem = RFData_scrambled(:, P.Trans.Connector, :); % Unscramble the RF data by using the channel to element map
% figure; imagesc(squeeze(RFData(:, :, 1)))

% Isolate RF for row-column and column-row TX/RX pairs
RFData = zeros(numSamples, P.numElements, ntaTX);
for ai = 1:numel(P.Receive)
    mask = P.Receive(ai).Apod > 0; % Get a mask of which elements have an apodization > 0 (should be only the rows or columns for each plane wave)
    RFData(:, :, ai) = RFData_allElem(:, mask, ai);
end
% genSliderV2(RFData)

%% Get additional parameters
c0 = P.Resource.Parameters.speedOfSound; % Speed of sound in medium [m/s]
RF_fs = P.Receive(1).ADCRate; % Sampling frequency of RF data
ElemPos = P.Trans.ElementPos(:, 1:2)' .* P.wl; % Element positions in [m]. Matrix of size [2 (x, y), # rows + # columns]
ratio = RF_fs/c0;

%% Define the receive angle configurations
naRX = 11; % # of receive angles in one direction

if mod(naRX, 2) ~= 1
    error('# of receive angles (naRX) must be odd')
end
ntaRX = naRX*2; % Total # of receive angles

% Configuration 1: the uniform square
anglesRXList = linspace(-maTX, maTX, naRX - 1);
anglesRXList = [anglesRXList(1:(naRX - 1)/2), 0, anglesRXList((naRX - 1)/2 + 1:end)];
anglesRXList = anglesRXList(:);
anglesRX = listToAnglesRCA(anglesRXList, 'RX');

delta_angles = plotAngleCombos_RCA_func(anglesTX, anglesRX);
figure; plot(delta_angles(:, 1), delta_angles(:, 2), 'o'); axis image

%% Compress the RF Data
% RawDataKK = zeros(numSamples, ntaTX, ntaRX); % Initialize KK-compressed RF data
RFData_hilbert = hilbert(RFData);
RawDataKK = DataCompressKK_RCA(RFData_hilbert, anglesRX, ratio, ElemPos, RF_fs);



