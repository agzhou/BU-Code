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
nta = P.na*2; % R-C and C-R total # of plane waves
numSamples = P.Resource.RcvBuffer.rowsPerFrame/nta; % # samples per element per acquisition (plane wave)
RFData_scrambled = permute(reshape(RcvData', [numChannels, numSamples, nta]), [2, 1, 3]);
% figure; imagesc(squeeze(RFData(:, :, 1)))
RFData_allElem = RFData_scrambled(:, P.Trans.Connector, :); % Unscramble the RF data by using the channel to element map
% figure; imagesc(squeeze(RFData(:, :, 1)))

% Isolate RF for row-column and column-row TX/RX pairs
RFData = zeros(numSamples, P.numElements, nta);
for ai = 1:numel(P.Receive)
    mask = P.Receive(ai).Apod > 0; % Get a mask of which elements have an apodization > 0 (should be only the rows or columns for each plane wave)
    RFData(:, :, ai) = RFData_allElem(:, mask, ai);
end
% genSliderV2(RFData)

%%