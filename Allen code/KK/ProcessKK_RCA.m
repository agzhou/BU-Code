%% Description: main script to run for processing RCA simulation data for KK beamforming

%% Add paths

%% Load data
% datapath = 'U:\Projects\Ultrasound\Datasets\Allen Data\RCA Verasonics simulations';
% ...

%% Reshape Verasonics RF data (CODE DOES NOT ACCOUNT FOR MULTIPLE FRAMES YET)
% Originally: [# samples * # angles, # elements]
% Reshape to: [# samples, # elements, # angles]
% numFrames = P.numFramesPerBuffer;
numChannels = P.Resource.Parameters.numRcvChannels;
ntaTX = P.na*2; % R-C and C-R total # of transmit plane waves
naTX = ntaTX/2;
anglesTXList = P.angles; anglesTXList = anglesTXList(:); % Make into a column vector
% anglesTX = [[anglesTXList; zeros(size(anglesTXList))], [zeros(size(anglesTXList)); anglesTXList]]; % List of transmit angles. Dimensions [ntaTX, 2 (x and y angle)]
anglesTX = listToAnglesRCA(anglesTXList, 'TX');
maTX = max(anglesTXList); % max TX angle
daTX = mean(diff(anglesTXList));
% if all(diff(anglesTXList) == daTX)
%     error('TX angle increment is not uniform')
% end
% numSamples = P.Resource.RcvBuffer.rowsPerFrame/ntaTX; % # samples per element per acquisition (plane wave)
numSamples = P.Receive(1).endSample; % # samples per element per acquisition (plane wave)
% RFData_scrambled = permute(reshape(RcvData', [numChannels, numSamples, ntaTX]), [2, 1, 3]);

RFData_scrambled = zeros(numSamples, numChannels, ntaTX);
% first half of rcvs
for ai = 1:naTX
    RFData_scrambled(:, :, ai)  = RcvData(P.Receive(ai).startSample:P.Receive(ai).endSample, :);
end

% second half of rcvs
for ai = 1:naTX
    RFData_scrambled(:, :, naTX + ai) = RcvData(P.Receive(naTX + ai).startSample:P.Receive(naTX + ai).endSample, :);
end
% RFData_scrambled = permute(reshape(RcvData', [numChannels, numSamples, ntaTX]), [2, 1, 3]);

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
numElements = P.numElements;
RF_fs = P.Receive(1).ADCRate * 1e6; % Sampling frequency of RF data [Hz]
ElemPos = P.Trans.ElementPos(:, 1:2)' .* P.wl; % Element positions in [m]. Matrix of size [2 (x, y), # rows + # columns]
ratio = RF_fs/c0;

%% Define the receive angle configurations
% naRX = 11; % # of receive angles in one direction
naRX = naTX

if mod(naRX, 2) ~= 1
    error('# of receive angles (naRX) must be odd')
end
ntaRX = naRX*2; % Total # of receive angles

% Configuration 1: the uniform square
% anglesRXList = linspace(-maTX, maTX, naRX - 1);
% anglesRXList = [anglesRXList(1:(naRX - 1)/2), 0, anglesRXList((naRX - 1)/2 + 1:end)];
daRX = daTX;
offset = 0.5
% offset = naRX/naTX/2
anglesRXList = [-((naRX-1)/2 - offset)*daRX:daRX: -offset*daRX, 0, offset*daRX:daRX:((naRX-1)/2 - offset)*daRX];

anglesRXList = anglesRXList(:);
anglesRX = listToAnglesRCA(anglesRXList, 'RX');

% % Testing
% anglesRX = fliplr(anglesTX);

delta_angles = plotAngleCombos_RCA_func(anglesTX, anglesRX);
% figure; plot(delta_angles(:, 1), delta_angles(:, 2), 'o'); axis image; title('Delta angles'); xlabel('y angle [deg]'); xlabel('y angle [deg]'); fontsize(20, 'points')
plotDeltaAngles_RCA(anglesTX, anglesRX)

%% Get TX time delays
spw = P.Receive(1).samplesPerWave; % # samples per wavelength
time_delays_TX = zeros(P.na * 2, numElements * 2); % In units of [samples]. Dimensions are [total # TX angles, total # elements]
for tai = 1:length(P.TX)
    time_delays_TX(tai, :) = P.TX(tai).Delay .* spw;
    % time_delays_TX(tai, :) = P.TX(tai).Delay .* 1;
end

%% Compress the RF Data
% RawDataKK = zeros(numSamples, ntaTX, ntaRX); % Initialize KK-compressed RF data
RFData_hilbert = hilbert(RFData);
% RFData_hilbert = RFData;

RawDataKK = DataCompressKK_RCA(RFData_hilbert, anglesRX, ratio, ElemPos, time_delays_TX);
% RawDataKK = DataCompressKK_RCA_nocompensation(RFData_hilbert, anglesRX, ratio, ElemPos);
% figure; imagesc(squeeze(abs(RawDataKK(:, 6, :))))

%% Beamforming and other key parameters' definitions
numElements = P.numElements; % # of elements in one dimension (# rows = # columns)
param.fs = RF_fs;                           % [Hz]   sampling frequency
param.pitch = P.Trans.spacingMm/1e3; % Element pitch [m]
param.fc = P.Trans.frequency*1e6;                       % [Hz]   center frequency
param.c = c0;                               % [m/s]  longitudinal sound speed
% param.fnumber = [0.6, 0.6];                        % [ul]   receive f-number
param.elements = ElemPos; % Element coordinates (x, y) [m]
wavelength = param.c/param.fc;              % [m] convert from wavelength to meters
% samplesPerWave = param.fs/param.fc;     % the number of samples per wavelength
% note: this is off by a factor of two because you also account for
% roundtrip time. In otherwords, there are 4 samples per wavelength, but in
% practice that becomes 8 since you are also accounting for time to go to
% and from the transducer.
param.t0 = 0; % not sure..................................................

% [~,I] = max(source_sig); % Find the index where the source input signal is maximum (wrt kgrid.dt)
% param.t0 = (kgrid.t_array(I))/param.fc; % Sequence start time (time offset) [s] --> convert the index wrt kgrid.dt to time wrt the RF sampling frequency
% % param.TXdelay = time_delays;
% param.DecimRate = 1;    % Decimation rateCreate beamforming grid

xCoord = ((-numElements/2):0.5:(numElements/2))*param.pitch;  % [m]   Beamformed points x coordinates
yCoord = xCoord;
zbounds_mm = [0, 5]; % Z bounds/extents [mm]
zbounds = zbounds_mm ./ 1e3; % Z bounds/extents [m]
zCoord = zbounds(1):0.5*P.wl:zbounds(2);   % [m]    Beamformed points z coordinates
% zCoord = (1:0.025:32)*wavelength;   % [m]    Beamformed points z coordinates
[X, Y, Z] = meshgrid(xCoord, yCoord, zCoord);

BFgrid = struct('X', X, 'Y', Y, 'Z', Z); % Struct for the beamforming grid

%% KK Beamforming

% % ReconKK = BeamformKK_RCA(RawDataKK, anglesRX, anglesTX, BFgrid, param);
% [ReconKK, LUTTX, LUTRX] = BeamformKK_RCA(RawDataKK, anglesRX, anglesTX, BFgrid, param);
% 
% % figure; imagesc(squeeze(max(abs(ReconKK), [], 1))')
% 
% % Testing: look at individual volumes for each plane wave
% figure; imagesc(squeeze(max(abs(ReconKK(:, :, :, 6, 6)), [], 1))')
% test = sum(ReconKK, 5);
% figure; imagesc(squeeze(max(abs(test(:, :, :, 6)), [], 1))')
% fulltest = squeeze(sum(ReconKK, [4, 5]));
% figure; imagesc(squeeze(max(abs(fulltest), [], 1))')
% 
% genSliderV2(squeeze(max(abs(ReconKK), [], [1, 5])))
% % temp = squeeze(max(abs(fulltest), [], 1))';
% % figure; imagesc([temp(:, 41:end), temp(:, 1:40)])
% 
% % Plot images from the CR and RC sums individually
% figure; imagesc(abs(squeeze(sum(ReconKK(:, :, :, 12:22, 12:22), [1, 4, 5]))))
% figure; imagesc(abs(squeeze(sum(ReconKK(:, :, :, 1:11, 1:11), [1, 4, 5])))')
% 
% [ReconKK] = BeamformKK_RCA_notsummingcorrectly(RawDataKK, anglesRX, anglesTX, BFgrid, param);
% ReconKKCR = sum(ReconKK(:, :, :, 1:naTX, 1:naRX), [4, 5]);
% ReconKKRC = sum(ReconKK(:, :, :, naTX + 1:2*naTX, naRX + 1:2*naRX), [4, 5]);
% properCPWC = ReconKKCR + ReconKKRC;
% figure; imagesc(squeeze(max(abs(properCPWC), [], 1))')


% [ReconKK, LUTTX, LUTRX] = BeamformKK_RCA(RawDataKK, anglesRX, anglesTX, BFgrid, param);
[ReconKK] = BeamformKK_RCA(RawDataKK, anglesRX, anglesTX, BFgrid, param);

%% Plot KK MIP
ylims = [0, 5];
KKMIP_fh = figure;
imagesc(xCoord*1e3, zCoord*1e3, squeeze(max(abs(ReconKK), [], 1))'); colormap gray
title('KK')
xlabel('x [mm]')
ylabel('z [mm]')
% KKMIP_fh.Position(4) = KKMIP_fh.Position(3)*( max(zCoord) - min(zCoord) )/( max(xCoord) - min(xCoord) );
fontsize(20, 'points')
ylim(ylims)

%% Plot Verasonics DAS
vs_numGridPts = PData.Size;
vs_xCoord = ((1:vs_numGridPts(1)).*PData.PDelta(1) - 1 + PData.Origin(1) ).*P.wl;
vs_zCoord = ((1:vs_numGridPts(3)).*PData.PDelta(3) - 0 + PData.Origin(3) ).*P.wl;

DASMIP_fh = figure;
% imagesc(vs_xCoord, vs_zCoord, squeeze(max(abs(IQ(:, :, 1:length(zCoord))), [], 1))')
imagesc(vs_xCoord*1e3, vs_zCoord*1e3, squeeze(max(abs(IQ(:, :, :)), [], 1))' .^ 0.5); colormap gray
title('DAS')
xlabel('x [mm]')
ylabel('z [mm]')
% DASMIP_fh.Position(4) = DASMIP_fh.Position(3)*( max(vs_zCoord) - min(vs_zCoord) )/( max(vs_xCoord) - min(vs_xCoord) );
fontsize(20, 'points')
ylim(ylims)