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

% Reshape the RcvData so the TX angles are a new dimension
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
% figure; imagesc(squeeze(RFData_scrambled(:, :, 1)))

% Unscramble the RF data by using the channel to element map
RFData_allElem = RFData_scrambled(:, P.Trans.Connector, :);
figure; imagesc(squeeze(RFData_allElem(:, :, 1)))

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
% param.t0 = 0;
param.t0 = (P.startDepth + P.Trans.lensCorrection + P.TW.peak)*P.wl/c0; % Time offset [s]

% [~,I] = max(source_sig); % Find the index where the source input signal is maximum (wrt kgrid.dt)
% param.t0 = (kgrid.t_array(I))/param.fc; % Sequence start time (time offset) [s] --> convert the index wrt kgrid.dt to time wrt the RF sampling frequency
% % param.TXdelay = time_delays;
% param.DecimRate = 1;    % Decimation rateCreate beamforming grid

xCoord = (-numElements/2*param.pitch):0.5*P.wl:(numElements/2*param.pitch);  % [m]   Beamformed points x coordinates
% xCoord = ((-numElements/2):1:(numElements/2))*param.pitch;  % [m]   Beamformed points x coordinates
% xCoord = ((-numElements/2):0.5:(numElements/2))*param.pitch;  % [m]   Beamformed points x coordinates
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

interp_method = 'linear';
% [ReconKK, LUTTX, LUTRX] = BeamformKK_RCA(RawDataKK, anglesRX, anglesTX, BFgrid, param);
% [ReconKK] = BeamformKK_RCA(RawDataKK, anglesRX, anglesTX, BFgrid, param, interp_method, 'compounded');
[ReconKKAllAngles] = BeamformKK_RCA(RawDataKK, anglesRX, anglesTX, BFgrid, param, interp_method, 'allAngles');

%% Testing with access to the individual TX-RX pairs' volumes

inds_CR = 1:naTX;
inds_RC = naTX + 1:2*naTX;
anglesTX_CR = anglesTX(inds_CR, :);
anglesTX_RC = anglesTX(inds_RC, :);
anglesRX_CR = anglesRX(inds_CR, :);
anglesRX_RC = anglesRX(inds_RC, :);

delta_angles_CR = zeros(naTX^2, 2);
delta_angles_RC = zeros(naTX^2, 2);
pair_inds_CR = zeros(naTX^2, 2);
pair_inds_RC = zeros(naTX^2, 2);
for ind = 1:naTX
    angleTX = anglesTX_CR(ind, :);
    delta_angles_CR((ind - 1)*naTX + 1:ind*naTX, :) = anglesRX_CR - angleTX;
    pair_inds_CR((ind - 1)*naTX + 1:ind*naTX, :) = [ones(naTX, 1).*ind, (1:naRX)'];
end
for ind = 1:naTX
    angleTX = anglesTX_RC(ind, :);
    % delta_angles_RC(naTX^2 + (ind - 1)*naTX + 1:naTX^2 + ind*naTX, :) = anglesRX_RC - angleTX;
    delta_angles_RC((ind - 1)*naTX + 1:ind*naTX, :) = anglesRX_RC - angleTX;
    pair_inds_RC((ind - 1)*naTX + 1:ind*naTX, :) = [naTX + ones(naTX, 1).*ind, naTX + (1:naRX)'];
end

% ==== Choose which points on the delta angle plot to remove from the final compounded beamformed image ==== %
% Scheme 1: remove the "extra" ks along the axes
% mask_remove_CR = xor(delta_angles_CR(:, 1) == 0, delta_angles_CR(:, 2) == 0);
% mask_remove_RC = xor(delta_angles_RC(:, 1) == 0, delta_angles_RC(:, 2) == 0);
% mask_remove_CR = delta_angles_CR(:, 2) == 0 & delta_angles_CR(:, 1) ~= 0;
% mask_remove_RC = delta_angles_RC(:, 1) == 0 & delta_angles_RC(:, 2) ~= 0;
% mask_remove_CR = delta_angles_CR(:, 2) == 0;
% mask_remove_RC = delta_angles_RC(:, 1) == 0;

% Scheme 2: Remove quadrants 2-4
% mask_remove_CR = ~(delta_angles_CR(:, 2) >= 0 & delta_angles_CR(:, 1) >= 0);
% mask_remove_RC = ~(delta_angles_RC(:, 1) >= 0 & delta_angles_RC(:, 2) >= 0);
mask_remove_CR = ~(delta_angles_CR(:, 2) > 0 & delta_angles_CR(:, 1) >= 0);
mask_remove_RC = ~(delta_angles_RC(:, 1) > 0 & delta_angles_RC(:, 2) >= 0);

% Scheme 3: Remove quadrants 1, 2, 4
mask_remove_CR = ~(delta_angles_CR(:, 2) < 0 & delta_angles_CR(:, 1) <= 0);
mask_remove_RC = ~(delta_angles_RC(:, 1) < 0 & delta_angles_RC(:, 2) <= 0);

% ==== Compound with the masked angles ==== %
ReconKKMasked = zeros(size(BFgrid.X));
% CR
for tai = 1:naTX
    for rai = 1:naRX
        ind = (tai - 1)*naTX + rai;
        if ~mask_remove_CR(ind)
            ReconKKMasked = ReconKKMasked + ReconKKAllAngles(:, :, :, tai, rai);
        end
    end
end
% RC
for tai = 1:naTX
    for rai = 1:naRX
% for tai = naTX + 1:2*naTX
%     for rai = naRX + 1:2*naRX
        ind = (tai - 1)*naTX + rai;
        if ~mask_remove_RC(ind)
            ReconKKMasked = ReconKKMasked + ReconKKAllAngles(:, :, :, naTX + tai, naRX + rai);
        end
    end
end


% Plot the masked delta angles
figure; hold on
plot(rad2deg(delta_angles_CR(~mask_remove_CR, 1)), rad2deg(delta_angles_CR(~mask_remove_CR, 2)), 'o', 'MarkerSize', 8, 'LineWidth', 2)
plot(rad2deg(delta_angles_RC(~mask_remove_RC, 1)), rad2deg(delta_angles_RC(~mask_remove_RC, 2)), 'x', 'MarkerSize', 8, 'LineWidth', 2)
xlim([-max(abs(delta_angles), [], 'all'), max(delta_angles, [], 'all')])
ylim([-max(abs(delta_angles), [], 'all'), max(delta_angles, [], 'all')])
% axis image
axis square
title('Delta angles'); xlabel('x angle [deg]'); ylabel('y angle [deg]'); fontsize(20, 'points')
hold off

% Plot MIP
ylims = [0, 5];
KKMaskedMIP_fh = figure;
% imagesc(xCoord*1e3, zCoord*1e3, squeeze(max(abs(ReconKKMasked), [], 1))'); colormap gray
imagesc(xCoord*1e3, zCoord*1e3, squeeze(max(abs(ReconKKMasked), [], 2))'); colormap gray
title('KK')
xlabel('x [mm]')
ylabel('z [mm]')
% KKMIP_fh.Position(4) = KKMIP_fh.Position(3)*( max(zCoord) - min(zCoord) )/( max(xCoord) - min(xCoord) );
fontsize(20, 'points')
ylim(ylims)

% % Plot MIP
% ylims = [0, 5];
% KKMaskedMIP_fh = figure;
% % imagesc(xCoord*1e3, zCoord*1e3, squeeze(max(abs(ReconKKMasked), [], 1))'); colormap gray
% imagesc(xCoord*1e3, zCoord*1e3, squeeze(max(abs(conj(ReconKKMasked)), [], 2))'); colormap gray
% title('KK')
% xlabel('x [mm]')
% ylabel('z [mm]')
% % KKMIP_fh.Position(4) = KKMIP_fh.Position(3)*( max(zCoord) - min(zCoord) )/( max(xCoord) - min(xCoord) );
% fontsize(20, 'points')
% ylim(ylims)

volumeViewer(abs(ReconKKMasked))
%% Plot KK MIP
ylims = [0, 5];
KKMIP_fh = figure;
% imagesc(xCoord*1e3, zCoord*1e3, squeeze(max(abs(ReconKK), [], 1))'); colormap gray
imagesc(xCoord*1e3, zCoord*1e3, squeeze(max(abs(ReconKK), [], 2))' .^ 1); colormap gray
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
imagesc(vs_xCoord*1e3, vs_zCoord*1e3, squeeze(max(abs(IQ(:, :, :)), [], 1))' .^ 1); colormap gray
title('DAS')
xlabel('x [mm]')
ylabel('z [mm]')
% DASMIP_fh.Position(4) = DASMIP_fh.Position(3)*( max(vs_zCoord) - min(vs_zCoord) )/( max(vs_xCoord) - min(vs_xCoord) );
fontsize(20, 'points')
ylim(ylims)