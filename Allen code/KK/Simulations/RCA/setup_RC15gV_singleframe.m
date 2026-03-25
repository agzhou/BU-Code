
%% Description
% Single acquisition and saving of RF data with the RC15gV probe
% Saves all defined # of frames in one file
% C-R and R-C pairs of TX-RX
% Uses saveRcvData external function for saving
% Note: update the savepath variable as needed

% Collects nbuf buffers of nf frames with some duty cycle for saving delays

%% Specify system parameters
clearvars

codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
addpath(AllenVerasonicsCodePath)

cd 'C:\Users\agzhou\Vantage-4.9.7-2505271400'
% cd 'C:\Users\BOAS-US\Desktop\Vantage-4.9.7-2505271400'
% cd 'C:\Users\BOAS-US\Desktop\Vantage-4.9.5-2409181500'
activate

%% Testing for TX angles
% % Use different TX angles for R-C and C-R
% maxAngle = 5; % [deg]
% angleRange = [-maxAngle, maxAngle].*pi/180; % Angle range in radians
% 
% % For now:
% %   - # RC transmits is odd, includes angle = 0
% %   - # CR transmits is even, does not include angle = 0
% naTX_RC = 5; % # of transmit angles in one dimension for row-column transmit
% naTX_CR = naTX_RC - 1;
% anglesTX_RC = linspace(angleRange(1), angleRange(2), naTX_RC); % R-C TX angles (could be X or Y angles)
% daTX = mean(diff(anglesTX_RC));
% anglesTX_CR = -(naTX_CR - 1) / 2 * daTX:daTX:(naTX_CR - 1) / 2 * daTX;
% 
% figure;
% plot(zeros(length(anglesTX_RC), 1), anglesTX_RC, '.')
% hold on
% plot(zeros(length(anglesTX_CR), 1), anglesTX_CR, 'x')
% hold off

%%
savepath = uigetdir('F:\', 'Select the save path');
savepath = [savepath, '\'];

parameterPrompt = {'Probe voltage [V]', 'Start depth [mm]', 'End depth [mm]', 'Pulse Repetition Frequency [Hz]', 'Frame rate [Hz]', 'Number of angles', 'Maximum angle [degrees]', 'Probe frequency [MHz]', 'Speed of sound [m/s]', 'Simulate Mode (0-off, 1-on, 2-RcvLoop)', 'Save RcvData (0-no, 1-yes)', 'Number of frames per superframe'}; % 'Save RF data (0-no, 1-yes)', 
parameterDefaults = {'20', '0', '8', '60000', '450', '11', '5', '13.6', '1540', '1', '1', '1'};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% Store the user inputs for parameters into the corresponding variables
initialVoltage = str2double(parameterUserInput{1});
startDepthMM = str2double(parameterUserInput{2});
endDepthMM = str2double(parameterUserInput{3});
PRF = str2double(parameterUserInput{4});
frameRate = str2double(parameterUserInput{5}); % subframe (volume) rate (Hz)
na = str2double(parameterUserInput{6});
maxAngle = str2double(parameterUserInput{7});
probe_freq = str2double(parameterUserInput{8});
speedOfSound = str2double(parameterUserInput{9});
simMode = str2double(parameterUserInput{10});
saveRcvDataFlag = str2double(parameterUserInput{11});
numFramesPerSF = str2double(parameterUserInput{12});

bufferIndex = 0;
runVSX = 1;
movePointsOrNot = 0;
numChannels = 256; % enable channels

numBuffers = 1;
bufferDutyCycle = 1/3;
% disp(num2str(numFramesPerBuffer / frameRate / bufferDutyCycle))

angleRange = [-maxAngle, maxAngle].*pi/180; % Angle range in radians

% Need at least 2 acquisitions to use multiple angles. 
% Otherwise, set angle to 0 degrees.
if na >= 2 
    angles = linspace(angleRange(1), angleRange(2), na);
else
    angles = 0;
end

% numAngles = length(angles);
pair = 2; % The R-C and C-R pair of acquisitions per angle

% Resource is a structure, define system parameters
Resource.Parameters.numTransmit = numChannels; % number of transmit channels
Resource.Parameters.numRcvChannels = numChannels; % number of receive channels
% Resource.Parameters.connector = 1; % transducer connector to use since the current plate for the 256 bit system is split into two 128 bit connectors. 1 is left and 2 is right
Resource.Parameters.speedOfSound = speedOfSound; % speed of sound in m/s, the 1540 is for average human tissue

%% Define Transducer

Trans.name = 'RC15gV'; 
Trans.frequency = probe_freq; % Not needed if using the default center frequency
Trans.units = 'wavelengths'; % or mm

Trans = computeTrans(Trans); % Generate required attributes for the probe into the Trans structure; e.g., the transducer element positions
% Trans.maxHighVoltage = ; % set maximum high voltage that is allowed to the transducer

L = Trans.spacingMm*Trans.numelements/2/1e3; % m
wl = Resource.Parameters.speedOfSound / Trans.frequency / 1e6; % m

startDepth = startDepthMM/1e3/wl; % start depth in wavelengths
endDepth = endDepthMM/1e3/wl;

%% angles
% angpitch = wl / (Trans.spacingMm*Trans.numelements / 2 / 1e3);
% angles = -(na - 1) / 2 * angpitch : angpitch : (na - 1) / 2 * angpitch
%% enable time tag
TimeTagEna = 0;
% 0: disable
% 1: enable but don't reset counter
% 2: enable and reset counter

%% Simulation things - Media structure (define scattering points and attenuation)
Resource.Parameters.simulateMode = simMode; % run script in simulate mode. Set to 0 if not

zd_mm = 4; % in mm
zd = zd_mm/wl/1e3; % z position of point scatterer, in [wavelengths]

% Set up Media model for the simulation, which generates the scattering points with 3D location
% and reflectivity. For 1D transducer arrays, they are aligned on the
% x-axis with the center at x = 0, and scan depth is in z.
Media.MP(1, :) = [0, 0, zd, 1.0]; % [x, y, z, reflectivity]. x, y, z are defined as # of wavelengths.
% Media.MP(2, :) = [30, 30, 70, 1.0]; % [x, y, z, reflectivity]. x, y, z are defined as # of wavelengths.
% Media.MP(3, :) = [20, -20, 100, 1.0]; % [x, y, z, reflectivity]. x, y, z are defined as # of wavelengths.
% Media.MP(1, :) = [30, 30, 70, 1.0]; % [x, y, z, reflectivity]. x, y, z are defined as # of wavelengths.

% new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vesselX = 30e-6;    % x dimension
% vesselY = 30e-6;    % x dimension
% vesselZ = endDepthMM * 1e-3; % z dimension
% flow_v_mm_s = 300;
%  
% Media.MP = randomPts3D_func(vesselX, vesselY, vesselZ, wl);

Media.attenuation = 0;
% Media.attenuation = -0.7; % media attenuation in dB/cm/MHz

Media.function = 'movePointsZ3D'; % move points in _ dimension after each frame

%% PData structure (Pixel Data --> image reconstruction range)
% For 2D scans and slices of 3D scans, it's always a rectangular area at a
% fixed location in the transducer coord system

numElements = Trans.numelements./2; % the structure gives # row elements + # column elements

PData.PDelta = [Trans.spacing, Trans.spacing, 0.5]; % Spacing between pixels in x, y, z, in wavelengths

PData.Coord = 'rectangular'; % rectangular coords, could change to polar or spherical
% Set PData array dimensions --> # of rows, columns, sections (planes
% parallel to the xy plane)
% For a 3D scan, rows - y axis, columns - x axis, sections - z axis
PData.Size(1) = ceil(numElements.*Trans.spacing./PData.PDelta(2)); % # rows
PData.Size(2) = ceil(numElements.*Trans.spacing./PData.PDelta(1)); % # cols
PData.Size(3) = ceil((endDepth - startDepth)./PData.PDelta(3)); % sections

% Define the location (x, y, z) of the upper left corner of the array
half_probe_dist = (numElements-1)./2.*Trans.spacing;
PData.Origin = [-half_probe_dist, half_probe_dist, startDepth];
% PData.Origin = [-half_probe_dist, -half_probe_dist, startDepth];

% Upper left corner if you look aligned with positive z

% Set a local region to view/use for processing
PData.Region(1) = struct('Shape',struct('Name','PData'));

PData.Region(2).Shape = struct('Name', 'Slice', 'Orientation', 'xz', ...
                            'oPAIntersect', PData.Origin(2) - (numElements-1).*Trans.spacing./2); % out of Plane Axis Intersection
PData.Region(3).Shape = struct('Name', 'Slice', 'Orientation', 'yz', ...
                            'oPAIntersect', PData.Origin(1) + (numElements-1).*Trans.spacing./2);
PData.Region(4).Shape = struct('Name', 'Slice', 'Orientation', 'xy', ...
                            'oPAIntersect', Media.MP(3)); % currently set to the plane intersecting the only scatter point

%     'Position', [0, 0, 10], ...
%                      'width', PData.Size(2), 'height', PData.Size(1)./2);
                    % Position is relative to the global coords
% PData.Region = computeRegions(PData);

%% Transmission Waveform (TW)
tw.A = Trans.frequency; % frequency of transmission pulse, sets half cycle period of the waveform...
tw.B = 0.67; % amount of time (0.1 - 1.0) that the transmission drivers are active in the half cycle period. Controsl how much power is delivered.
             % Apparently using B = 0.67 approximates a sine wave.
tw.C = 2; % number of half cycles in the transmission waveform. 2 half cycles = 1 full cycle burst
tw.D = 1; % initial polarity of the first half cycle (1 = +, 0 = -)
TW(1).type = 'parametric';
TW(1).Parameters = [tw.A, tw.B, tw.C, tw.D];

% Note: can modify B for transmission apodization......... See tutorial, "The modification of the B parameters for Transmit apodization can be
% performed automatically by VSX by setting the weighting values in the TX.Apod attribute of a TX structure..."

% In the case where one really would like independent parametric waveforms on individual
% channels, the TW.Parameters array can be expanded to a 128 (or 256) row, two
% dimensional array, where each row specifies the waveform for the transmitter of the same
% number as the row index.
% The Vantage system supports other methods for defining transmit waveforms, including the
% TW.types of ‘envelope’, ‘pulseCode’, and ‘function’. These additional methods are
% described in the Sequence Programming Manual.

TPC.hv = initialVoltage;

%% Transmit action - TX structure

% Need a TX structure for each unique transmit action in the imaging
% sequence

% na*2 transmissions of a plane wave in pairs, one by all row elements and then one by all
% column elements
TX = repmat(struct('waveform', 1, ...
                   'focus', 0, ... % plane wave
                   'Steer', [0.0, 0.0], ... % theta, alpha (beam angle projected in xz from +z axis, beam angle wrt xz)
                   'Apod', zeros(1, Trans.numelements)), 1, na*2);
for n = 1:na
    TX(n).Apod(1:Trans.numelements/2) = ones(1, Trans.numelements/2); % Turn on columns (y)
    TX(n).Steer = [angles(n), 0];
    TX(n).Delay = computeTXDelays(TX(n));
end

for n = 1:na
    TX(na + n).Apod(Trans.numelements/2 + 1 : end) = ones(1, Trans.numelements/2); % Turn on rows (x)
    TX(na + n).Steer = [0, angles(n)];
    TX(na + n).Delay = computeTXDelays(TX(na + n));
end


%% Define Time Gain Control waveform (TGC)
% Accounts for decrease in amplitude of echoes for longer distance traveled

% TGC curve definition
% TGC.CntrlPts = [0 785.2216 1023 1023 1023 1023 1023 1023];
TGC.CntrlPts = [1023 1023 1023 1023 1023 1023 1023 1023];
% TGC(1).CntrlPts = [500,590,650,710,770,830,890,950]; % 0 to 1023, minimum to maximum gain
                                                     % Values represent the
                                                     % gain at increasing
                                                     % depth in the
                                                     % acquisition period.
                                                     % They are equally
                                                     % distributed over the
                                                     % 0 to rangeMax depth
                                                     % (in wavelengths)
TGC(1).rangeMax = endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC); % Parameters can be adjusted later with GUI sliders

%% Receiver array object
maxAcqLength = ceil(sqrt(endDepth^2 + 2*(numElements*Trans.spacing)^2)); % account for the longest distance an echo could travel
Receive = repmat(struct('Apod', zeros(1, Trans.numelements), ... 
                        'startDepth', startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0, ...
                        'LowPassCoef', [], ...
                        'InputFilter', []), 1, pair * numBuffers * numFramesPerSF * na);
j = 1;
% an = 0;
for nbuf = 1:numBuffers
    
    for nf = 1:numFramesPerSF
        an = 0; % acquisition number
        
        % Move points after all the acquisitions for one frame
        Receive(j).callMediaFunc = movePointsOrNot;
    %     Receive(j).mode = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for n = 1:na
            an = an + 1;
            Receive(j).bufnum = nbuf;
            Receive(j).framenum = nf;
            Receive(j).acqNum = an;
            Receive(j).Apod(Trans.numelements/2 + 1 : end) = ones(1, Trans.numelements/2);
            j = j + 1;
        end
    
        for n = 1:na
            an = an + 1;
            Receive(j).bufnum = nbuf;
            Receive(j).framenum = nf;
            Receive(j).acqNum = an;
            Receive(j).Apod(1:Trans.numelements/2) = ones(1, Trans.numelements/2);
            j = j + 1;
        end
        
    end
end

%% Allocate storage space for RF data acquisition
% RcvBuffer dimensions are (samples, channels, frames)
% RcvBuffer is accessible in Matlab as a cell array

% If not set by the user, VSX will set the effective ADC sample rate to 4x
% the transducer center frequency (system max seems to be 250 MHz). There
% are preset values for sampling rate, and the system will choose the
% closest one above whatever the user sets.

% RcvBuffer dimensions: (samples, channels, frames, pages)


%%%% from Nikunj's SetUpCustomIntegratedRecon.m code
if strcmp(Receive(1).sampleMode,'custom')
    error('No handling of condition for custom Receive sampling. Refer to VsUpdate line 712 to implement');
else
    fs = 4*Trans.frequency;
    samplesPerWave = 4;
end

% if statement included to match verasonics automatic extension to
% multiples of 128 samples
nSmpls = 2*(maxAcqLength - startDepth) * samplesPerWave; % maxAcqLength is the Receive(1).endDepth
% nSmpls = 2*(Receive(1).endDepth - Receive(1).startDepth) * samplesPerWave;
if abs(round(nSmpls/128) - nSmpls/128) < .01
    numRcvSamples = 128*round(nSmpls/128);
else
    numRcvSamples = 128*ceil(nSmpls/128);
end

% startSample = (0:(na-1))*numRcvSamples + 1;
% endSample = startSample + numRcvSamples - 1;
%%%%

% spw = 3.6765; % samples per wave, it isn't always exactly 4... check p107
% nspa = spw*(2*(Receive(1).endDepth - Receive(1).startDepth));
% nspa = 128 * ceil(nspa/128); % # samples per acquisition
% maxAcqLength_adjusted = nspa / spw / 2;

maxAcqLength_adjusted = numRcvSamples / samplesPerWave / 2;

for nbuf = 1:numBuffers
    Resource.RcvBuffer(nbuf).rowsPerFrame = numRcvSamples * na * 2;
    Resource.RcvBuffer(nbuf).colsPerFrame = Resource.Parameters.numRcvChannels; % Usually 1:1 to # of receive channels available in the system. Can change to 256 with the 2D probe and new connector plate.
    Resource.RcvBuffer(nbuf).numFrames = numFramesPerSF; % minimum # frames of RF data to acquire; RcvBuffer contains all the data needed for a whole frame, including multiple acquisition passes needed for reconstruction. Software can re-process RcvBuffer frames
    Resource.RcvBuffer(nbuf).datatype = 'int16'; % 16 bit signed integers are the only supported datatype
end
% Commenting below section because it doesn't work for the second set of
% TXs
% for lss = 1:length(startSample)
%     Receive(lss).startSample = startSample(lss);
%     Receive(lss).endSample = endSample(lss);    
% %     Receive(lss).decimSampleRate = samplesPerWave * Trans.frequency;
%     Receive(lss).decimSampleRate = 62.5;
% 
% end


Resource.Parameters.verbose = 2; % Describe errors in varying levels
% Resource.InterBuffer(1).pagesPerFrame = pair*na*numSubFrames; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numSamplesInBuffer = Resource.RcvBuffer(1).rowsPerFrame * Resource.RcvBuffer(1).colsPerFrame * Resource.RcvBuffer(1).numFrames
numGBInBuffer = numSamplesInBuffer ./ 1024^3 * 2 % # samples * (2 bytes per int16 sample) 
numSamplesPerBufferFrame = Resource.RcvBuffer(1).rowsPerFrame * Resource.RcvBuffer(1).colsPerFrame
numGBPerBufferFrame = numSamplesPerBufferFrame ./ 1024^3 * 2 % # samples * (2 bytes per int16 sample) 

if numGBPerBufferFrame > 2
    warning('Buffer size per frame is too large, exiting')
    return
end

if ((maxAcqLength_adjusted + (endDepth))*wl / speedOfSound) > 1/PRF
    error('Error: the PRF is too high, it will send the next transmission before the previous transmission reflects from the deepest part of the region')
end

%% Reconstruction
% numRegions = 3;
% 
% Resource.ImageBuffer(1).numFrames = numSupFrames; % Define an ImageBuffer with a # of frames
% Resource.InterBuffer(1).numFrames = numSupFrames; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Resource.InterBuffer(1).numFrames = 1;
% 
% % Recon = struct('senscutoff', 0.6, ... % Threshold for which the reconstruction doesn't consider an element's contribution due to directivity of the element, for a certain pixel (whose echoes are at an angle to the element). Should be in radians.
% %                'pdatanum', 1, ... % Which PData structure to use
% %                'rcvBufFrame', -1, ... % Use the most recently transferred frame
% %                'IntBufDest', [1, 1], ... % idk but it's for the IQ (complex) data
% %                'ImgBufDest', [1, -1], ... % [buffer #, frame #] Auto-increment ImageBuffer for each reconstruction???? % something is [first/oldest frame, last/newest frame]
% %                'RINums', [1:2*na]); % The ReconInfo structure #(s). Each Recon must have its own unique set of ReconInfo #s
% 
% sco = 0.6; %%%%
% % sco = 0.4;
% Recon = struct('senscutoff', sco, ... % Threshold for which the reconstruction doesn't consider an element's contribution due to directivity of the element, for a certain pixel (whose echoes are at an angle to the element). Should be in radians.
%                'pdatanum', 1, ... % Which PData structure to use
%                'rcvBufFrame', -1, ... % Use the most recently transferred frame
%                'IntBufDest', [1, -1], ... % IQ (complex) data, Auto-increment for every frame
%                'ImgBufDest', [1, -1], ... % [buffer #, frame #] Auto-increment ImageBuffer for each reconstruction???? % something is [first/oldest frame, last/newest frame]
%                'RINums', [1:2*na]); % The ReconInfo structure #(s). Each Recon must have its own unique set of ReconInfo #s
% 
% % Recon = repmat(Recon, 1, numFrames);
% % for nf = 1:numFrames
% %     Recon(nf).IntBufDest = [1, nf];
% %     Recon(nf).ImgBufDest = [1, nf];
% % end
% 
% ReconInfo = repmat(struct('mode', 'accumIQ_replaceIntensity', ... % reconstruct, and replace intensity data in ImageBuffer and IQ data in InterBuffer (see Table 12.4 in Tutorial)
%                    'txnum', 1, ...                 % TX structure to use
%                    'rcvnum', 1, ...                % RX structure to use
%                    'regionnum', 1), 1, 2*na);                % PData Region to process in
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % for nf = 1
% for n = 1:2*na % need to change this and above for more than 1 frame
%     % - Set specific ReconInfo attributes.
%     % ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
%     ReconInfo(n).txnum = n;
%     ReconInfo(n).rcvnum = n;
%     ReconInfo(n).pagenum = n; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     ReconInfo(1).regionnum = 1; %1 for the whole volume, 5 for the slices
% 
% end


%% Process structures

Process(1).classname = 'External';
Process(1).method = 'saveTimetag';
Process(1).Parameters = {'srcbuffer', 'none', ...
                             'dstbuffer', 'none'};
nprevproc = 1; % number of previous Processes

for nbuf = 1:numBuffers
    Process(nbuf + nprevproc).classname = 'External';
    Process(nbuf + nprevproc).method = 'saveRcvData_timetag'; % Function name
    % Process(1).Parameters = {'srcbuffer', 'bufferName', ...
    %                          'srcbufnum', 1, ... % # of buffer to process
    %                          'srcframenum', 1, ... % starting frame #
    %                          'srcsectionnum', 1, ...
    %                          ' srcpagenum', 1, ...
    %                          'dstbuffer', 'bufferName', ... % destination buf
    %                          'dstbufnum', 1, ...
    %                          'dstframenum', 1, ...
    %                          'dstsectionnum', 1, ...
    %                          'dstpagenum', 1};
    
    Process(nbuf + nprevproc).Parameters = {'srcbuffer', 'receive', ...
                             'srcbufnum', nbuf, ... % # of buffer to process
                             'dstbuffer', 'none'};
end

%%
makeParameterStructureSmall_ULM;
%% New Event structure

Resource.VDAS.dmaTimeout = 1000;
% Flow:
% 1. Transmit (TX)
% 2. Receive (Receive)
% 3. Reconstruction (Recon)
% 4. Processing (Process)
% 5. Control (SeqControl)

scInd = 1; % sequence control index
SeqControl(scInd).command = 'timeToNextAcq'; % In us, allowed range is from 10 - 4190000
                                         % Very useful if you are switching
                                         % the TPC (voltage) across acqs,
                                         % since it takes 800 us - 8 ms to
                                         % switch
% SeqControl(scInd).condition = 'ignore';  % don't print the warning message

timePerAcq = 1 / PRF * 1e6; % PRF in us

timePerAcqLimits = [10, 4190000];
if timePerAcq < timePerAcqLimits(1)
    warning('Acquisition time too short, setting to minimum of 10 us')
    SeqControl(scInd).argument = timePerAcqLimits(1); 
elseif timePerAcq > timePerAcqLimits(2)
    warning('Acquisition time too long, setting to maximum of 4190000 us')
    SeqControl(scInd).argument = timePerAcqLimits(2);
else
    SeqControl(scInd).argument = timePerAcq;
end

scInd = scInd + 1;
SeqControl(scInd).command = 'returnToMatlab';

scInd = scInd + 1;
SeqControl(scInd).command = 'jump'; % jump to
SeqControl(scInd).argument = 1;     % first event

% frame/volume rate
timePerFrame = SeqControl(scInd - 2).argument * na * 2;     % Time to acquire all the acquisitions for one frame/volume based on PRF (us)
frameTimeGap = 1 / frameRate * 1e6 - timePerFrame + SeqControl(scInd-2).argument;      % Add delays to account for the frame/volume rate set above. % Add delays to account for the frame/volume rate set above. Add the PRF time because this value replaces one of those delays too.

timePerTransferSeconds = numGBPerBufferFrame / 6.6; % PCie speed is 6.6 GB/s
frameTimeGapDMALimit = timePerTransferSeconds * 1e6; % [us]
frameRateDMALimit = 1/(frameTimeGapDMALimit * 1e-6);

% Check if the inputted frame rate is faster than what the DMA time allows
if frameRateDMALimit < frameRate
    warning("**** DMA time limits the frame rate to " + num2str(frameRateDMALimit) + " Hz, exiting ****")
    return
end

scInd = scInd + 1;
SeqControl(scInd).command = 'timeToNextAcq';

if frameTimeGap < timePerAcqLimits(1)
    warning('Frame delay time too short, setting to minimum of 10 us')
    SeqControl(scInd).argument = timePerAcqLimits(1); 
elseif frameTimeGap > timePerAcqLimits(2)
    warning('Frame delay time too long, setting to maximum of 4190000 us')
    SeqControl(scInd).argument = timePerAcqLimits(2);
else
    SeqControl(scInd).argument = frameTimeGap;
end

% frame/volume rate noop
scInd = scInd + 1;
SeqControl(scInd).command = 'noop';                     % no operation
frame_noop_time_us = SeqControl(scInd - 1).argument;
SeqControl(scInd).argument = frame_noop_time_us / 200 * 1e3;  % (value*200nsec; max. value is 2^25 - 1 for 6.7 sec)
SeqControl(scInd).condition = 'Hw&Sw';                  % need to enable the noop in hardware

% buffer rate

% need to change this to be consistent with the if blocks above
timePerBuffer = 1 / frameRate * numFramesPerSF * 1e6;                 % Time to acquire all the frames within one buffer (us)
bufferTimeGap = timePerBuffer / bufferDutyCycle - timePerBuffer;          % Add delay to account for the buffer rate duty cycle set above

scInd = scInd + 1;
SeqControl(scInd).command = 'timeToNextAcq';

if bufferTimeGap < timePerAcqLimits(1)
    warning('Buffer delay time too short, setting to minimum of 10 us')
    SeqControl(scInd).argument = timePerAcqLimits(1); 
elseif bufferTimeGap > timePerAcqLimits(2)
    warning('Buffer delay time too long, setting to maximum of 4190000 us')
    SeqControl(scInd).argument = timePerAcqLimits(2);
else
    SeqControl(scInd).argument = bufferTimeGap;
end

% buffer rate noop
scInd = scInd + 1;
SeqControl(scInd).command = 'noop'; % jump to
buffer_noop_time_us = SeqControl(scInd - 1).argument;
SeqControl(scInd).argument = buffer_noop_time_us / 200 * 1e3; % (value*200nsec; max. value is 2^25 - 1 for 6.7 sec)
SeqControl(scInd).condition = 'Hw&Sw'; % need to enable the noop in hardware

% Sync for aligning the hardware to when the data is done saving
scInd = scInd + 1;
SeqControl(scInd).command = 'sync';
SeqControl(scInd).argument = 30000000; % 30 s
% SeqControl(scInd).argument = 1000000 * vts.delay_s*5; % Timeout set to 5x the input delay just in case


n = 0;
for nbuf = 1:numBuffers
    for nf = 1:numFramesPerSF
    
        for a = 1:na % go through all the angles for each frame
            n = n + 1;
            Event(n).info = 'Transmit all columns and receive all rows';
            Event(n).tx = a.*2 - 1; % Use ath TX structure
            Event(n).rcv = (nbuf - 1) .* numFramesPerSF .* pair .* na + (nf - 1).*pair.*na + a.*2 - 1; % Use nth Receive structure % need to make this alternate between (1 and 2) * numframes or something
            Event(n).recon = 0; % 0 means no reconstruction
            Event(n).process = 0; % 0 means no processing
            Event(n).seqControl = 1;
            
            n = n + 1;
            Event(n).info = 'Transmit all rows and receive all columns';
            Event(n).tx = a.*2; 
            Event(n).rcv = (nbuf - 1) .* numFramesPerSF .* pair .* na + (nf - 1).*pair.*na + a.*2; 
            Event(n).recon = 0; 
            Event(n).process = 0; 
            Event(n).seqControl = 1;  
        
        end
        scInd = scInd + 1; 
        SeqControl(scInd).command = 'transferToHost'; % Transfer every frame
%         Event(n).seqControl = [4, 5, scInd];
%         Event(n).seqControl = [4, scInd];
            
        scInd = scInd + 1;
        SeqControl(scInd).command = 'waitForTransferComplete';
        SeqControl(scInd).argument = scInd - 1;
        Event(n).seqControl = [4, scInd - 1, scInd];

    end

    % Don't worry about timeToNextAcq while saving the superframe
    Event(n).seqControl = [scInd - 1, scInd];
%     Event(n).seqControl = [scInd];
%     Event(n).seqControl = [6, scInd];
    
    if saveRcvDataFlag
        n = n + 1;
        Event(n).info = 'Save data - ext proc func';
        Event(n).tx = 0; 
        Event(n).rcv = 0; 
        Event(n).recon = 0;
        Event(n).process = nbuf + nprevproc; 
        Event(n).seqControl = 8;
    end

end

% n = n + 1;
% 
% Event(n).info = 'Jump';
% Event(n).tx = 0; 
% Event(n).rcv = 0; 
% Event(n).recon = 0;
% Event(n).process = 0; 
% Event(n).seqControl = 3; 


% %% User specified UI Control Elements
% 
% import vsv.seq.uicontrol.VsSliderControl
% 
% % - Time Tag
% UI(1).Control = VsSliderControl('LocationCode', 'UserB5',...
%                                 'Label', 'Time Tag', ...
%                                 'SliderMinMaxVal', [0, 2, TimeTagEna],...
%                                 'SliderStep', [0.5, 0.5], ...
%                                 'ValueFormat', '%1.0f',...
%                                 'Callback', @TimeTagCallback);
% 
% 
% % External function definitions.
% 
% import vsv.seq.function.ExFunctionDef
% 
% EF(1).Function = vsv.seq.function.ExFunctionDef('readTimeTag',@readTimeTag);

%% Save all the data/structures to a .mat file.
currentDir = cd; currentDir = regexp(currentDir, filesep, 'split');
filename = 'RC15gV_Allen_singleframe.mat';

save(fullfile(currentDir{1:find(contains(currentDir,"Vantage"),1)})+"\MatFiles\"+filename);

% %% Initialize time tagging if enabled
% import com.verasonics.hal.hardware.*
% switch TimeTagEna
%     case 0
%         % disable time tag
%         rc = Hardware.enableAcquisitionTimeTagging(false);
%         if ~rc
%             error('Error from enableAcqTimeTagging')
%         end
%         tagstr = 'off';
%     case 1
%         % enable time tag
%         rc = Hardware.enableAcquisitionTimeTagging(true);
%         if ~rc
%             error('Error from enableAcqTimeTagging')
%         end
%         tagstr = 'on';
%         disp('**** Time tagging enabled on mode 1 ****')
%     case 2
%         % enable time tag and reset counter
%         rc = Hardware.enableAcquisitionTimeTagging(true);
%         if ~rc
%             error('Error from enableAcqTimeTagging')
%         end
%         rc = Hardware.setTimeTaggingAttributes(false, true); % reset hardware counter to 0 (otherwise, it continuously counts up from system bootup until it gets to 107,000s - see p37 of User Manual
%         if ~rc
%             error('Error from setTimeTaggingAttributes')
%         end
%         tagstr = 'on, reset';
%         disp('**** Time tagging enabled on mode 2 ****')
% end

%% Run VSX automatically and make parameter structure for RF file naming

if runVSX
    disp("running VSX")
    VSX
end

%% Save post-acquisition parameters in a structure P

% save([savepath, 'params.mat'], 'angles', 'startDepth', 'startDepthMM', 'endDepth', 'endDepthMM', 'Event', 'fps_target', 'maxAcqLength_adjusted', 'maxAngle', 'Media', 'na', 'nf', 'Receive', 'Resource', 'SeqControl', 'TGC', 'Trans', 'TW', 'TX', 'wl', 'numElements', '-v7.3')
% save([savepath, 'params.mat'], 'P')

makeParameterStructure_ULM;
savefast([savepath, 'params.mat'], 'P')
% saveRcvData(RcvData{1})
clearvars RcvData
% save([savepath, 'workspace.mat'], '-v7.3', '-nocompression')

%% **** Callback routines used by UIControls (UI) ****
% %% Time tag callback test
% 
% function TimeTagCallback(~, ~, UIValue)
%     import com.verasonics.hal.hardware.*
%     TimeTagEna = round(UIValue);
%     VDAS = evalin('base', 'VDAS');
%     switch TimeTagEna
%         case 0
%             if VDAS % can't execute this command if HW is not present
%                 % disable time tag
%                 rc = Hardware.enableAcquisitionTimeTagging(false);
%                 if ~rc
%                     error('Error from enableAcqTimeTagging')
%                 end
%             end
%             tagstr = 'off';
%         case 1
%             if VDAS
%                 % enable time tag
%                 rc = Hardware.enableAcquisitionTimeTagging(true);
%                 if ~rc
%                     error('Error from enableAcqTimeTagging')
%                 end
%             end
%             tagstr = 'on';
%         case 2
%             if VDAS
%                 % enable time tag and reset counter
%                 rc = Hardware.enableAcquisitionTimeTagging(true);
%                 if ~rc
%                     error('Error from enableAcqTimeTagging')
%                 end
%                 rc = Hardware.setTimeTaggingAttributes(false, true); % reset hardware counter to 0 (otherwise, it continuously counts up from system bootup until it gets to 107,000s - see p37 of User Manual
%                 if ~rc
%                     error('Error from setTimeTaggingAttributes')
%                 end
%             end
%             tagstr = 'on, reset';
%     end
%     % display at the GUI slider value
%     h = findobj('Tag', 'UserB5Edit');
%     set(h,'String', tagstr);
%     assignin('base', 'TimeTagEna', TimeTagEna);
% end
% 
% %% **** Callback routines used by External function definition (EF) ****
% 
% function readTimeTag(RDatain)
%     persistent frmCount
%     if isempty(frmCount)
%         frmCount = 0;
%     end
%     % get time tag from first two samples
%     % time tag is 32 bit unsigned interger value, with 16 LS bits in sample 1
%     % and 16 MS bits in sample 2.  Note RDatain is in signed INT16 format so must
%     % convert to double in unsigned format before scaling and adding
%     W = zeros(2, 1);
%     for i=1:2
%         W(i) = double(RDatain(i, 1));
%         if W(i) < 0
%             % translate 2's complement negative values to their unsigned integer
%             % equivalents
%             W(i) = W(i) + 65536;
%         end
%     end
%     timeStamp = W(1) + 65536 * W(2);
%     % the 32 bit time tag counter increments every 25 usec, so we have to scale
%     % by 25 * 1e-6 to convert to a value in seconds
%     frmCount = frmCount + 1;
%     if mod(frmCount, 25) == 1
%         TimeTagEna = evalin('base', 'TimeTagEna');
%         if TimeTagEna
%             disp(['Time tag value in seconds ', num2str(timeStamp/4e4,'%2.3f')]);
%         end
%     end
% end