% **** USE MULTIPLE BUFFERS **** %

warning('Need to fix data racing, aperture')
%% 0. Description
% Continuous acquisition and saving of RF data with the RC15gV probe
% CPWC, stacks all frames per superframe in one transfer/file
% C-R and R-C pairs of TX-RX
% Uses saveRcvData external function for saving
% Starts on an external trigger

% Collects a stack of nf frames with one transfer for all at once

% This script sets up multiple buffers for possibly-continuous acquisition:
%   1. Acquire one buffer
%   2. Acquire the next buffer while transferring and saving the previous buffer
%   3. Repeat step 2
%% 1. Specify system parameters
clearvars

codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Verasonics");
AcquisitionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Ultrasound data acquisition");
addpath(AcquisitionCodePath)

addpath('C:\Users\BOAS-US\Documents\GitHub\BU-Code\Allen code\Air Puff\')

cd 'C:\Users\BOAS-US\Desktop\Vantage-5.0.0-p1'
% cd 'C:\Users\BOAS-US\Desktop\Vantage-4.9.5-2409181500'
% cd 'G:\My Drive\Verasonics files\Vantage-4.9.2-2308102000'
activate

savepath = uigetdir('G:\', 'Select the save path');
savepath = [savepath, '\'];

parameterPrompt = {'Probe voltage [V]', 'Start depth [mm]', 'End depth [mm]', 'Pulse Repetition Frequency [Hz]', 'Frame rate [Hz]', 'Number of angles', 'Maximum angle [degrees]', 'Probe frequency [MHz]', 'Speed of sound [m/s]', 'Simulate Mode (0-off, 1-on, 2-RcvLoop)', 'Save RcvData (0-no, 1-yes)', 'Number of frames per superframe', 'Use air puff (0-no, 1-yes)', 'Probe connector', 'SSD write speed [GB/s]', 'Probe aperture [mm]', 'Time per superframe [s]', 'ADC Sampling Mode (50, 67, 100, 200% of center frequency)', 'Number of buffers'}; % 'Save RF data (0-no, 1-yes)', 
parameterDefaults = {'30', '0', '8', '60000', '2500', '11', '5', '13.6', '1540', '0', '1', '500', '0', 'UTA-408GE', '1.45', '8.8', '1.5', '200', '3'};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% ADC_sampleMode = 'BS67BW';
% spw_guess = 1.3333;

% ADC_sampleMode = 'BS100BW';
% spw_guess = 2;

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
if mod(numFramesPerSF, 2) ~= 0
    error('# of frames per SF must be even')
end
useTriggers = str2double(parameterUserInput{13});
connectorPlate = parameterUserInput{14};
% Maximum PCIe DMA rate for the Vantage 256 is 6.6 GB/s, but the connector
% type can affect this
DMARate = getDMARate(connectorPlate);
SSDWriteRate = str2double(parameterUserInput{15});
apertureMM = str2double(parameterUserInput{16});
% apertureMM = 12.8; % Use some subset of the probe elements [mm]
% apertureMM = 8;
TimePerSF = str2double(parameterUserInput{17});
sfRate = 1/TimePerSF; % Superframe rate [Hz]

ADC_sampleModeNumber = str2double(parameterUserInput{18}); % Options: 50, 67, 100, 200 (% of center frequency)
[ADC_sampleMode, spw_guess] = getADCSampleMode(ADC_sampleModeNumber);

numBuffers = str2double(parameterUserInput{19}); % Default of 3 buffers. The basic idea is to have two buffers so we can do continuous acquisition/saving, but adding a 3rd buffer to give us some *buffer room*

% tagtest = Hardware.enableAcquisitionTimeTagging(1);
bufferIndex = 0;
runVSX = 1;
movePointsOrNot = 0;
numChannels = 256; % enable channels

% Angles for plane waves are equally distributed over the defined range/# angles
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

% Resource.Parameters.waitForProcessing = 1;
%  the hardware will wait before each 'transferToHost' for the software to finish 
% processing the previous frame transferred to the host. When software asks for the most 
% recently transferred frame, it gets the frame previously transferred and can start 
% processing immediately.  This also releases the hardware to start transferring the 
% previously acquired frame and to start acquiring the next.  Acquisition, data transfer and 
% processing are then all occurring at the same time, but using different data - a process 
% known as pipelining.

%% 1.5. Specify the functional stimulus parameters
if useTriggers
    [apis, vts, daqrate, numTrials] = functionalParameterInputPrompt;
end

%% 2. Define Transducer structure

Trans.name = 'RC15gV'; 
Trans.frequency = probe_freq; % Not needed if using the default center frequency
Trans.units = 'wavelengths'; % or mm

Trans = computeTrans(Trans); % Generate required attributes for the probe into the Trans structure; e.g., the transducer element positions
% Trans.maxHighVoltage = ; % set maximum high voltage that is allowed to the transducer

L = Trans.spacingMm*Trans.numelements/2/1e3; % Probe width, in m
wl = Resource.Parameters.speedOfSound / Trans.frequency / 1e6; % Wavelength, in m

startDepth = startDepthMM/1e3/wl; % start depth in wavelengths
endDepth = endDepthMM/1e3/wl; % end depth in wavelengths

%% Set the active aperture
apertureTotalMM = Trans.spacingMm * (Trans.numelements/2); % Total available probe aperture [mm]
% apertureElem = floor(apertureMM/apertureTotalMM * Trans.numelements);
apertureElem = floor(apertureMM / Trans.spacingMm); % Number of active elements
if mod(apertureElem, 2) ~= 0, error('Number of active elements must be even'), end

%% Modify the angles to reduce grating lobes (see Sauvage et al., 2020)
% angpitch = wl / (Trans.spacingMm*Trans.numelements / 2 / 1e3);
% angles = -(na - 1) / 2 * angpitch : angpitch : (na - 1) / 2 * angpitch
% maxAngle = max(angles);
% warning('**** MODIFYING THE ANGULAR PITCH AND RANGE TO MINIMIZE GRATING LOBES ****')
% disp('Actual maxAngle [deg]: ')
% disp(max(angles * 180 / pi))

%% enable time tag
TimeTagEna = 1;
% 0: disable
% 1: enable but don't reset counter
% 2: enable and reset counter

%% Simulation things - Media structure (define scattering points and attenuation)
Resource.Parameters.simulateMode = simMode; % run script in simulate mode. Set to 0 if not

% xd_mm = 5; % in mm
% xd = xd_mm/wl/1e3;

% Set up Media model for the simulation, which generates the scattering points with 3D location
% and reflectivity. For 1D transducer arrays, they are aligned on the
% x-axis with the center at x = 0, and scan depth is in z.
Media.MP(1, :) = [0, 0, 50, 1.0]; % [x, y, z, reflectivity]. x, y, z are defined as # of wavelengths.
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

% Media.attenuation = 0;
Media.attenuation = -0.7; % media attenuation in dB/cm/MHz

Media.function = 'movePointsZ3D'; % move points in _ dimension after each frame


%% Transmission Waveform (TW)
tw.A = Trans.frequency; % frequency of transmission pulse, sets half cycle period of the waveform...
tw.B = 0.67; % amount of time (0.1 - 1.0) that the transmission drivers are active in the half cycle period. Controsl how much power is delivered.
             % Apparently using B = 0.67 approximates a sine wave.
% tw.C = 2; % number of half cycles in the transmission waveform. 2 half cycles = 1 full cycle burst
tw.C = 3; % number of half cycles in the transmission waveform. 2 half cycles = 1 full cycle burst
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

% Define transmit element apodization (for the rows or columns)
activeElem = kaiser(apertureElem).';
emitElem = [zeros(1, ((Trans.numelements/2) - apertureElem)/2), activeElem, zeros(1, ((Trans.numelements/2) - apertureElem)/2)];

% na*2 transmissions of a plane wave in pairs, one by all row elements and then one by all
% column elements
TX = repmat(struct('waveform', 1, ...
                   'focus', 0, ... % plane wave
                   'Steer', [0.0, 0.0], ... % theta, alpha (beam angle projected in xz from +z axis, beam angle wrt xz)
                   'Apod', zeros(1, Trans.numelements)), 1, na*2); % Initialize apod to all zeros, set the active ones below
for n = 1:na
    % TX(n).Apod(1:Trans.numelements/2) = ones(1, Trans.numelements/2); % Turn on columns (y)
    TX(n).Apod(1:Trans.numelements/2) = emitElem; % Turn on columns (y)
    TX(n).Steer = [angles(n), 0];
    TX(n).Delay = computeTXDelays(TX(n));
end

for n = 1:na
    % TX(na + n).Apod(Trans.numelements/2 + 1 : end) = ones(1, Trans.numelements/2); % Turn on rows (x)
    TX(na + n).Apod(Trans.numelements/2 + 1 : end) = emitElem; % Turn on rows (x)
    TX(na + n).Steer = [0, angles(n)];
    TX(na + n).Delay = computeTXDelays(TX(na + n));
end

%% Define Time Gain Control waveform (TGC)
% Accounts for decrease in amplitude of echoes for longer distance traveled

% TGC curve definition
% TGC.CntrlPts = [0 785.2216 1023 1023 1023 1023 1023 1023];
% TGC.CntrlPts = [1023 1023 1023 1023 1023 1023 1023 1023];
TGC.CntrlPts = [590,650,710,770,830,890,950,1010];
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

%% RcvProfile adjustment (8/7/25 change)
RcvProfile.antiAliasCutoff = 20; % Low pass filter at 20 MHz (RC15gV bandwidth goes to 19 MHz)
% RcvProfile.LnaZinSel = 25;

%% Receiver array object

% Define receive element apodization (for the rows or columns)
rcvElem = [zeros(1, ((Trans.numelements/2) - apertureElem)/2), activeElem, zeros(1, ((Trans.numelements/2) - apertureElem)/2)];

maxAcqLength = ceil(sqrt(endDepth^2 + 2*(numElements*Trans.spacing)^2)); % account for the longest distance an echo could travel
Receive = repmat(struct('Apod', zeros(1, Trans.numelements), ... % Initialize all elements to 0 apod initially, set below
                        'startDepth', startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ... % This is changed in the loop below
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', ADC_sampleMode, ...
                        'mode', 0, ...
                        'callMediaFunc', 0, ...
                        'LowPassCoef', [], ...
                        'InputFilter', []), 1, pair * numFramesPerSF * na * numBuffers);
j = 1;
% an = 0;
for nbuf = 1:numBuffers
    an = 0; % acquisition number

    for nf = 1:numFramesPerSF
%         an = 0; % acquisition number
        
        % Move points after all the acquisitions for one frame
        Receive(j).callMediaFunc = movePointsOrNot;
    %     Receive(j).mode = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for n = 1:na
            an = an + 1;
            Receive(j).bufnum = nbuf;
%             Receive(j).framenum = nf;
%             Receive(j).framenum = nbuf;
            Receive(j).framenum = 1; % 1 stacked frame per buffer
            Receive(j).acqNum = an;
            % Receive(j).Apod(Trans.numelements/2 + 1 : end) = ones(1, Trans.numelements/2);
            Receive(j).Apod(Trans.numelements/2 + 1 : end) = rcvElem;
            j = j + 1;
        end
    
        for n = 1:na
            an = an + 1;
            Receive(j).bufnum = nbuf;
%             Receive(j).framenum = nf;
%             Receive(j).framenum = nbuf;
            Receive(j).framenum = 1; % 1 chunked frame per buffer
            Receive(j).acqNum = an;
            % Receive(j).Apod(1:Trans.numelements/2) = ones(1, Trans.numelements/2);
            Receive(j).Apod(1:Trans.numelements/2) = rcvElem;
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
    samplesPerWave = spw_guess;
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
    Resource.RcvBuffer(nbuf).rowsPerFrame = numRcvSamples * na * 2 * numFramesPerSF;
    Resource.RcvBuffer(nbuf).colsPerFrame = Resource.Parameters.numRcvChannels; % Usually 1:1 to # of receive channels available in the system. Can change to 256 with the 2D probe and new connector plate.
%     Resource.RcvBuffer(nbuf).colsPerFrame = 160; % Usually 1:1 to # of receive channels available in the system. Can change to 256 with the 2D probe and new connector plate.
%     Resource.RcvBuffer(nbuf).numFrames = numFramesPerSF; % minimum # frames of RF data to acquire; RcvBuffer contains all the data needed for a whole frame, including multiple acquisition passes needed for reconstruction. Software can re-process RcvBuffer frames
    Resource.RcvBuffer(nbuf).numFrames = 1; % minimum # frames of RF data to acquire; RcvBuffer contains all the data needed for a whole frame, including multiple acquisition passes needed for reconstruction. Software can re-process RcvBuffer frames
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

numSamplesPerSubFrame = Resource.RcvBuffer(1).rowsPerFrame / numFramesPerSF * Resource.RcvBuffer(1).colsPerFrame; % This should be all frames within a superframe
numSamplesInBuffer = Resource.RcvBuffer(1).rowsPerFrame * Resource.RcvBuffer(1).colsPerFrame * Resource.RcvBuffer(1).numFrames
numGBInBuffer = numSamplesInBuffer ./ 1024^3 * 2 % # samples * (2 bytes per int16 sample) 
numSamplesPerBufferFrame = Resource.RcvBuffer(1).rowsPerFrame * Resource.RcvBuffer(1).colsPerFrame
numGBPerBufferFrame = numSamplesPerBufferFrame ./ 1024^3 * 2 % # samples * (2 bytes per int16 sample) 

if numGBPerBufferFrame > 2
    warning('Buffer size per frame is too large (> 2 GB), exiting')
    return
end

% if ((maxAcqLength_adjusted + (endDepth-startDepth))*wl / speedOfSound) > 1/PRF
% if ((maxAcqLength_adjusted + (endDepth))*wl / speedOfSound) > 1/PRF
if ((2*maxAcqLength_adjusted)*wl / speedOfSound) > 1/PRF
    error('Error: the PRF is too high, it will send the next transmission before the previous transmission reflects from the deepest part of the region')
end

%% Check for if the superframe size and superframe rate are incompatible
DMATime = numGBInBuffer/DMARate; % Time to DMA one superframe [s]
AcqTime = numFramesPerSF / frameRate; % Time to acquire on superframe [s]
SFTime = 1/sfRate; % Nominal time per superframe [s]
SSDWriteTime = numGBInBuffer/SSDWriteRate; % Theoretical time to write one superframe to disk [s]

disp("======== Timing checks ========")
disp("Time to acquire: " + num2str(AcqTime) + "s")
disp("Time to DMA and write to disk: " + num2str(DMATime + SSDWriteTime) + "s")

% disp("Time to acquire and DMA: " + num2str(AcqTime + DMATime) + "s")
% disp("Time to acquire, DMA, and write to disk: " + num2str(AcqTime + DMATime + SSDWriteTime) + "s")
numFramesPerSFPossible = SFTime / [numSamplesPerSubFrame / 1024^3 * 2 *(1/DMARate + 1/SSDWriteRate)]; % # of possible frames per superframe under these acquisition settings
if (SFTime - DMATime) < 0
    error("There is not enough time to acquire and DMA the amount of frames specified, according to the superframe rate")
elseif (SFTime - DMATime - SSDWriteTime) < 0
    disp("**** Number of frames possible with the set acquisition rates: " + num2str(numFramesPerSFPossible) + "****")
    error("There is not enough time to acquire, DMA, and save to disk the amount of frames specified, according to the superframe rate")
end

disp("===============================")

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

%% Store a select couple of parameters into a structure for updating the save data's filename
makeParameterStructureSmall_functional;

%% Event structure

Resource.VDAS.dmaTimeout = 100000; % [ms]

% Set the shot-to-shot (each angle) timing according to the PRF
scInd = 1; % sequence control index
SeqControl(scInd).command = 'timeToNextAcq'; % In us, allowed range is from 10 - 4190000
                                         % Very useful if you are switching
                                         % the TPC (voltage) across acqs,
                                         % since it takes 800 us - 8 ms to
                                         % switch
SeqControl(scInd).condition = 'ignore';  % don't print the warning message

timePerAcq = 1 / PRF * 1e6; % time step according to the PRF [us]

timePerAcqLimits = [10, 4190000]; % Verasonics hardware limits for timetoNextAcq [us]
if timePerAcq < timePerAcqLimits(1)
    warning('Shot acquisition time too short, setting to minimum of 10 us')
    SeqControl(scInd).argument = timePerAcqLimits(1); 
elseif timePerAcq > timePerAcqLimits(2)
    warning('Shot acquisition time too long, setting to maximum of 4190000 us')
    SeqControl(scInd).argument = timePerAcqLimits(2);
else
    SeqControl(scInd).argument = timePerAcq;
end

% 2. Return to Matlab SeqControl
scInd = scInd + 1;
SeqControl(scInd).command = 'returnToMatlab';

% 3. Jump to some event to keep the acquisition looping
scInd = scInd + 1;
SeqControl(scInd).command = 'jump'; % jump to
if useTriggers
    SeqControl(scInd).argument = 2;     % second event
else
    SeqControl(scInd).argument = 1;     % first event
end
SeqControl(scInd).condition = 'exitAfterJump'; % Normally, jumping auto returns to Matlab if it returns to the first event, but not for other events

% 4. Set the frame/volume rate
timePerFrame = SeqControl(scInd-2).argument * na * 2;     % Time to acquire all the acquisitions for one frame/volume based on the PRF [us]
% frameTimeGap = 1 / frameRate * 1e6 - timePerFrame;      % Add delays to account for the frame/volume rate set above
frameTimeGap = 1 / frameRate * 1e6 - timePerFrame + SeqControl(scInd-2).argument;      % Add delays to account for the frame/volume rate set above. Add the PRF time because this value replaces one of those delays too.

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

% 5. frame/volume rate noop (deprecated)
scInd = scInd + 1;
SeqControl(scInd).command = 'noop';                     % no operation
frame_noop_time_us = SeqControl(scInd - 1).argument;
SeqControl(scInd).argument = frame_noop_time_us / 200 * 1e3;  % (value*200nsec; max. value is 2^25 - 1 for 6.7 sec)
SeqControl(scInd).condition = 'Hw&Sw';                  % need to enable the noop in hardware

% 6. buffer rate (deprecated)

% need to change this to be consistent with the if blocks above
timePerBuffer = 1 / frameRate * numFramesPerSF * 1e6;                 % Time to acquire all the frames within one buffer (us)
% bufferTimeGap = timePerBuffer / bufferDutyCycle - timePerBuffer;          % Add delay to account for the buffer rate duty cycle set above

scInd = scInd + 1;
SeqControl(scInd).command = 'timeToNextAcq';
SeqControl(scInd).argument = timePerAcqLimits(2); % dummy statement
% if bufferTimeGap < timePerAcqLimits(1)
%     warning('Buffer delay time too short, setting to minimum of 10 us')
%     SeqControl(scInd).argument = timePerAcqLimits(1); 
% elseif bufferTimeGap > timePerAcqLimits(2)
%     warning('Buffer delay time too long, setting to maximum of 4190000 us')
%     SeqControl(scInd).argument = timePerAcqLimits(2);
% else
%     SeqControl(scInd).argument = bufferTimeGap;
% end

% 7. buffer rate noop (deprecated)
scInd = scInd + 1;
SeqControl(scInd).command = 'noop';
buffer_noop_time_us = SeqControl(scInd - 1).argument;
SeqControl(scInd).argument = buffer_noop_time_us / 200 * 1e3; % (value*200nsec; max. value is 2^25 - 1 for 6.7 sec)
SeqControl(scInd).condition = 'Hw&Sw'; % need to enable the noop in hardware

% 8. Trigger input
scInd = scInd + 1;
SeqControl(scInd).command = 'triggerIn';
SeqControl(scInd).argument = 0; % 0-255. Each increment of 1 corresponds to 250 ms. The default is 0 and means to wait indefinitely.
SeqControl(scInd).condition = 'Trigger_2_Rising'; % Which trigger in port and type to use
% SeqControl(scInd).command = 'pause';
% SeqControl(scInd).argument = 19; % see p137
% SeqControl(scInd).condition = 'extTrigger'; % need to enable the noop in hardware

% 9. Trigger output
% "Generates external 1 microsecond active low output on the TRIG OUT BNC
%  connector. A delay can be set in the argument field." (p138)
scInd = scInd + 1;
SeqControl(scInd).command = 'triggerOut';
% SeqControl(scInd).argument = 0; % 0-255. Each increment of 1 corresponds to 250 ms. The default is 0 and means to wait indefinitely.
SeqControl(scInd).condition = 'syncNone'; % syncNone -> generate the trigger asap after the scheduled time

% 10. Sync to make the software sequencer also wait for the trigger input
scInd = scInd + 1;
SeqControl(scInd).command = 'sync';
SeqControl(scInd).argument = 10000000; % 10 s

% 11. Sync for aligning the hardware to when the data is done saving
scInd = scInd + 1;
SeqControl(scInd).command = 'sync';
if useTriggers
    SeqControl(scInd).argument = 1e6 * vts.delay_s*5; % Timeout set to 5x the input delay just in case
else
    SeqControl(scInd).argument = 20*1e6; % 20 s
end

% 12. Control superframe rate
scInd = scInd + 1;
SeqControl(scInd).command = 'timeToNextAcq';
SeqControl(scInd).argument = (1/sfRate - 1/frameRate*numFramesPerSF) * 1e6 + SeqControl(1).argument; % [us]

if useTriggers
    n = 1;
    Event(n).info = 'Wait for external trigger to start the acquisition sequence';
    Event(n).tx = 1; % It seems to not work properly if there isn't some acquisition event combined here
    Event(n).rcv = 0; 
    Event(n).recon = 0;
    Event(n).process = 1; % save the initial timetag
    Event(n).seqControl = [8, 11];
else
    n = 0;
end

% for nbuf = 1
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
%             Event(n).seqControl = 11;
%             Event(n).seqControl = [1, 11];

%               if mod(n, 90) == 0 & n > 0
%                 scInd = scInd + 1; 
%                 SeqControl(scInd).command = 'transferToHost'; % sub-DMA
%                 Event(n).seqControl = [1, scInd];
%               end

            n = n + 1;
            Event(n).info = 'Transmit all rows and receive all columns';
            Event(n).tx = a*2; 
            Event(n).rcv = (nbuf - 1) .* numFramesPerSF .* pair .* na + (nf - 1).*pair.*na + a*2; 
            Event(n).recon = 0; 
            Event(n).process = 0; 
            Event(n).seqControl = 1;
%             Event(n).seqControl = 11;
%             Event(n).seqControl = [1, 11];

%               if mod(n, 90) == 0 & n > 0
% %               if mod(n, n) == 0 & n > 0
%                 scInd = scInd + 1; 
%                 SeqControl(scInd).command = 'transferToHost'; % sub-DMA
%                 Event(n).seqControl = [1, scInd];
%               end
      
        end

        Event(n).seqControl = [4]; % set the frame rate control
%         scInd = scInd + 1; 
%         SeqControl(scInd).command = 'transferToHost'; % sub-DMA
%         Event(n).seqControl = [1, scInd];

%         % Transfer the previously acquired frame
%         scInd = scInd + 1; 
%         SeqControl(scInd).command = 'transferToHost'; % Transfer every frame
% %         Event(n).seqControl = [4, 5, scInd]; % includes some noop
% %         Event(n).seqControl = [4, scInd];
% 
%         % includes the waitForTransferComplete
%         scInd = scInd + 1;
%         SeqControl(scInd).command = 'waitForTransferComplete';
%         SeqControl(scInd).argument = scInd - 1;
%         Event(n).seqControl = [4, scInd - 1, scInd];

    end

    % Transfer the previously acquired frame
    scInd = scInd + 1; 
    SeqControl(scInd).command = 'transferToHost'; % Transfer every frame
%         Event(n).seqControl = [4, 5, scInd]; % includes some noop
%         Event(n).seqControl = [4, scInd];

    % Pause only the software sequencer, until the DMA is complete
    scInd = scInd + 1;
    SeqControl(scInd).command = 'waitForTransferComplete';
    SeqControl(scInd).argument = scInd - 1;
%     Event(n).seqControl = [4, scInd - 1, scInd];
    Event(n).seqControl = [12, scInd - 1, scInd]; % Superframe rate control

    if saveRcvDataFlag
        n = n + 1;
    
        Event(n).info = 'Save data - ext proc func';
        Event(n).tx = 0; 
        Event(n).rcv = 0; 
        Event(n).recon = 0;
        Event(n).process = nbuf + nprevproc; 
%         Event(n).seqControl = 7; 
        Event(n).seqControl = 0; 
        % Event(n).seqControl = 11;  % Sync hardware and software... Don't know if this is necessary

        n = n + 1;
    
        Event(n).info = 'Make sure the data is done saving before transferring the next - reset the waitForTransferComplete flag';
        Event(n).tx = 0; 
        Event(n).rcv = 0; 
        Event(n).recon = 0;
        Event(n).process = 0;
        scInd = scInd + 1;
        SeqControl(scInd).command = 'markTransferProcessed'; % Clear the waitForTransferComplete flag once the data for one buffer is saved
        SeqControl(scInd).argument = scInd - 2; % Refer to the last transferToHost command
        % Event(n).seqControl = [scInd, 11];
        Event(n).seqControl = [scInd]; % Do not sync, since we want asynchronous hardware and software sequencers to use the multiple buffer method
    end

end

n = n + 1;
Event(n).info = 'Jump';
Event(n).tx = 0; 
Event(n).rcv = 0; 
Event(n).recon = 0;
Event(n).process = 0; 
Event(n).seqControl = 3; 
% Event(n).seqControl = [3, 11]; 

% Add trigger out to the first frame within a superframe or buffer group
% Event(2).seqControl = [1, 9];

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
filename = 'RC15gV_Allen_loop_functional_buffers.mat';

save(fullfile(currentDir{1:find(contains(currentDir,"Vantage"),1)})+"\MatFiles\"+filename);

%% Run the air puff script before running VSX
if useTriggers
    [Mcr_d, Mcr_fcp] = controlAirPuff_func(apis, vts, daqrate, numTrials); % Need to use Mcr_ because VSX will autoclear most variables
    daqStartTimetag = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
    savefast([savepath, 'daqStartTimetag'], 'daqStartTimetag')
end

%% Initialize time tagging if enabled
import com.verasonics.hal.hardware.*
switch TimeTagEna
    case 0
        % disable time tag
        rc = Hardware.enableAcquisitionTimeTagging(false);
        if ~rc
            error('Error from enableAcqTimeTagging')
        end
        tagstr = 'off';
    case 1
        % enable time tag
        rc = Hardware.enableAcquisitionTimeTagging(true);
        if ~rc
            error('Error from enableAcqTimeTagging')
        end
        tagstr = 'on';
        disp('**** Time tagging enabled on mode 1 ****')
    case 2
        % enable time tag and reset counter
        rc = Hardware.enableAcquisitionTimeTagging(true);
        if ~rc
            error('Error from enableAcqTimeTagging')
        end
        rc = Hardware.setTimeTaggingAttributes(false, true); % reset hardware counter to 0 (otherwise, it continuously counts up from system bootup until it gets to 107,000s - see p37 of User Manual
        if ~rc
            error('Error from setTimeTaggingAttributes')
        end
        tagstr = 'on, reset';
        disp('**** Time tagging enabled on mode 2 ****')
end

%% Run VSX automatically and make parameter structure for RF file naming

if runVSX
    disp("running VSX")
    VSX
end

%% Read the air puff data - may need to put this in the saveRcvData Processing...
if useTriggers
    [inScanData, timeStamp, triggerTime] = read(Mcr_d, seconds(Mcr_d.NumScansAvailable / Mcr_d.Rate), "OutputFormat", "Matrix");
end

%% Save post-acquisition parameters in a structure P

makeParameterStructure_functional;
save([savepath, 'params.mat'], 'P')
if useTriggers
    save([savepath, 'triggerData.mat'], 'inScanData', 'timeStamp', 'triggerTime')
end
% saveRcvData(RcvData{1})
clear RcvData
% save([savepath, 'workspace.mat'], '-v7.3', '-nocompression')

%% **** Callback routines used by UIControls (UI) ****

%% Time tag callback test

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

%% **** Callback routines used by External function definition (EF) ****

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