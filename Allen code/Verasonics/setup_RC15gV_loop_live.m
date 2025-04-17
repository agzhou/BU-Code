
%% Description
% Continuous acquisition and saving of RF data with the RC15gV probe
% Saves all defined # of frames in one file
% C-R and R-C pairs of TX-RX
% Uses saveRcvData external function for saving
% Note: update the savepath variable as needed

%% Specify system parameters
clear

cd 'C:\Users\BOAS-US\Desktop\Vantage-4.9.5-2409181500'
% cd 'G:\My Drive\Verasonics files\Vantage-4.9.2-2308102000\Allen code\12-18-2024 RCA acq testing'

activate
% numElements = 80;

savepath = "G:\Allen\Data\04-15-2025 AZ02 fUS RC15gV\run 0\";
savepath = char(savepath);
mkdir(savepath)

ReconRegion = 1;

% tagtest = Hardware.enableAcquisitionTimeTagging(1);
supFrameIndex = 0;

runVSX = 1;
simOrNot = 0;
movePointsOrNot = 0;

initialVoltage = 20; % V

startDepthMM = 0; % start depth in mm
endDepthMM = 10;

fps_target = 50;   % Intended (sub)frame rate
supFrameBurstRate = .5; % Defines spacing between end of superframe burst and the next burst after jumping

numChannels = 256; % enable channels

numSupFrames = 10; % # of superframes, MUST BE ONE OR EVEN FOR VSX
numSubFrames = 1; % # of subframes
na = 21; % # of acquisitions per frame (acquisition pairs)
maxAngle = 5; % degrees
angleRange = [-maxAngle, maxAngle].*pi/180; % Angle range in radians

% Need at least 2 acquisitions to use multiple angles. 
% Otherwise, set angle to 0 degrees.
if na >= 2 
    angles = linspace(angleRange(1), angleRange(2), na);
else
    angles = 0;
end

pair = 2; % The R-C and C-R pair of acquisitions per angle

% Resource is a structure, define system parameters
Resource.Parameters.numTransmit = numChannels; % number of transmit channels
Resource.Parameters.numRcvChannels = numChannels; % number of receive channels
% Resource.Parameters.connector = 1; % transducer connector to use since the current plate for the 256 bit system is split into two 128 bit connectors. 1 is left and 2 is right
Resource.Parameters.speedOfSound = 1540; % speed of sound in m/s, the 1540 is for average human tissue
% Resource.Parameters.speedOfSound = 1481;
Resource.Parameters.verbose = 2; % Describe errors in varying levels

%% Define Transducer

Trans.name = 'RC15gV'; 
% Trans.frequency = 18.5; % Not needed if using the default center frequency
% Trans.frequency = 15.625;

Trans.units = 'wavelengths'; % or mm
% Trans.units = 'mm';

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
Resource.Parameters.simulateMode = simOrNot; % run script in simulate mode. Set to 0 if not

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
vesselX = 30e-6;    % x dimension
vesselY = 30e-6;    % x dimension
vesselZ = endDepthMM * 1e-3; % z dimension
flow_v_mm_s = 300;
%  
% Media.MP = randomPts3D_func(vesselX, vesselY, vesselZ, wl);

% Media.attenuation = 0;
Media.attenuation = -0.7; % media attenuation in dB/cm/MHz

Media.function = 'movePointsZ3D'; % move points in _ dimension after each frame

%% PData structure (Pixel Data --> image reconstruction range)
% For 2D scans and slices of 3D scans, it's always a rectangular area at a
% fixed location in the transducer coord system
xyplane = 50;

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

% PData.Region(2).Shape = struct('Name', 'Slice', 'Orientation', 'xz', ...
%                             'oPAIntersect', PData.Origin(2) - (numElements-1).*Trans.spacing./2); % out of Plane Axis Intersection
% PData.Region(3).Shape = struct('Name', 'Slice', 'Orientation', 'yz', ...
%                             'oPAIntersect', PData.Origin(1) + (numElements-1).*Trans.spacing./2);
% PData.Region(4).Shape = struct('Name', 'Slice', 'Orientation', 'xy', ...
%                             'oPAIntersect', Media.MP(3)); % currently set to the plane intersecting the only scatter point
PData(1).Region(2) = struct('Shape',struct('Name','Slice','Orientation','yz','oPAIntersect',PData.Origin(1)+PData.PDelta(2)*(PData.Size(2)-1)/2,'andWithPrev',1));
PData(1).Region(3) = struct('Shape',struct('Name','Slice','Orientation','xz','oPAIntersect',PData.Origin(2)-PData.PDelta(1)*((PData.Size(1)-1)/2),'andWithPrev',1));
PData(1).Region(4) = struct('Shape',struct('Name','Slice','Orientation','xy','oPAIntersect',xyplane,'andWithPrev',1));

PData.Region = computeRegions(PData);

PData.Region(5).PixelsLA = unique([PData.Region(2).PixelsLA; PData.Region(3).PixelsLA; PData.Region(4).PixelsLA]);
PData.Region(5).Shape.Name = 'Custom';
PData.Region(5).numPixels = length(PData.Region(5).PixelsLA);

%     'Position', [0, 0, 10], ...
%                      'width', PData.Size(2), 'height', PData.Size(1)./2);
                    % Position is relative to the global coords

% Define Super Region coordinates in 'P' structure in units of wavelengths
P.zCoord = (PData(1).Origin(3) + (0:PData(1).Size(3)-1)*PData(1).PDelta(3));
P.xCoord = (PData(1).Origin(1) + (0:PData(1).Size(2)-1)*PData(1).PDelta(1));
P.yCoord = (PData(1).Origin(2) - (0:PData(1).Size(1)-1)*PData(1).PDelta(2));

RegionIndices = reshape(PData.Region(1).PixelsLA,PData.Size);

RegXY = reshape(PData.Region(4).PixelsLA,[PData.Size(1),PData.Size(2),1]);
RegYZ = reshape(PData.Region(2).PixelsLA,[PData.Size(1),1,PData.Size(3)]);
RegXZ = reshape(PData.Region(3).PixelsLA,[1,PData.Size(2),PData.Size(3)]);

find(squeeze(any(bsxfun(@eq,RegionIndices,RegXY),[1,2])))
find(squeeze(any(bsxfun(@eq,RegionIndices,RegXZ),[2,3])))
find(squeeze(any(bsxfun(@eq,RegionIndices,RegYZ),[1,3])))

% Display window

xd = 70;

Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).Title = 'XZ plane';
Resource.DisplayWindow(1).pdelta = 0.4;
Resource.DisplayWindow(1).Position = [0,40, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).Orientation = 'xz';
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),PData(1).Region(3).Shape.oPAIntersect,0.0];
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).numFrames = numSubFrames;

Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).Title = 'YZ plane';
Resource.DisplayWindow(2).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(2).Position = [430,40, ...
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
Resource.DisplayWindow(2).Orientation = 'yz';
Resource.DisplayWindow(2).ReferencePt = [PData(1).Region(2).Shape.oPAIntersect,-PData(1).Origin(2),0];
Resource.DisplayWindow(2).Colormap = gray(256);
Resource.DisplayWindow(2).AxesUnits = 'mm';
Resource.DisplayWindow(2).numFrames = numSubFrames;

Resource.DisplayWindow(3).Type = 'Verasonics';
Resource.DisplayWindow(3).Title = 'XY plane';
Resource.DisplayWindow(3).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(3).Position = [860,40, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(3).pdelta), ... % width
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(3).pdelta)];    % height
Resource.DisplayWindow(3).Orientation = 'xy';
Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),-PData(1).Origin(2),PData.Region(4).Shape.oPAIntersect];
Resource.DisplayWindow(3).Colormap = gray(256);
Resource.DisplayWindow(3).AxesUnits = 'mm';
Resource.DisplayWindow(3).numFrames = numSubFrames;

% % xz
% Resource.DisplayWindow(1).Type = 'Verasonics';
% Resource.DisplayWindow(1).Title = 'Slice xz plane';
% Resource.DisplayWindow(1).pdelta = 0.3; % pixel spacing (in wavelengths) on the display window, for all dimensions
% llx = 100; % lower left corner x on screen
% xmult = 200;
% lly = 150; % lower left corner y
% Resource.DisplayWindow(1).Position = [llx, lly, ...
%                                       ceil(PData.Size(2).* PData.PDelta(1) ./ Resource.DisplayWindow(1).pdelta), ... % width (x)
%                                       ceil(PData.Size(3).* PData.PDelta(3) ./ Resource.DisplayWindow(1).pdelta)]; % height (z)
% Resource.DisplayWindow(1).ReferencePt = [PData.Origin(1), 0, PData.Origin(3)]; % Display Window location wrt transducer coords
% Resource.DisplayWindow(1).AxesUnits = 'mm'; % can change to mm
% Resource.DisplayWindow(1).Colormap = gray(256);
% Resource.DisplayWindow(1).Orientation = 'xz';
% Resource.DisplayWindow(1).numFrames = numSubFrames; % Define buffer size for a history of displayed frames
% 
% % xy
% Resource.DisplayWindow(2).Type = 'Verasonics';
% Resource.DisplayWindow(2).Title = 'Slice xy plane';
% Resource.DisplayWindow(2).pdelta = 0.3; % pixel spacing (in wavelengths) on the display window, for all dimensions
% 
% Resource.DisplayWindow(2).Position = [llx + 2.*xmult, lly, ...
%                                       ceil(PData.Size(2).* PData.PDelta(1) ./ Resource.DisplayWindow(1).pdelta), ... % width (x)
%                                       ceil(PData.Size(1).* PData.PDelta(2) ./ Resource.DisplayWindow(1).pdelta)]; % height (z)
% Resource.DisplayWindow(2).ReferencePt = [PData.Origin(1), -PData.Origin(2), xd]; % Display Window location wrt transducer coords
% Resource.DisplayWindow(2).AxesUnits = 'mm'; % can change to mm
% Resource.DisplayWindow(2).Colormap = gray(256);
% Resource.DisplayWindow(2).Orientation = 'xy';
% Resource.DisplayWindow(2).numFrames = numSubFrames; % Define buffer size for a history of displayed frames
% 
% % yz
% Resource.DisplayWindow(3).Type = 'Verasonics';
% Resource.DisplayWindow(3).Title = 'Slice yz plane';
% Resource.DisplayWindow(3).pdelta = 0.3; % pixel spacing (in wavelengths) on the display window, for all dimensions
% 
% Resource.DisplayWindow(3).Position = [llx + 4.*xmult, lly, ...
%                                       ceil(PData.Size(1).* PData.PDelta(2) ./ Resource.DisplayWindow(1).pdelta), ... % width (x)
%                                       ceil(PData.Size(3).* PData.PDelta(3) ./ Resource.DisplayWindow(1).pdelta)]; % height (z)
% Resource.DisplayWindow(3).ReferencePt = [0, -PData.Origin(2), PData.Origin(3)]; % Display Window location wrt transducer coords
% Resource.DisplayWindow(3).AxesUnits = 'mm'; % can change to mm
% Resource.DisplayWindow(3).Colormap = gray(256);
% Resource.DisplayWindow(3).Orientation = 'yz';
% Resource.DisplayWindow(3).numFrames = numSubFrames; % Define buffer size for a history of displayed frames


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
                        'InputFilter', []), 1, pair*numSupFrames*numSubFrames*na);
j = 1;
% an = 0;
for nsupf = 1:numSupFrames
    an = 0;
    for nsubf = 1:numSubFrames
        % Move points after all the acquisitions for one frame
        Receive(j).callMediaFunc = movePointsOrNot;
    %     Receive(j).mode = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for n = 1:na
            an = an + 1;
            Receive(j).framenum = nsupf;
            Receive(j).acqNum = an;
            Receive(j).Apod(Trans.numelements/2 + 1 : end) = ones(1, Trans.numelements/2);
            j = j + 1;
        end
    
        for n = 1:na
            an = an + 1;
            Receive(j).framenum = nsupf;
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

Resource.RcvBuffer(1).datatype = 'int16'; % 16 bit signed integers are the only supported datatype

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

startSample = (0:(na-1))*numRcvSamples + 1;
endSample = startSample + numRcvSamples - 1;
%%%%

% spw = 3.6765; % samples per wave, it isn't always exactly 4... check p107
% nspa = spw*(2*(Receive(1).endDepth - Receive(1).startDepth));
% nspa = 128 * ceil(nspa/128); % # samples per acquisition
% maxAcqLength_adjusted = nspa / spw / 2;
Resource.RcvBuffer(1).rowsPerFrame = numRcvSamples * na * 2 * numSubFrames;
maxAcqLength_adjusted = numRcvSamples / samplesPerWave / 2;

% Commenting below section because it doesn't work for the second set of
% TXs
% for lss = 1:length(startSample)
%     Receive(lss).startSample = startSample(lss);
%     Receive(lss).endSample = endSample(lss);    
% %     Receive(lss).decimSampleRate = samplesPerWave * Trans.frequency;
% %     Receive(lss).decimSampleRate = 62.5;
% 
% end

Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels; % Usually 1:1 to # of receive channels available in the system. Can change to 256 with the 2D probe and new connector plate.

Resource.RcvBuffer(1).numFrames = numSupFrames; % minimum # frames of RF data to acquire; RcvBuffer contains all the data needed for a whole frame, including multiple acquisition passes needed for reconstruction. Software can re-process RcvBuffer frames
% Resource.InterBuffer(1).pagesPerFrame = pair*na*numSubFrames; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numSamplesInBuffer = Resource.RcvBuffer(1).rowsPerFrame * Resource.RcvBuffer(1).colsPerFrame * Resource.RcvBuffer(1).numFrames
numGBInBuffer = numSamplesInBuffer ./ 1024^3 * 2 % # samples * (2 bytes per int16 sample) 
numSamplesPerBufferFrame = Resource.RcvBuffer(1).rowsPerFrame * Resource.RcvBuffer(1).colsPerFrame
numGBPerBufferFrame = numSamplesPerBufferFrame ./ 1024^3 * 2 % # samples * (2 bytes per int16 sample) 

if numGBPerBufferFrame > 2
    warning('Buffer size per frame is too large, exiting')
    return

end


%% Reconstruction
numRegions = 3;

Resource.ImageBuffer(1).numFrames = numSupFrames; % Define an ImageBuffer with a # of frames
Resource.InterBuffer(1).numFrames = numSupFrames; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resource.InterBuffer(1).numFrames = 1;

% Recon = struct('senscutoff', 0.6, ... % Threshold for which the reconstruction doesn't consider an element's contribution due to directivity of the element, for a certain pixel (whose echoes are at an angle to the element). Should be in radians.
%                'pdatanum', 1, ... % Which PData structure to use
%                'rcvBufFrame', -1, ... % Use the most recently transferred frame
%                'IntBufDest', [1, 1], ... % idk but it's for the IQ (complex) data
%                'ImgBufDest', [1, -1], ... % [buffer #, frame #] Auto-increment ImageBuffer for each reconstruction???? % something is [first/oldest frame, last/newest frame]
%                'RINums', [1:2*na]); % The ReconInfo structure #(s). Each Recon must have its own unique set of ReconInfo #s

sco = 0.6; %%%%
% sco = 0.4;
Recon = struct('senscutoff', sco, ... % Threshold for which the reconstruction doesn't consider an element's contribution due to directivity of the element, for a certain pixel (whose echoes are at an angle to the element). Should be in radians.
               'pdatanum', 1, ... % Which PData structure to use
               'rcvBufFrame', -1, ... % Use the most recently transferred frame
               'IntBufDest', [1, -1], ... % IQ (complex) data, Auto-increment for every frame
               'ImgBufDest', [1, -1], ... % [buffer #, frame #] Auto-increment ImageBuffer for each reconstruction???? % something is [first/oldest frame, last/newest frame]
               'RINums', [1:2*na]); % The ReconInfo structure #(s). Each Recon must have its own unique set of ReconInfo #s

% Recon = repmat(Recon, 1, numFrames);
% for nf = 1:numFrames
%     Recon(nf).IntBufDest = [1, nf];
%     Recon(nf).ImgBufDest = [1, nf];
% end

ReconInfo = repmat(struct('mode', 'accumIQ', ... % reconstruct, and replace intensity data in ImageBuffer and IQ data in InterBuffer (see Table 12.4 in Tutorial)
                   'txnum', 1, ...                 % TX structure to use
                   'rcvnum', 1, ...                % RX structure to use
                   'regionnum', ReconRegion), 1, 2*na);                % PData Region to process in

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for nf = 1
for n = 1:2*na % need to change this and above for more than 1 frame
    % - Set specific ReconInfo attributes.
    % ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
    ReconInfo(n).txnum = n;
    ReconInfo(n).rcvnum = n;
%     ReconInfo(n).pagenum = n; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ReconInfo(1).regionnum = 1; %1 for the whole volume, 5 for the slices

end

ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
ReconInfo(2*na).mode = 'accumIQ_replaceIntensity';

%% Process the Reconstructed data
% e.g., scaling and compression to make an image look good on the screen

pgainValue = 1.0; % Image processing gain
persValue = 0;
rejectLevel = 2;
% rejectLevel = 300;
compFac = 40;

% First image (xz)
Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';             % To not overwrite orig data while processing, system uses another ImageP buffer as output.
Process(2).Parameters = {'imgbufnum', 1, ...             % which ImageBuffer to process
                         'framenum', -1, ...              % -1 means use last frame in ImageBuffer
                         'pdatanum', 1, ...              % PData structure which was used in Reconstruction
                         'srcData', 'intensity3D', ...
                         'pgain', pgainValue, ...
                         'reject', rejectLevel, ...                % Make intensity values below this threshold appear as black (reduce low intensity noise)
                         'persistMethod', 'simple', ...   % simple: Add a fraction of the previous weighted average frames' invensity values to the current one. See manual for 'dynamic' option, which is good when there is a lot of motion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         'persistLevel', persValue, ...   % NewAvg = PL * PrevAvg + (1 - PL) * NewFrame; PL = persistLevel/100
                         'interpMethod', '4pt', ...
                         'grainRemoval', 'medium', ...    % low, medium, high. Remove pixels that differ significantly from their neighbors
                         'processMethod', 'none', ...     % see manual, reduces variation in line structures detected within the filter kernel???
                         'averageMethod', 'none', ...     % None or can do Running averages (2 or 3), can do things like spatial compounding...
                         'compressMethod', 'log', ...     % log or power (x^a fraction) compression
                         'compressFactor', compFac, ...        % Higher compressFactor means smaller powers for the power option (more compression), or a more rapid rise to the log curve (raise brightness of low intensity values). Not a real log bc intensities of 0 need to be mapped to 0, not -inf
                         'mappingMethod', 'full', ...     % Portion of the colormap to use. lowerHalf and upperHalf would be used to do the combined B-mode and Doppler imaging, for example.
                         'display', 1, ...                % 1: show processed image on screen, 0: don't but still tto useransfer the processed data to the DisplayData buffer
                         'displayWindow', 1};             % which displayWindow 

% Second image (xy)
Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';             % To not overwrite orig data while processing, system uses another ImageP buffer as output.
Process(3).Parameters = {'imgbufnum', 1, ...             % which ImageBuffer to process
                         'framenum', -1, ...              % -1 means use last frame in ImageBuffer
                         'pdatanum', 1, ...              % PData structure which was used in Reconstruction
                         'srcData', 'intensity3D', ...
                         'pgain', pgainValue, ...
                         'reject', rejectLevel, ...                % Make intensity values below this threshold appear as black (reduce low intensity noise)
                         'persistMethod', 'simple', ...   % simple: Add a fraction of the previous weighted average frames' invensity values to the current one. See manual for 'dynamic' option, which is good when there is a lot of motion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         'persistLevel', persValue, ...   % NewAvg = PL * PrevAvg + (1 - PL) * NewFrame; PL = persistLevel/100
                         'interpMethod', '4pt', ...
                         'grainRemoval', 'medium', ...    % low, medium, high. Remove pixels that differ significantly from their neighbors
                         'processMethod', 'none', ...     % see manual, reduces variation in line structures detected within the filter kernel???
                         'averageMethod', 'none', ...     % None or can do Running averages (2 or 3), can do things like spatial compounding...
                         'compressMethod', 'log', ...     % log or power (x^a fraction) compression
                         'compressFactor', compFac, ...        % Higher compressFactor means smaller powers for the power option (more compression), or a more rapid rise to the log curve (raise brightness of low intensity values). Not a real log bc intensities of 0 need to be mapped to 0, not -inf
                         'mappingMethod', 'full', ...     % Portion of the colormap to use. lowerHalf and upperHalf would be used to do the combined B-mode and Doppler imaging, for example.
                         'display', 1, ...                % 1: show processed image on screen, 0: don't but still tto useransfer the processed data to the DisplayData buffer
                         'displayWindow', 2};             % which displayWindow 

% Third image (yz)
Process(4).classname = 'Image';
Process(4).method = 'imageDisplay';             % To not overwrite orig data while processing, system uses another ImageP buffer as output.
Process(4).Parameters = {'imgbufnum', 1, ...             % which ImageBuffer to process
                         'framenum', -1, ...              % -1 means use last frame in ImageBuffer
                         'pdatanum', 1, ...              % PData structure which was used in Reconstruction
                         'srcData', 'intensity3D', ...
                         'pgain', pgainValue, ...
                         'reject', rejectLevel, ...                % Make intensity values below this threshold appear as black (reduce low intensity noise)
                         'persistMethod', 'simple', ...   % simple: Add a fraction of the previous weighted average frames' invensity values to the current one. See manual for 'dynamic' option, which is good when there is a lot of motion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         'persistLevel', persValue, ...   % NewAvg = PL * PrevAvg + (1 - PL) * NewFrame; PL = persistLevel/100
                         'interpMethod', '4pt', ...
                         'grainRemoval', 'medium', ...    % low, medium, high. Remove pixels that differ significantly from their neighbors
                         'processMethod', 'none', ...     % see manual, reduces variation in line structures detected within the filter kernel???
                         'averageMethod', 'none', ...     % None or can do Running averages (2 or 3), can do things like spatial compounding...
                         'compressMethod', 'log', ...     % log or power (x^a fraction) compression
                         'compressFactor', compFac, ...        % Higher compressFactor means smaller powers for the power option (more compression), or a more rapid rise to the log curve (raise brightness of low intensity values). Not a real log bc intensities of 0 need to be mapped to 0, not -inf
                         'mappingMethod', 'full', ...     % Portion of the colormap to use. lowerHalf and upperHalf would be used to do the combined B-mode and Doppler imaging, for example.
                         'display', 1, ...                % 1: show processed image on screen, 0: don't but still tto useransfer the processed data to the DisplayData buffer
                         'displayWindow', 3};             % which displayWindow 

Process(1).classname = 'External';
Process(1).method = 'saveRcvData'; % Function name
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

Process(1).Parameters = {'srcbuffer', 'receive', ...
                         'srcbufnum', 1, ... % # of buffer to process
                         'dstbuffer', 'none'};
%                          'srcframenum', -1, ... % last frame transferred
%%
makeParameterStructureSmall;
%% New Event structure

Resource.VDAS.dmaTimeout = 1000;
% Flow:
% 1. Transmit (TX)
% 2. Receive (Receive)
% 3. Reconstruction (Recon)
% 4. Processing (Process)
% 5. Control (SeqControl)


% SeqControl(1).argument = 150;  
scInd = 1; % sequence control index
SeqControl(scInd).command = 'timeToNextAcq'; % In us, allowed range is from 10 - 4190000
                                         % Very useful if you are switching
                                         % the TPC (voltage) across acqs,
                                         % since it takes 800 us - 8 ms to
                                         % switch
timePerAcq = 1 / fps_target / (na * 2) * 1e6; % frame rate converted to acq time step (us)
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

timePerSubframe = SeqControl(scInd).argument * na * 2; % us
timeGapBetweenLastSuperframeBurstAndNextSuperframeBurst = 1 / supFrameBurstRate * 1e6 - timePerSubframe * numSubFrames * numSupFrames;

scInd = scInd + 1;

SeqControl(scInd).command = 'returnToMatlab';

scInd = scInd + 1;

SeqControl(scInd).command = 'jump'; % jump to
SeqControl(scInd).argument = 1;     % first event

% superframe burst rate
scInd = scInd + 1;
SeqControl(scInd).command = 'timeToNextAcq'; % jump to

if timeGapBetweenLastSuperframeBurstAndNextSuperframeBurst < timePerAcqLimits(1)
    warning('Superframe burst delay time too short, setting to minimum of 10 us')
    SeqControl(scInd).argument = timePerAcqLimits(1); 
elseif timeGapBetweenLastSuperframeBurstAndNextSuperframeBurst > timePerAcqLimits(2)
    warning('Superframe burst delay time too long, setting to maximum of 4190000 us')
    SeqControl(scInd).argument = timePerAcqLimits(2);
else
    SeqControl(scInd).argument = timeGapBetweenLastSuperframeBurstAndNextSuperframeBurst;
end

% Transfer data to host, needed for hardware but not simulation, which writes directly to RcvBuffer
% NEED A UNIQUE SEQCONTROL FOR EACH TRANSFERTOHOST COMMAND!!!!!!!!!!!!!!!!!!

n = 0;
for nsupf = 1:numSupFrames
    for nsubf = 1:numSubFrames
    
        for a = 1:na % go through all the angles for each frame
            n = n + 1;
            Event(n).info = 'Transmit all columns and receive all rows';
            Event(n).tx = a.*2 - 1; % Use ath TX structure
            Event(n).rcv = (nsupf - 1) .* numSubFrames .* pair .* na + (nsubf-1).*2.*na + a.*2 - 1; % Use nth Receive structure % need to make this alternate between (1 and 2) * numframes or something
            Event(n).recon = 0; % 0 means no reconstruction
            Event(n).process = 0; % 0 means no processing
            Event(n).seqControl = 1;
    %         Event((nsubf-1).*2.*na + a.*2 - 1 + i).seqControl = 0;
            
            n = n + 1;
            Event(n).info = 'Transmit all rows and receive all columns';
            Event(n).tx = a.*2; 
            Event(n).rcv = (nsupf - 1) .* numSubFrames .* pair .* na + (nsubf-1).*pair.*na + a.*2; 
            Event(n).recon = 0; 
            Event(n).process = 0; 
            Event(n).seqControl = 1;  
    %         Event((nsubf-1).*2.*na + a.*2 + i).seqControl = 0;
        
        end
    end

    scInd = scInd + 1;
%     Event(n).seqControl = [Event(n).seqControl, scInd]; 
    SeqControl(scInd).command = 'transferToHost';

% test 11/6/24
%     scInd = scInd + 1;
    
%     SeqControl(scInd).command = 'waitForTransferComplete';
%     SeqControl(scInd).argument = scInd - 1;
% 
%     Event(n).seqControl = [1, scInd, scInd-1]; % transfer and wait for the transfer
    Event(n).seqControl = [1, scInd];
    
    n = n + 1;

    Event(n).info = ['Reconstruction'];
    Event(n).tx = 0; 
    Event(n).rcv = 0; 
    Event(n).recon = 1;  %%
    Event(n).process = 0; 
    Event(n).seqControl = 0; 

    n = n + 1;

    Event(n).info = ['Processing xz'];
    Event(n).tx = 0; 
    Event(n).rcv = 0; 
    Event(n).recon = 0; 
    Event(n).process = 2; 
    Event(n).seqControl = 0; 

    n = n + 1;
    
    Event(n).info = ['Processing xy'];
    Event(n).tx = 0; 
    Event(n).rcv = 0; 
    Event(n).recon = 0; 
    Event(n).process = 3; 
    Event(n).seqControl = 0; 

    n = n + 1;
    
    Event(n).info = ['Processing yz'];
    Event(n).tx = 0; 
    Event(n).rcv = 0; 
    Event(n).recon = 0; 
    Event(n).process = 4; 
    Event(n).seqControl = 0;
% 
%     n = n + 1;

end

% Event(n).seqControl = [4, scInd, scInd-1]; % transfer and wait for the transfer

% n = n + 1;
% 
% Event(n).info = 'Save data - ext proc func';
% Event(n).tx = 0; 
% Event(n).rcv = 0; 
% Event(n).recon = 0;
% Event(n).process = 1; 
% Event(n).seqControl = 0; 

% Test for time tag
% n = n+1;
% Event(n).info = 'ext function call';
% Event(n).tx = 0;         % no transmit
% Event(n).rcv = 0;        % no rcv
% Event(n).recon = 0;      % reconstruction
% Event(n).process = 2;    % external processing function
% Event(n).seqControl = 0; % none

n = n + 1;

Event(n).info = 'Jump';
Event(n).tx = 0; 
Event(n).rcv = 0; 
Event(n).recon = 0;
Event(n).process = 0; 
Event(n).seqControl = 3; 

%  n = n + 1;
% 
%     Event(n).info = 'Save data - ext proc func';
%     Event(n).tx = 0; 
%     Event(n).rcv = 0; 
%     Event(n).recon = 0;
%     Event(n).process = 1; 
%     Event(n).seqControl = 2; 


%% User specified UI Control Elements

%% User specified UI Control Elements
import vsv.seq.uicontrol.VsSliderControl;

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f', ...
                  'Callback', @SensCutoffCallback );


UI(2).Control = VsSliderControl('LocationCode','UserB4','Label','XY Slider',...
              'SliderMinMaxVal',[startDepth,endDepth,xyplane],...
              'SliderStep',[0.0125,0.1],'ValueFormat','%1.3f', ...
              'Callback', @ChangeXYCallback );
        
UI(3).Control = VsSliderControl('LocationCode','UserB5','Label','XZ Slider',...
              'SliderMinMaxVal',[min(P.yCoord),max(P.yCoord),PData(1).Region(3).Shape.oPAIntersect],...
              'SliderStep',[0.0125,0.1],'ValueFormat','%1.3f', ...
              'Callback', @ChangeXZCallback );
          
UI(4).Control = VsSliderControl('LocationCode','UserB6','Label','YZ Slider',...
              'SliderMinMaxVal',[min(P.xCoord),max(P.xCoord),PData(1).Region(2).Shape.oPAIntersect],...
              'SliderStep',[0.0125,0.1],'ValueFormat','%1.3f', ...
              'Callback', @ChangeYZCallback );


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

%% Save all the data/structures to a .mat file.
currentDir = cd; currentDir = regexp(currentDir, filesep, 'split');
filename = 'RC15gV_Allen_loop_live.mat';

save(fullfile(currentDir{1:find(contains(currentDir,"Vantage"),1)})+"\MatFiles\"+filename);


%% Run VSX automatically and make parameter structure for RF file naming

if runVSX
    disp("running VSX")
    VSX
end

%% Save post-acquisition parameters in a structure P

% save([savepath, 'params.mat'], 'angles', 'startDepth', 'startDepthMM', 'endDepth', 'endDepthMM', 'Event', 'fps_target', 'maxAcqLength_adjusted', 'maxAngle', 'Media', 'na', 'nf', 'Receive', 'Resource', 'SeqControl', 'TGC', 'Trans', 'TW', 'TX', 'wl', 'numElements', '-v7.3')
% save([savepath, 'params.mat'], 'P')

% makeParameterStructure;
% savefast([savepath, 'params.mat'], 'P')
% saveRcvData(RcvData{1})
% save([savepath, 'workspace.mat'])

%% **** Callback routines used by UIControls (UI) ****
%% Time tag callback test

function SensCutoffCallback(~,~,UIValue)
    %SensCutoff - Sensitivity cutoff change
    ReconL = evalin('base', 'Recon');
    for i = 1:size(ReconL,2)
        ReconL(i).senscutoff = UIValue;
    end
    assignin('base','Recon',ReconL);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Recon'};
    assignin('base','Control', Control);
end


function ChangeXYCallback(~,~,UIValue)    
    % Import relevant structs
    tic;
    Resource = evalin('base','Resource');
    PData = evalin('base','PData');
    P = evalin('base','P');
    
    % Adjust PData Region
    [~,I] = min(abs(UIValue-P.zCoord));
    RegionIndices = reshape(PData.Region(1).PixelsLA,PData.Size);
    PData.Region(4).PixelsLA = reshape(RegionIndices(:,:,I),[],1,1);
    PData(1).Region(4).Shape.oPAIntersect = P.zCoord(I);
    
    % Adjust Display Window
    Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),-PData(1).Origin(2),PData.Region(4).Shape.oPAIntersect];
    
    % Update Overlapped region
    PData.Region(5).PixelsLA = unique([PData.Region(2).PixelsLA; PData.Region(3).PixelsLA; PData.Region(4).PixelsLA]);
    PData.Region(5).numPixels = length(PData.Region(5).PixelsLA);

    % Update relevant controls and commands
    assignin('base','Resource',Resource);
    assignin('base','PData',PData);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','DisplayWindow','Recon'};
    assignin('base','Control', Control);
    toc
end

function ChangeXZCallback(~,~,UIValue)    
    % Import relevant structs
    tic;
    Resource = evalin('base','Resource');
    PData = evalin('base','PData');
    P = evalin('base','P');
    
    % Adjust PData Region
    [~,I] = min(abs(UIValue-P.yCoord));
    RegionIndices = reshape(PData.Region(1).PixelsLA,PData.Size);
    PData.Region(3).PixelsLA = reshape(RegionIndices(I,:,:),[],1,1);
    PData(1).Region(3).Shape.oPAIntersect = P.yCoord(I);
    
    % Adjust Display Window
    Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),PData(1).Region(3).Shape.oPAIntersect,0.0];
    
    % Update Overlapped region
    PData.Region(5).PixelsLA = unique([PData.Region(2).PixelsLA; PData.Region(3).PixelsLA; PData.Region(4).PixelsLA]);
    PData.Region(5).numPixels = length(PData.Region(5).PixelsLA);

    % Update relevant controls and commands
    assignin('base','Resource',Resource);
    assignin('base','PData',PData);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','DisplayWindow','Recon'};
    assignin('base','Control', Control);

    toc
end

function ChangeYZCallback(~,~,UIValue)    
    % Import relevant structs
    tic;
    Resource = evalin('base','Resource');
    PData = evalin('base','PData');
    P = evalin('base','P');
    
    % Adjust PData Region
    [~,I] = min(abs(UIValue-P.xCoord));
    RegionIndices = reshape(PData.Region(1).PixelsLA,PData.Size);
    PData.Region(2).PixelsLA = reshape(RegionIndices(:,I,:),[],1,1);
    PData(1).Region(2).Shape.oPAIntersect = P.xCoord(I);
    
    % Adjust Display Window
    Resource.DisplayWindow(2).ReferencePt = [PData(1).Region(2).Shape.oPAIntersect,-PData(1).Origin(2),0];
    
    % Update Overlapped region
    PData.Region(5).PixelsLA = unique([PData.Region(2).PixelsLA; PData.Region(3).PixelsLA; PData.Region(4).PixelsLA]);
    PData.Region(5).numPixels = length(PData.Region(5).PixelsLA);

    % Update relevant controls and commands
    assignin('base','Resource',Resource);
    assignin('base','PData',PData);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','DisplayWindow','Recon'};
    assignin('base','Control', Control);
    toc
end