%
% File name: SetUpVermon1024_8MHz_FlashAngles_light_superframes.m
%            Example of plane wave imaging with multiple steered angle transmits in a light synthetic
%            aperture sequence.  This scripts implements superframes to allow for fast DMA transfers
%            as well as realtime reconstruction to monitor the acquisitions.
%
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage NXT Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% Copyright © 2013-2025 Verasonics, Inc.

clear all

fsFileName = 'testFileName';

if isunix()
    fsDriveNames = {'/media/verasonics/WD0','/media/verasonics/WD1','/media/verasonics/WD2','/media/verasonics/WD3','/media/verasonics/WD4'};
    %fsDriveNames = {'/media/verasonics/WD1','/media/verasonics/WD2','/media/verasonics/WD3','/media/verasonics/WD4'};
    %fsDriveNames = {'/media/verasonics/WD0','/media/verasonics/WD1'};
    %fsDriveNames = {'/media/verasonics/WD0'};
elseif ispc()
    fsDriveNames = {'D:\','E:\','F:\','G:\'};
end

Resource.VDAS.acqTimeTagEnable = true; %use time tags to validate timing of RF acqusitions

P.startDepth = 0;   % start of Acquisition depth range in wavelengths
P.endDepth = 128; %
P.viewAngle = 20 * pi/180; % angle between z-axis and surface of cone.

%  Setup 5 plane wave angles for spatial compounding
%NOTE: length of AZ must be the same as length of EL
AZ = [0 -5 5 0 0]*pi/180; % azimuthal tilt
EL = [0 0 0 -5 5]*pi/180; % elevational tilt
na = length(AZ);


% "Light" acquisition sequence for 2D matrix with 4:1 mux, just receiving
% on neighboring apertures + same aperture
TXapertures = [1 1 2 2 2 3 3 3 4 4];
numTxap = length(TXapertures);
Rcvapertures = [1 2 1 2 3 2 3 4 3 4];
numRcvap = length(Rcvapertures);
P.numAcqPerFrame = numTxap * na;

acqPeriod_usec = 220;  % 220us is shortest possible ttna for the imaging depth

% For picking the number of subframes: pick a large enough number so that
% the rcvdata frame size is >100MB but too large of a value will lower
% recon frame rate.

P.nSubFrames = 24;     % reconstruct 1 frame out of nSubFrames,
demodFrequency = 15.62; % 8 samples/wavelength oversample
Resource.RcvBuffer(1).rowsPerFrame = P.nSubFrames*P.numAcqPerFrame*2816;  % <-- 2816 can be found by running VSX and looking at the calculated value for Receive(1).endSample
% half the data
%demodFrequency = demodFrequency/2;
%Resource.RcvBuffer(1).rowsPerFrame = Resource.RcvBuffer(1).rowsPerFrame/2;

reconTime = 240e-3; % measured time to perform reconstruction on 1 subframe % Determined using vsProfile
P.numFrames = 8;  % min 2 for ping-pong, increment by 2, increase number so that frames aren't lost if the writes need to be buffered  2 x # of drives seems to work well

Resource.VDAS.watchdogTimeout = 0;  % Disable the watchdog timer
Resource.VDAS.dmaTimeout = 4000; % 4 seconds
acqDuration_s = 60*5; % 5 minutes
Resource.VDAS.acqTimeTagEnable = true; %use time tags to validate timing of RF acqusitions

RcvProfile.AntiAliasCutoff = 20;

ReconRegion = 4; %Set ReconRegion=4 to reconstruct only 2 orthogonal planes, or ReconRegion =1 to reconstruct the whole volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specify system parameters.
Resource.Parameters.numTransmit = 256;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;   % number of receive channels.
Resource.Parameters.verbose = 2; % 1 enables and 0 disables
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.fakeScanhead = 0; % allow system HW operation with nothing connected
Resource.Parameters.Connector = 1;
Resource.Parameters.simulateMode = 0;

%% Specify Trans structure array.
Trans.name = 'Matrix1024-8'; %or 'Matrix1024-3';
Trans.units = 'mm';
Trans = computeTrans(Trans);
Trans.HVMux = computeUTAMux1024; % add the HV Mux programming for the UTA 1024-MUX

Papod = [ ones(1,256) zeros(1,256) zeros(1,256) zeros(1,256); ...
    zeros(1,256)  ones(1,256) zeros(1,256) zeros(1,256); ...
    zeros(1,256) zeros(1,256)  ones(1,256) zeros(1,256); ...
    zeros(1,256) zeros(1,256) zeros(1,256)  ones(1,256)];

for i = 1:size(Papod, 1)
    Paper(i) = computeMuxAperture(Papod(i, :), Trans);
end


%Intermediate Variables
waveLength = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
if strcmp(Trans.units,'mm')
    Trans.ElementPosMm = Trans.ElementPos;
    Trans.ElementPosWL = Trans.ElementPos./waveLength;
else
    Trans.ElementPosMm = Trans.ElementPos.*waveLength;
    Trans.ElementPosWL = Trans.ElementPos;
end

%% PData
PData(1).PDelta = [Trans.spacing, Trans.spacing, 1];
if ~exist('extent','var'), extent = max(max(Trans.ElementPosWL(:,1)),max(Trans.ElementPosWL(:,2))); end
zApex = -extent/tan(P.viewAngle);

PData(1).Size(1) = ceil(2.0*(P.endDepth-zApex)*tan(P.viewAngle)/PData(1).PDelta(2));  if mod(PData(1).Size(1),2)==0, PData(1).Size(1) =  PData(1).Size(1)+1; end
PData(1).Size(2) = PData(1).Size(1);
PData(1).Size(3) = ceil((P.endDepth)/PData(1).PDelta(3));
PData(1).Origin = [-((PData(1).Size(2)-1)/2)*PData(1).PDelta(1), ((PData(1).Size(1)-1)/2)*PData(1).PDelta(1), 0];

PData(1).Region(1) = struct('Shape',struct('Name','Pyramid','Position',[0,0,zApex],'angle',P.viewAngle,'z1',P.startDepth,'z2',P.endDepth));
PData(1).Region(2) = struct('Shape',struct('Name','Slice','Orientation','yz','oPAIntersect',PData.Origin(1)+PData.PDelta(2)*(PData.Size(2)-1)/2,'andWithPrev',1));
PData(1).Region(3) = struct('Shape',struct('Name','Slice','Orientation','xz','oPAIntersect',PData.Origin(2)-PData.PDelta(1)*(PData.Size(1)-1)/2,'andWithPrev',1));
[PData(1).Region] = computeRegions(PData(1));

PData(1).Region(4).PixelsLA = unique([PData.Region(2).PixelsLA; PData.Region(3).PixelsLA]);
PData(1).Region(4).Shape.Name = 'Custom';
PData(1).Region(4).numPixels = length(PData.Region(4).PixelsLA);

%% Specify Media.  Use point targets in middle of PData.
Media.MP(1,:) = [0,0,80,1.0];
Media.attenuation = -0.5;
Media.numPoints = size(Media.MP,1);

%% Resources
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numFrames;

Resource.ImageBuffer(1).numFrames = P.numFrames;
Resource.InterBuffer(1).numFrames = 1;

Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).Title = '3D light - FlashAngles Image - XZ plane';
Resource.DisplayWindow(1).pdelta = 0.25;
Resource.DisplayWindow(1).Position = [0,580, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).Orientation = 'xz';
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0.0,0.0];
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).numFrames = 20;

Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).Title = '3D light - FlashAngles Image - YZ plane';
Resource.DisplayWindow(2).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(2).Position = [430,580, ...
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
Resource.DisplayWindow(2).Orientation = 'yz';
Resource.DisplayWindow(2).ReferencePt = [0,-PData(1).Origin(2),0];
Resource.DisplayWindow(2).Colormap = gray(256);
Resource.DisplayWindow(2).AxesUnits = 'mm';
Resource.DisplayWindow(2).numFrames = 20;

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,.67,2,1];  % A, B, C, D

TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', zApex, ...
    'Steer', [0.0,0.0], ...
    'Apod', zeros(1,1024), ...
    'Delay', zeros(1,1024),...
    'peakCutOff', 2,...
    'peakBLMax', 20,...
    'aperture',0), 1,na*numTxap);


for i = 1:na
    kk = numTxap*(i-1);
    for j=1:numTxap
        TX(j+kk).Steer= [AZ(i), EL(i)];
        TX(j+kk).Apod = squeeze(Papod(TXapertures(j),:));
        TX(j+kk).aperture = Paper(TXapertures(j));
        TX(j+kk).Delay = computeTXDelays(TX(j+kk));
    end
end

%Specify TGC Waveform structure.
TGC(1).CntrlPts = [0 298 416 595 773 892 1012 1023];
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC);

temp = (P.endDepth-zApex)*tan(P.viewAngle)+ max([max(Trans.ElementPosWL(:,1)), max(Trans.ElementPosWL(:,2))]);
maxAcqLength = sqrt(P.endDepth^2 + temp^2);

%Receive
Receive = repmat(struct(...
    'Apod', zeros(1,Trans.numelements), ...
    'startDepth', 0, ...
    'endDepth', 128/(4*2)*ceil(maxAcqLength/(128/(4*2))), ...
    'TGC', 1, ...
    'mode', 0, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'aperture',0, ...
    'demodFrequency', demodFrequency,...
    'callMediaFunc', 0), 1,  P.numFrames*numRcvap*na*P.nSubFrames);

rcvI = 1;
for frameI = 1:P.numFrames
    acqIndex = 1;
    for subFrameI = 1:P.nSubFrames
        for angleI = 1:na
            for w = 1:numRcvap
                Receive(rcvI).acqNum = acqIndex;
                Receive(rcvI).framenum = frameI;
                Receive(rcvI).aperture = Paper(Rcvapertures(w));
                Receive(rcvI).Apod = squeeze(Papod(Rcvapertures(w),:));
                acqIndex = acqIndex + 1;
                rcvI = rcvI + 1;
            end
        end
    end
end

%%
senscutoff = 0.6;
% Specify Recon structure arrays.

Recon = repmat(struct('senscutoff', senscutoff, ...
    'pdatanum', 1, ...
    'rcvBufFrame', -1, ...
    'IntBufDest', [1,1],...
    'ImgBufDest', [1,-1], ...
    'RINums', 1:na*numTxap), 1, 1);


% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...
    'txnum', 1, ...
    'rcvnum', 1, ...
    'scaleFactor', 1, ...
    'regionnum', ReconRegion), 1, na*numTxap);

% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';

for w=1:na
    k= numTxap*(w-1);
    for i = 1:numTxap
        ReconInfo(i+k).txnum = i+k;
        ReconInfo(i+k).rcvnum = i+k;
    end
end

ReconInfo(i+k).Post = 'IQ2IntensityImageBuf';


%%
pgain = 88;
% Specify Process structure arrays.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...
    'framenum',-1,...
    'pdatanum',1,...
    'srcData','intensity3D',...
    'pgain', pgain,...
    'persistMethod','none',...
    'persistLevel',30,...
    'interpMethod','4pt',...
    'compressMethod','log',...
    'compressFactor',55,...
    'mappingMethod','full',...
    'display',1,...
    'displayWindow',1};

Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',1,...
    'framenum',-1,...
    'pdatanum',1,...
    'srcData','intensity3D',...
    'pgain', pgain,...
    'persistMethod','none',...
    'persistLevel',30,...
    'interpMethod','4pt',...
    'compressMethod','log',...
    'compressFactor',55,...
    'mappingMethod','full',...
    'display',1,...
    'displayWindow',2};

Process(3).classname = 'External';
Process(3).method = 'fastsave';
Process(3).Parameters = { ...
    'srcbuffer','receive',... % name of buffer to export.
    'srcbufnum',    1,    ... % buffer number of buffer to export
    'srcframenum',  -1,   ... % process the last frame in a RcvBuffer
    'dstbuffer','none'
    };

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = acqPeriod_usec;
SeqControl(3).command = 'timeToNextAcq';  % time between superframes
SeqControl(3).argument = acqPeriod_usec;
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'triggerOut'; %use the trigger out to validate the actual HW frame rate of each superframe
SeqControl(5).argument = 2;
SeqControl(5).condition = {'SignalLevel', 5, 'Polarity', 'ActiveHigh'};

frameRateFactor = 1; % Factor for converting sequenceRate to frameRate.
nsc = length(SeqControl)+1; % Count of SeqControl objects.

n = 1; % Event structure arrays.
saveEvents = [];
for j = 1:Resource.RcvBuffer(1).numFrames
    for jj = 1:P.nSubFrames
        v = P.nSubFrames*numTxap*na*(j-1) + numTxap*na*(jj-1);
        for w = 1:na
            k = numTxap*(w-1);
            for i = 1:numTxap
                Event(n).info = sprintf('TX/RCV SuperFrame:%i SubFrame:%i Angle:%i Aperture:%i', j, jj, w, i);
                Event(n).tx = k+i;     % Loop over the TX structure.
                Event(n).rcv = i+k+v;  % Rcv structure of frame.
                Event(n).recon = 0;    % no reconstruction.
                Event(n).process = 0;  % no processing
                if (jj == 1)&&(w == 1)&&(i == 1)
                    Event(n).seqControl = [2,5]; % trigger out on the first transmit of the first subframe
                else
                    Event(n).seqControl = 2; %just TTNA between all other acqusitions
                end
                n = n+1;
            end
        end
    end
    Event(n-1).seqControl = 3;  %TTNA between frames

    Event(n).info = 'Transfer To Host';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [nsc, nsc + 1]; % set wait time and transfer data
    SeqControl(nsc).command = 'transferToHost';
    SeqControl(nsc+1).command = 'waitForTransferComplete';  % required or else an incomplete frame could be saved
    SeqControl(nsc+1).argument = nsc;
    nsc = nsc + 2;
    n = n+1;

    saveEvents = [saveEvents n];
    Event(n).info = 'save all subframes';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'markTransferProcessed';  %only required if not using live reconstruction
    SeqControl(nsc).argument = nsc - 2;
    nsc = nsc + 1;

    n = n+1;

    %if mod(j,2) %use this to reconstruct every other frame
    if true  %set to false if live monitoring is not needed
        Event(n).info = ['3D Reconstruct Frame ' num2str(j)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 1;    % Reconstruction of the slices in real time or the whole volume depending on ReconRegion
        Event(n).process = 0;
        Event(n).seqControl = 0;
        n = n+1;

        Event(n).info = ['Process XZ - Frame ' num2str(j)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 1; % Display of XZ plane in real time
        Event(n).seqControl = 0;
        n = n+1;

        Event(n).info = ['Process YZ - Frame ' num2str(j)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 2; % Display of YZ plane in real time
        Event(n).seqControl = 0;
        n = n+1;
    end
end

Event(n).info = 'Return to Matlab';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 4;
n=n+1;

Event(n).info = 'Jump';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;


%% User specified UI Control Elements
import vsv.seq.uicontrol.VsToggleButtonControl

UI(1).Control = VsToggleButtonControl('LocationCode','UserB1',...
    'Label','Start Fast Save',...
    'Callback',@saveAcquisition);

%% Calculations to report:
totalSpace = 0;
for a = 1:length(fsDriveNames)
    fileObj = java.io.File(fsDriveNames{a});
    totalSpace = totalSpace + double(fileObj.getTotalSpace())/1e9;
end
availableSpace_GB = totalSpace; %available space in GB = 2TB

superFrameSize_GB = Resource.RcvBuffer.rowsPerFrame*Resource.RcvBuffer.colsPerFrame*2/1e9;
%assert(superFrameSize_GB<2);  %2GB DMA LIMIT (no longer need to enforce for 2.2.0 and later
volumeRate = 1/(((na * numTxap - 1) * SeqControl(2).argument + SeqControl(3).argument ) * 1e-6);
superFrameRate = volumeRate/P.nSubFrames;
avgDataFlow_GBpS = superFrameSize_GB*superFrameRate;  %< this should be less than 7GB/s, the maxmimum dma rate ???
maxNumFiles = floor(availableSpace_GB/superFrameSize_GB);
maxAcqDuration_s = maxNumFiles/superFrameRate;
maxAcqDuration_m = maxAcqDuration_s/60;
requiredSaveSpeed = superFrameSize_GB*superFrameRate;

disp("=====================================================")
fprintf(" Available Storage space: %g GB\n",availableSpace_GB);
fprintf(" SuperFrame Size: %.2g GB\n", superFrameSize_GB)
fprintf(" Acquisition Rate: %.2g vol/sec\n", volumeRate)
fprintf(" Super-Volume Rate: %.2g superVol/sec\n", superFrameRate)
fprintf(" Average Data Throughput: %.4g GB/sec\n", avgDataFlow_GBpS)
fprintf(" Max # of files: %g\n", maxNumFiles)
fprintf(" Max Acq. Duration: %i sec == %i minutes\n", round(maxAcqDuration_s), round(maxAcqDuration_m))
disp("=====================================================")

% Save all the structures to a .mat file.
filename = 'MatFiles/Vermon1024_8MHz-FlashAngles-Light_superframe.mat';
save superFrameFrameRate superFrameRate
save(filename);

% **** Callback routines used by UI Controls. ****
function saveAcquisition(~,~,UIValue)
saveEvents = evalin('base', 'saveEvents');
Event = evalin('base', 'Event');
for i = 1:length(saveEvents)
    Event(saveEvents(i)).process = UIValue*3;  %for external
end
assignin('base','Event',Event);
Control = evalin('base','Control');
if isempty(Control(1).Command), n=1; else, n=length(Control)+1; end
Control(n).Command = 'update&Run';
Control(n).Parameters = {'Event'};
assignin('base','Control', Control);
end
