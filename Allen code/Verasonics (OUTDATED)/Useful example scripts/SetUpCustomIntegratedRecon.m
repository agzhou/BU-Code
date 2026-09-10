% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpGE9LDFlashAngles.m - Example of curved array flash imaging
%                                      with steering angle transmits
% Description:
%   Sequence programming file for GE9LD curved array, using flash transmits
%   with multiple steering angles. 192 transmit and receive channels
%   are active for each acquisition. Processing is asynchronous with
%   respect to acquisition.
%
% Testing: Tested with software release 3.2
%
% Last update:
%   08-16-2016 Modified for SW 3.2

clearvars
% close all

%% Initial Parameters and Structure Definition
P.startDepth = 0;  % P.startDepth and P.endDepth are in wavelength
P.endDepth = 450;

na = 15;
if (na > 1)
    dtheta = (48*pi/180)/(na-1);
    startAngle = -24*pi/180;
else
    dtheta = 0;
    startAngle=0;
end % set dtheta to range over +/- 12 degrees.

TXangle = linspace(startAngle,-startAngle,na);

% Specify system parameters.
Resource.Parameters.numTransmit = 256;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
Resource.Parameters.waitForProcessing = 1;

% Specify Trans structure array.
Trans.name = 'GE9LD';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 35;  % set maximum high voltage limit for pulser supply.

lambda = Resource.Parameters.speedOfSound/(Trans.frequency*1e6);

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.

% Specify Media object. 'pt1.m' script defines array of point targets.
% pt1;
USTBPointTargets;
Media.attenuation = -0.5;
% Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = na*4096;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
% Resource.RcvBuffer(1).numFrames = 20;     % 20 frames for rcvDataLoop buffer.
Resource.RcvBuffer(1).numFrames = 1;     % 20 frames for rcvDataLoop buffer.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.

Resource.ImageBuffer(1).numFrames = 1;
Resource.ImageBuffer(2).numFrames = 1;

% First Display window for Verasonics Recon
Resource.DisplayWindow(1).Title = 'Vsonics Reconstruction';
Resource.DisplayWindow(1).pdelta = 0.45;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [150,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Matlab';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Second Display window for custom recon
Resource.DisplayWindow(2).Title = 'Custom Reconstruction';
Resource.DisplayWindow(2).pdelta = 0.45;
Resource.DisplayWindow(2).Position = [150+DwWidth*1.5,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(2).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(2).Type = 'Matlab';
Resource.DisplayWindow(2).numFrames = 20;
Resource.DisplayWindow(2).AxesUnits = 'mm';
Resource.DisplayWindow(2).Colormap = gray(256);

%% Specify Transmit Structures
% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

[Wvfm2Wy, peak, numsamples, rc, ~] = computeTWWaveform(TW(1));

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
% - Set event specific TX attributes.
for n = 1:na   % na transmit events
    TX(n).Steer = [(startAngle+(n-1)*dtheta),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [153,308,410,520,605,665,705,760];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

TPC(1).hv = 20;

%% Specify Receive structure arrays.
% - We need na Receives for every frame.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, na*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(na*(i-1)+1).callMediaFunc = 1;
    for j = 1:na
        Receive(na*(i-1)+j).framenum = i;
        Receive(na*(i-1)+j).acqNum = j;      % two acquisitions per frame
    end
end


%% Verasonics Reconstruction
% - We need one Recon structure which will be used for each frame.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [2,-1], ...
               'RINums',1:na);

% Define ReconInfo structures.
% We need na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 4, ...  % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, na);
% - Set specific ReconInfo attributes.
if na>1
    ReconInfo(1).mode = 'replaceIQ';
    for j = 1:na  % For each row in the column
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
    end
    ReconInfo(na).mode = 'accumIQ_replaceIntensity';  % accumulate and detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end

%% Custom Reconstruction
% - We need a process structure for the external beamforming class
Process(1).classname = 'External';
Process(1).method = 'externalBeamforming';
Process(1).Parameters = {'srcbuffer','receive',...  % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,...
                         'dstbuffer','image',...
                         'dstbufnum',1,...
                         'dstframenum',-1};

% Define parameter structure for beamforming class
xCoord = (PData(1).Origin(1) + (0:PData(1).Size(2)-1)*PData(1).PDelta(1))*lambda;
zCoord = (PData(1).Origin(3) + (0:PData(1).Size(1)-1)*PData(1).PDelta(3))*lambda;

% for UTA 408-GE, use code [1 7 1 0] in computeUTA to find UTA.TransConnector
chkUTA = computeUTA([1 7 1 0]);

if strcmp(Receive(1).sampleMode,'custom')
    error('No handling of condition for custom Receive sampling. Refer to VsUpdate line 712 to implement');
else
    fs = 4*Trans.frequency;
    samplesPerWave = 4;
end

% if statement included to match verasonics automatic extension to
% multiples of 128 samples
nSmpls = 2*(maxAcqLength - P.startDepth) * samplesPerWave;
if abs(round(nSmpls/128) - nSmpls/128) < .01
    numRcvSamples = 128*round(nSmpls/128);
else
    numRcvSamples = 128*ceil(nSmpls/128);
end

startSample = (0:(na-1))*numRcvSamples + 1;
endSample = startSample + numRcvSamples - 1;

param = struct('fs',fs*1e6,...
        'pitch', Trans.spacingMm*1e-3,...
        'fc', Trans.frequency*1e6,...
        'c', Resource.Parameters.speedOfSound,...
        'fnumber', 0.6,...
        't0',(Receive(1).startDepth + peak + 2*Trans.lensCorrection)/(Trans.frequency*1e6),...
        'TXangle',TXangle,...
        'ElemPos',Trans.ElementPos(:,1).'*lambda,...
        'xCoord',xCoord,...
        'zCoord',zCoord,...
        'numEl',int32(Trans.numelements),...
        'szRF',Resource.RcvBuffer.rowsPerFrame,...
        'szRFframe',numRcvSamples - 1,...
        'szX',length(xCoord),...
        'szZ',length(zCoord),...
        'na',na,...
        'nc',Resource.Parameters.numRcvChannels,...
        'ConnMap',chkUTA.TransConnector(Trans.ConnectorES),...
        'startSample',int32(startSample),...
        'endSample',int32(endSample),...
        'initFlag',false);
    
beamformer = vsv.reconraw.DASBMode(param);

%% Specify Process structure array.
pers = 20;
cmpFactor = 40;
Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...     % pgain is image processing gain
                         'reject',2,...
                         'grainRemoval','none',...
                         'persistMethod','none',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',cmpFactor,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

pers = 20;
cmpFactor = 40;
Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...     % pgain is image processing gain
                         'reject',2,...
                         'grainRemoval','none',...
                         'persistMethod','none',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',cmpFactor,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',2};
                     
                     
%% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 290*2;  % 290 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 20000;  % 10 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

%% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:na                      % Acquire frame
        Event(n).info = 'Full aperture.';
        Event(n).tx = j;
        Event(n).rcv = na*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'Verasonics Recon';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n+1;
    
    Event(n).info = 'Custom Recon';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;
    
    Event(n).info = 'VSonics process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;
    
    Event(n).info = 'Custom Recon process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    if (floor(i/5) == i/5)&&(i ~= Resource.RcvBuffer(1).numFrames)  % Exit to Matlab every 5th frame
        Event(n).seqControl = 4;
    else
        Event(n).seqControl = 0;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;


% External function
EF(1).Function = vsv.seq.function.ExFunctionDef('externalBeamforming', @externalBeamforming);


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
filename = 'CustomIntegratedRecon.mat';
currentDir = cd;currentDir = regexp(currentDir, filesep, 'split');
save(fullfile(currentDir{1:find(contains(currentDir,"Vantage"),1)})+"\MatFiles\"+filename);
VSX



%% Further testing
nPoints = beamformer.szZ*beamformer.szX;
tic; cRF = computeCRF(RcvData{1}(:,:,1),beamformer); toc
tic; idxtMTX = beamformer.computeDASFullKSpace(cRF); toc
idxtMTX = reshape(idxtMTX,[beamformer.numEl,beamformer.na,nPoints]);
imgCustomFull = abs(reshape(sum(sum(idxtMTX,1),2),[beamformer.szZ,beamformer.szX]));

imgCustomCRF = ImgDataP{1}(:,:,1);
tic; imgCustomCRF2 = beamformer.computeDAScrfBMode(cRF); toc
imgCustomCRF2 = abs(imgCustomCRF2);
imgVsonics = ImgDataP{2}(:,:,1);

% imgCustomFull = abs(reshape(squeeze(sum(idxtMTX)),[beamformer.szZ,beamformer.szX]));
%% Post Plotting comparison
figure
plotLogScaleImage(beamformer.xCoord*1e3,beamformer.zCoord*1e3,imgCustomCRF)
title('Custom Recon')
axis image

figure
plotLogScaleImage(beamformer.xCoord*1e3,beamformer.zCoord*1e3,imgVsonics)
title('Verasonics Recon')
axis image

figure
plotLogScaleImage(beamformer.xCoord*1e3,beamformer.zCoord*1e3,imgCustomCRF2)
title('External Custom Recon')
axis image

figure
plotLogScaleImage(beamformer.xCoord*1e3,beamformer.zCoord*1e3,imgCustomFull)
title('External Full Recon')
axis image
%% **** Callback routines used by External Function. ****

function ImageData = externalBeamforming(ReceiveData)
%     persistent customRecon
%     persistent VsonicsRecon
    
%     if isempty(customRecon)||~ishandle(customRecon)
%         customRecon = figure('name','Custom','NumberTitle','off');
%     end
%     
%     if isempty(VsonicsRecon)||~ishandle(VsonicsRecon)
%         VsonicsRecon = figure('name','VsonicsRecon','NumberTitle','off');
%     end

    % Get external beamformer from the base workspace
    beamformer = vsv.seq.getBaseComp('beamformer');
    
    ImgData = vsv.seq.getBaseComp('ImgData');
    
    % if any initial parameters need to be initialized
    if (~beamformer.initFlag)
        beamformer.initFlag = true;
    end
    
    cRF = zeros(size(ReceiveData));
    for i = 1:beamformer.na
        cRF(beamformer.startSample(i):beamformer.endSample(i),:) = hilbert(ReceiveData(beamformer.startSample(i):beamformer.endSample(i),:));
    end
    cRF([beamformer.startSample],:,:) = 0;
    
    % Calculate the beamformed image
    ImageData = beamformer.computeDAScrfBMode(cRF);
    
    
%     set(0,'CurrentFigure',customRecon)
%     colormap gray
%     imagesc(beamformer.xCoord*1e3,beamformer.zCoord*1e3,10*log10(abs(ImageData)))
%     axis image
    
    ImageData = abs(ImageData);
    
%     set(0,'CurrentFigure',VsonicsRecon)
%     colormap gray
%     imagesc(10*log10(abs(ImgData{1}(:,:,1,1))))
end