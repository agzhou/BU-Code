%
% File name: SetUpL22_14vXFlashAngles_fastSave.m - Example of plane wave imaging with multiple steered angle transmits
%
% Description:
%   Sequence programming file for L22-14vX Linear array, using a plane wave
%   transmits with multiple steering angles. All 128 transmit and receive
%   channels are active for each acquisition. This script uses 4X sampling with a 17.8571 MHz
%   processing center frequency. Processing is asynchronous with respect to acquisition.
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage NXT Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% Copyright © 2013-2025 Verasonics, Inc.

clear  all

P.nSubFrames = 50;
Resource.VDAS.acqTimeTagEnable = true; %use time tags to validate timing of RF acqusitions
fsFileName = 'testFileName';
if isunix()
    fsDriveNames = {'/media/verasonics/WD1','/media/verasonics/WD2','/media/verasonics/WD3','/media/verasonics/WD4'};
elseif ispc()
    fsDriveNames = {'E:\','F:\','G:\','H:\'};
end

P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 256;   % This should preferrably be a multiple of 128 samples.

na = 7;      % Set na = number of angles.
P.numAcqPerFrame = na;
if (na > 1)
    dtheta = (36*pi/180)/(na-1);
    P.startAngle = -36*pi/180/2;
else
    dtheta = 0;
    P.startAngle=0;
end % set dtheta to range over +/- 18 degrees.

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;

% Specify Trans structure array.
Trans.name = 'L22-14vX';
Trans.frequency = 18;
Trans = computeTrans(Trans);
%Trans.id = 0;

% Specify PData structure array.
PData(1).PDelta = [0.4, 0, 0.25];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).rowsPerFrame = P.nSubFrames*na*2432; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 30;    % 30 frames stored in RcvBuffer.
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L22-14vXFlashAngles';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [18, 0.67, 3, 1];
TW.clampHoldTimeus = 0.2;

RcvProfile.AntiAliasCutoff = 40;
RcvProfile.LnaGain         = 21;
RcvProfile.PgaGain         = 27;

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    P.startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
else
    P.startAngle = -fix(na/2)*dtheta;
end
for n = 1:na   % na transmit events
    TX(n).Steer = [(P.startAngle+(n-1)*dtheta),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [330 560 780 1010 1023 1023 1023 1023]; % [0,511,716,920,1023,1023,1023,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify RcvFilter structure arrays.
% Receive filter with a center frequency and demodFrequency of 17.8571 MHz and a 67%
% bandwidth
RcvFilter.sampleMode = 'NS200BW';
RcvFilter.demodFrequency = 17.8571;
RcvFilter.Bandwidth = [1-0.67/2, 1+0.67/2] * RcvFilter.demodFrequency; %With 0.67 as 67% bandwidth
RcvFilter = computeRcvFilter(RcvFilter);

% Specify Receive structure arrays.
% - We need na Receives for every frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'demodFrequency', 17.8571, ...
                        'rcvFilter', 1, ...
                        'callMediaFunc', 0), 1, na*P.nSubFrames*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes for each frame.
rcvI = 1;
for frameI = 1:Resource.RcvBuffer(1).numFrames
    Receive(rcvI).callMediaFunc = 1;
    acqI = 1;
    for subFrameI = 1:P.nSubFrames
        for angleI = 1:na
            Receive(rcvI).framenum = frameI;
            Receive(rcvI).acqNum = acqI;
            acqI = acqI + 1;
            rcvI = rcvI + 1;
        end
    end
end

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:na);

% Define ReconInfo structures.
% We need na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, na);
% - Set specific ReconInfo attributes.
if na>1
    ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
    for j = 1:na  % For each row in the column
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
    end
    ReconInfo(na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

fastSaveProcess = 2;
Process(fastSaveProcess).classname = 'External';
Process(fastSaveProcess).method = 'fastsave';
Process(fastSaveProcess).Parameters = { ...
    'srcbuffer','receive',... % name of buffer to export.
    'srcbufnum',    1,    ... % buffer number of buffer to export
    'srcframenum',  -1,   ... % process the last frame in a RcvBuffer
    'dstbuffer','none'
    };

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 160;  % 160 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 160;%20000 - (na-1)*160;  % 20 msec
SeqControl(4).command = 'returnToMatlab';

frameRateFactor = 5; % Factor for converting sequenceRate to frameRate.
nsc = length(SeqControl)+1; % Count of SeqControl objects.
n = 1; % Event structure arrays.
saveEvents = [];
acqI = 1;
for frameI = 1:Resource.RcvBuffer(1).numFrames
    for subFrameI = 1:P.nSubFrames
        for angleI = 1:na                      % Acquire frame
            Event(n).info = 'Full aperture.';
            Event(n).tx = angleI;
            Event(n).rcv = acqI;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = 2;
            n = n+1;
            acqI = acqI + 1;
        end
    end
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

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if floor(frameI/frameRateFactor) == frameI/frameRateFactor && frameI ~= Resource.RcvBuffer(1).numFrames
        Event(n).seqControl = 4;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;

% User specified UI Control Elements
import vsv.seq.uicontrol.VsSliderControl
import vsv.seq.uicontrol.VsToggleButtonControl

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7',...
    'Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
    'Callback',@SensCutoffCallback);

% - Range Change
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
MinMaxVal = [64,300,P.endDepth] ./ scaleToWvl;
AxesUnit = 'mm';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits) && ...
    strcmp(Resource.DisplayWindow(1).AxesUnits,'wavelengths')
    AxesUnit = 'wls';
    MinMaxVal = [64,300,P.endDepth];
end

UI(2).Control = VsSliderControl('LocationCode','UserA1',...
    'Label',['Range (',AxesUnit,')'],...
    'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],...
    'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

UI(3).Control = VsToggleButtonControl('LocationCode','UserB1',...
    'Label','FastSave',...
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
frameRate = 1/(((na - 1) * SeqControl(2).argument + SeqControl(3).argument ) * 1e-6);
superFrameRate = frameRate/P.nSubFrames;
avgDataFlow_GBpS = superFrameSize_GB*superFrameRate;  %< this should be less than 7GB/s, the maximum dma rate ???
maxNumFiles = floor(availableSpace_GB/superFrameSize_GB);
maxAcqDuration_s = maxNumFiles/superFrameRate;
maxAcqDuration_m = maxAcqDuration_s/60;
requiredSaveSpeed = superFrameSize_GB*superFrameRate;

disp("=====================================================")
fprintf(" Available Storage space: %g GB\n",availableSpace_GB);
fprintf(" SuperFrame Size: %.2g GB\n", superFrameSize_GB)
fprintf(" Acquisition Rate: %.2g fps\n", frameRate)
fprintf(" Super-Frame Rate: %.2g super-fps\n", superFrameRate)
fprintf(" Average Data Throughput: %.4g GB/sec\n", avgDataFlow_GBpS)
fprintf(" Max # of files: %g\n", maxNumFiles)
fprintf(" Max Acq. Duration: %i sec == %i minutes\n", round(maxAcqDuration_s), round(maxAcqDuration_m))
disp("=====================================================")

% Save all the structures to a .mat file.
save superFrameFrameRate.mat superFrameRate
save('MatFiles/L22-14vXFlashAngles_fastsave');

% **** Callback routines used by UI Controls. ****
function SensCutoffCallback(~,~,UIValue)
    ReconL = evalin('base', 'Recon');
    for i = 1:size(ReconL,2)
        ReconL(i).senscutoff = UIValue;
    end
    assignin('base','Recon',ReconL);
    Control = evalin('base','Control');
    if isempty(Control(1).Command), n=1; else, n=length(Control)+1; end
    Control(n).Command = 'update&Run';
    Control(n).Parameters = {'Recon'};
    assignin('base','Control', Control);
end

function RangeChangeCallback(hObject,~,UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No range change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.endDepth'));
        return
    end
    Trans = evalin('base','Trans');
    Resource = evalin('base','Resource');
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

    P = evalin('base','P');

    P.endDepth = UIValue * scaleToWvl;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits) && ...
       strcmp(Resource.DisplayWindow(1).AxesUnits,'wavelengths')
       P.endDepth = UIValue;
   end

    P.startDepthMm = P.startDepth/scaleToWvl;
    P.endDepthMm = P.endDepth/scaleToWvl;
    assignin('base','P',P);

    evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
    evalin('base','PData(1).Region = computeRegions(PData(1));');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength;
    end
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    Control = evalin('base','Control');
    if isempty(Control(1).Command), n=1; else, n=length(Control)+1; end
    Control(n).Command = 'update&Run';
    Control(n).Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end

function saveAcquisition(~,~,UIValue)
    saveEvents = evalin('base', 'saveEvents');
    Event = evalin('base', 'Event');
    for i = 1:length(saveEvents)
        Event(saveEvents(i)).process = UIValue*2;  %for external
    end
    assignin('base','Event',Event);
    Control = evalin('base','Control');
    if isempty(Control(1).Command), n=1; else, n=length(Control)+1; end
    Control(n).Command = 'update&Run';
    Control(n).Parameters = {'Event'};
    assignin('base','Control', Control);
end
