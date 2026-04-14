% File name: SetUpGEL8_18iD_FlashAngles_PulseInvDualDisp.m - Example of 
%            harmonic imaging using pulse inversion with multiple plane 
%            waves.
%
% Description:
%   This script provides an example of Pulse Inversion (PI) imaging using 
%   multiple angle unfocused waves for contrast agent imaging or tissue 
%   harmonic imaging. A broadband probe is required in order to transmit a 
%   one frequency and receive at twice that frequency. This script also
%   demostrates the capability of the system to accumulate RF signals in
%   the hardware by using the Receive.mode attribute. 
%
%   This script provides the posibility to do Bmode imaging or pulse
%   inversion imaging through a GUI control in VSX_Gui. 
%   For the pulse inversion sequence the script performs two transmits  
%   pulses with opposite polarity. Both pulses are done at TWfreq. During
%   the receive state, the two pulses are added in the hardware using the 
%   Receive.mode=1. This essentially cancels all the linear contributions 
%   coming back from the medium leaving the nonlinear component. 
%   The Bmode imaging is done at 2xTWfreq to be able to compare with the
%   pulse inversion.  
%   The scripts transmit and receives using all 168 channels available on 
%   the GEL8-18iD probe 
%
%   This sequence was developed and tested under the following conditions:
%       VDAS Configuration: Vantage NXT 256 HF System (Bmode Imaging
%       TX-freq:13.3 MHz, PI TX-freq:6.667 MHz )
%       .  Software Configuration: VantageNXT 2.0.0
%       .  Transducer: GE8-18iD
%       .  Matlab Version: R2022b
%
% Notice:
%   This file is based on a collaboration between the Ultrasound Lab of 
%   Polytechnique Montreal led by Jean Provost, and Verasonics Inc.  
%   This script servers users as a programming example for the Verasonics 
%   NXT Research Ultrasound System. Polytechnique or Verasonics makes no 
%   claims as to the functionality or intended application of this program 
%   and the user assumes all responsibility for its use.
% 
% Copyright © 2024 Verasonics, Inc.

clear all


P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 128;   % This should preferrably be a multiple of 128 samples.
na = 11;             % Set na = number of angles.
angleRange = 12;    % total angular range for plane wave angles (% set dtheta to range over +/- 18 degrees.)
Polarity =  -1;     % +1 or -1: Polarity adjusts the sign of the second pulse with respect to the first one

TWfreq = 6.667;       % transmit frequency (set lower than Trans.frequency for harmonic imaging)

if (na > 1)
    dtheta = (angleRange*pi/180)/(na-1);
else
    dtheta=0;
end

% Define system parameters.
Resource.Parameters.numTransmit = 256;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;

%% Specify Media object.
pt1;
Media.function = 'movePoints';

%% Specify Trans structure array.
Trans.name = 'GEL8-18iD';
Trans.units = 'wavelengths';  % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % GEL8-18iD transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 30;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [0.5, 0, 0.25];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

%% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = na*3072; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;    % 30 frames stored in RcvBuffer.
Resource.RcvBuffer(2).rowsPerFrame = 2*na*3072; % this size allows for maximum range
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = 10;    % 30 frames stored in RcvBuffer.

Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;

Resource.InterBuffer(2).datatype = 'complex';
Resource.InterBuffer(2).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(2).datatype = 'double';
Resource.ImageBuffer(2).numFrames = 10;

Resource.DisplayWindow(1).Title = 'GEL8-18iDFlashAngles_PulseInv';...
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [770,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).numFrames = 100;

Resource.DisplayWindow(2).Title = 'GEL8-18iDFlashAngles_Bmode';
Resource.DisplayWindow(2).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta);
Resource.DisplayWindow(2).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(2).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).AxesUnits = 'mm';
Resource.DisplayWindow(2).Colormap = gray(256);
Resource.DisplayWindow(2).numFrames = 100;
%% Specify Transmit waveform structure.
% Waveform for Bmode imaging twice the frequency used for Pulse Inversion
TW(1).type = 'parametric';
TW(1).Parameters = [2*TWfreq,1,1,1];   % Waveform for Bmode
TW(1).equalize = 1;

% use a relatively low frequency (default 5 MHz) to be better able to detect the second harmonic (default 10 MHz)
TW(2).type = 'parametric';
TW(2).Parameters = [TWfreq,1,2,1];   % positive polarity
TW(2).equalize = 0;

TW(3) = TW(2);
TW(3).Parameters(4) = Polarity*TW(1).Parameters(4);   % negative polarity
TW(3).equalize = 0;

%% Specify TX structure array.
%                  'Apod', kaiser(Trans.numelements)', ...
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, 2*na+na);
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
        startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
    else
        startAngle = -fix(na/2)*dtheta;
end

for n = 1:na  % Transmit for Bmode imaging
    TX(n).Steer = [(startAngle+(n-1)*dtheta),0.0]; % TW(1) for Bmode Imaging
    TX(n).Delay = computeTXDelays(TX(n));
end

j=0;
for n = na+2:2:3*na   % 2*na transmit events
    TX(n-1).Steer = [(startAngle + j*dtheta),0.0]; % TW(2) for Pulse Inversion (+ polarity)
    TX(n-1).Delay = computeTXDelays(TX(n-1));
    TX(n-1).waveform = 2;

    TX(n).Steer = TX(n-1).Steer;                   % TW(3)for Pulse Inversion (- polarity)
    TX(n).Delay = TX(n-1).Delay;
    TX(n).waveform = 3;
    j=j+1;
end

TPC(2).hv=1.6;

%% Specify Receive structure arrays.
% - We need 3*na Receives for every frame (na for Bmode and 2*na for Pulse inversion).
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'demodFrequency', TWfreq*2, ...
                        'callMediaFunc', 0), 1, 3*na*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = na*(i-1);
    Receive(k+1).callMediaFunc = 1;    % only move media points once per frame, not between angles
    for j = 1:na
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    end
end

lastBmodeRcv=k+j;

for i = 1:Resource.RcvBuffer(1).numFrames
    k = lastBmodeRcv+2*na*(i-1);
    Receive(k+1).callMediaFunc = 1;    % only move media points once per frame, not between angles

    m = 1;
    for j = 1:2:2*na
        Receive(k+j).Apod(1:Trans.numelements) = 1.0; % TW(2)
        Receive(k+j).framenum = i;
        Receive(k+j).bufnum = 2;
        Receive(k+j).acqNum = m;
        Receive(k+j+1).Apod(1:Trans.numelements) = 1.0; % TW(3)
        Receive(k+j+1).framenum = i;
        Receive(k+j+1).bufnum = 2;
        Receive(k+j+1).acqNum = m;
        Receive(k+j+1).mode = 1;
        m = m + 1;
    end
end




%% Specify TGC Waveform structure.
TGC.CntrlPts =  [535 1023 1023 1023 1023 1023 1023 1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Recon structure arrays.
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

Recon(2) = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [2,1], ...
               'ImgBufDest', [2,-1], ...
               'RINums',na+1:2*na);

% Define ReconInfo structures.
% We need na ReconInfo structures for given that the pulse inversion acquisitions were accumulated in HW
ReconInfo(na+1:2*na) = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, na);

% - Set specific ReconInfo attributes.
ReconInfo(na+1).mode = 'replaceIQ'; % replace IQ data
jj=0;
for j = 1:na  
    ReconInfo(j+na).txnum = j+na+jj;     
    ReconInfo(j+na).rcvnum = j+lastBmodeRcv+jj;
    jj=jj+1;
end
if na>1
    ReconInfo(2*na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo(2*na).mode = 'replaceIntensity';
end


%% Specify Process structure array.
pers = 0;
%Process(1)  is used for the Bmode,  however is displayed on the
%DisplayWindow(2) in order for the cineloop to save the Pulse Inversion
%data
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',63,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','log',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',2};

Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',88.0,...            % pgain is image processing gain
                         'reject',0,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','log',...
                         'compressFactor',50,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};


Process(3).classname = 'External';
Process(3).method = 'UIControl';
Process(3).Parameters = {'srcbuffer','none'};

%% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 2;

SeqControl(2).command = 'timeToNextAcq';  % time between PI acquisitions
twoWayTravel = 2.5 * Receive(1).endDepth  / Trans.frequency; % microsecs
SeqControl(2).argument = 100;...twoWayTravel;

SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 50000;  % 50 msec
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).condition = 'immediate';
SeqControl(5).argument = 1;
SeqControl(6).command = 'setTPCProfile';
SeqControl(6).condition = 'immediate';
SeqControl(6).argument = 2;


nsc = 7; % nsc is next count of SeqControl objects

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

%% Specify Event structure arrays.

n = 1;

Event(n).info = 'Set the Labels for the Voltage controls';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 3;
Event(n).seqControl = 0;
n = n+1;

% Switch to TPC profile 1 Bmode
Event(n).info = 'Switch to profile 1.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 5;
n = n+1;

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

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame
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
n=n+1;

%%

nStartPI = n;

% Switch to TPC profile 2 Pulse Inversion
Event(n).info = 'Switch to profile 1.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 6;
n = n+1;

for i = 1:Resource.RcvBuffer(1).numFrames
    k = lastBmodeRcv+2*na*(i-1);
    for j = 1:2:2*na                      % Acquire frame
        Event(n).info = 'First Waveform';
        Event(n).tx = na+j;
        Event(n).rcv = k+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Second Waveform';
        Event(n).tx = na+j+1;
        Event(n).rcv = k+j+1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 2;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame
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
      SeqControl(nsc).command = 'jump'; % 
      SeqControl(nsc).argument = nStartPI;
      nsc = nsc+1;
Event(n).seqControl = nsc-1;



%% User specified UI Control Elements

import vsv.seq.uicontrol.VsSliderControl
import vsv.seq.uicontrol.VsButtonControl

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                   'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                   'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
                   'Callback', @SensCutoffCallback);


% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                    'SliderMinMaxVal',MinMaxVal,...
                    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
                    'Callback', @RangeChangeCallback);

UI(4).Control  = VsButtonControl( 'LocationCode',   'UserB4', 'Label','Bmode ','Callback', @BmodeStart );
UI(5).Control  = VsButtonControl( 'LocationCode',   'UserB4', 'Label','Pulse Inv ','Callback', @PIStart );

EF(1).Function = vsv.seq.function.ExFunctionDef('UIControl', @UIControl);


% Save all the structures to a .mat file.
save('MatFiles/GEL8-18iDFlashAngles_PulseInvDualDisp');

%% **** Callback routines used by UIControls (UI) ****

function SensCutoffCallback(~,~,UIValue)
%Sensitivity cutoff change
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
%Range change
    UI = evalin('base','UI');
    set(UI(2).handle,'Interruptible','off');

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
    P.endDepth = UIValue;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            P.endDepth = UIValue*scaleToWvl;
        end
    end
    assignin('base','P',P);

    evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
    evalin('base','PData(1).Region = computeRegions(PData(1));');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
     
    evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
    evalin('base','PData(1).Region = computeRegions(PData(1));');
    evalin('base','Resource.DisplayWindow(2).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta);');
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength;
    end
    twoWayTravel = 2.5 * Receive(1).endDepth  / Trans.frequency; % microsecs
    SeqControl = evalin('base', 'SeqControl');
    SeqControl(2).argument = twoWayTravel;  % formerly 160 usec
    assignin('base','SeqControl',SeqControl);
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    Control = evalin('base','Control');
    if isempty(Control(1).Command), n=1; else, n=length(Control)+1; end
    Control(n).Command = 'update&Run';
    Control(n).Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon','SeqControl'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end

function PIStart(varargin)
    UI = evalin('base','UI');
 
    set(UI(4).handle,'Visible','on'); %UI(4) Bmode
    set(UI(5).handle,'Visible','off');  %UI(5) PI

    nStart = evalin('base','nStartPI');
    Control = evalin('base','Control');
    
    if isempty(Control(1).Command)
        n=1; 
    else 
        n=length(Control)+1;
    end

    Control(n).Command = 'set&Run';
    Control(n).Parameters = {'Parameters',1,'startEvent',nStart};
    evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
    assignin('base','Control',Control);

end

function BmodeStart(varargin)
    UI = evalin('base','UI');
    
    set(UI(5).handle,'Visible','on');   %UI(5) PI
    set(UI(4).handle,'Visible','off'); %UI(4) Bmode

    nStart = 2; % needs to be event 2 because Event(1) is to set the labels on the voltage buttons
    Control = evalin('base','Control');
    
    if isempty(Control(1).Command)
        n=1; 
    else 
        n=length(Control)+1;
    end

    Control(n).Command = 'set&Run';
    Control(n).Parameters = {'Parameters',1,'startEvent',nStart};
    evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
    assignin('base','Control',Control);

end

%% **** Callback routines used by External function definition (EF) ****

function UIControl(varargin)

if isempty(findobj('String','Bmode Voltage'))==1
    hv1Handle(1) = findobj('String','High Voltage P1');
    set(hv1Handle,'String','Bmode Voltage');
   
    hv2Handle(1) = findobj('String','High Voltage P2');
    set(hv2Handle,'String','PI Voltage');
end
end
