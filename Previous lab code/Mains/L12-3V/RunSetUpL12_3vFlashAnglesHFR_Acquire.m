% Notice: 
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility 
%   for its use
%
% File name: SetUpL12_3vFlashAngles.m - Example of 2-1 synthetic aperture plane
%                                       wave imaging with steering angle transmits
% Description: 
%   Sequence programming file for L12-3v Linear array, using 2-1 synthetic   
%   aperture plane wave transmits with multiple steering angles on a 128  
%   channel system. 128 transmit channels and 96 receive channels are  
%   active and positioned as follows (each char represents 4 elements) for  
%   each of the 2 synthetic apertures.
%
%   Element Nos.                                1              1
%                       3       6       9       2              9
%               1       3       5       7       9              2
%   Aperture : |        |       |       |       |              |
%               --------tttttttttttttttttttttttttttttttt-------
%               --------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr-------
%
%   The receive data from each of these apertures are stored under  
%   different acqNums in the Receive buffer. The reconstruction sums the 
%   IQ data from the 2 aquisitions and computes intensity values to produce 
%   the full frame. Processing is asynchronous with respect to acquisition.
%
% Last update:
% 12/06/2015 - modified for SW 3.0

clear all
cd C:\Users\BOAS-US\Documents\Vantage-3.4.0-1711281030

%% 1. Image parameters
P.startDepth = 60;   % Acquisition depth in wavelengths
P.endDepth = 140;   % This should preferrably be a multiple of 128 samples.
simulateMode = 1;   % acquire data with hardware
%  simulateMode = 1 forces simulate mode, even if hardware is present.
%  simulateMode = 2 stops sequence and processes RcvData continuously.
%%%% Define the super frame
P.CCangle=0; % imaging compounding angle, in degree
P.numAngles=1; % no. of angles for a CC US image
P.numCCframes = 200;      % no. of CC frames in a super frame
P.numFrames = 1;      % no. of super frames (real-time images are produced 1 per frame)
P.numAcqs =P.numAngles* P.numCCframes ;      % no. of Acquisitions in a "superframe"

if (P.numAngles > 1)
    dtheta = (P.CCangle*pi/180)/(P.numAngles-1);
    startAngle = -P.CCangle*pi/180/2;
else
    dtheta = 0;
    startAngle=0;
end % set dtheta to range over +/- P.CCangle/2 degrees.
%%%% RF data save name
CCFR=0.5; % khz
TintPlane=round(1000/P.numAngles/CCFR);
RFdataFilename = ['/Agl',num2str(P.CCangle),'-',num2str(P.numCCframes),'CC-',num2str(CCFR),'Khz-20v-nA',num2str(P.numAngles),'XZ40'];

%% 2. System parameters
filename = mfilename; % used to launch VSX automatically
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.connector = 1;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = simulateMode;

%% 3. Transducer parameters
Trans.name = 'L12-3v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans.frequency = 6;   % The center frequency for the A/D 4xFc sampling.
% note nominal center frequency in computeTrans is 7.813 MHz
Trans = computeTrans(Trans);  % L12-3v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

%% 4. PData structure array
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Resource.Parameters.numRcvChannels*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Resource.Parameters.numRcvChannels-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

%% 	5. Media object (for simulation)
pt1; % Specify Media object. 'pt1.m' script defines array of point targets.
Media.attenuation = -0.5;
Media.function = 'movePoints';

%% 	6. Maximum acquisition length and number of z samples
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Resource.Parameters.numRcvChannels-1)*Trans.spacing)^2));
maxZsamples=maxAcqLength *2*4; % *2 for round trip; *4 for 4 sampling points per wavelength

%% 	7. Resource parameters
%%%% ? Receive buffer
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = maxZsamples*P.numAcqs;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numFrames;       % number of 'super frames'
%%%% ? IQ buffer and Image buffer
Resource.InterBuffer(1).numFrames = 1;  % only one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = P.numFrames; % for reduced online visualization
%%%%%% Resource.ImageBuffer(1).numFrames = P.numFrames*P.numAcqs; % for all acquired data reconstruction
%%%% ? Display window parameters
Resource.DisplayWindow(1).Title = 'L12-3v-7AnglesFlashHFR';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
    DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
%%%%%% Resource.DisplayWindow(1).numFrames = P.numFrames*P.numAcqs; % for displaying all acquired images
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

%% 	8. Transmit waveform definition
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,20,1];

%% 	9. TX structures
%%%% ? Set TX structure array for P.numAngles transmits
TX = repmat(struct('waveform', 1, ...
        'Origin', [0.0,0.0,0.0], ...
        'aperture', 33, ... % use the center 128 elements for signal transmit
        'Apod', ones(1,Resource.Parameters.numTransmit), ...
        'focus', 0.0, ...
        'Steer', [0.0,0.0], ...
        'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, P.numAngles);
%%%% ? Set specific TX event for the P.numAngles transmits
if fix(P.numAngles/2) == P.numAngles/2       % if P.numAngles even
    P.startAngle = (-(fix(P.numAngles/2) - 1) - 0.5)*dtheta;
else
    P.startAngle = -fix(P.numAngles/2)*dtheta;
end
for iAngle = 1:P.numAngles   % P.numAngles transmit events
    TX(iAngle).Steer = [(P.startAngle+(iAngle-1)*dtheta),0.0];
    TX(iAngle).Delay = computeTXDelays(TX(iAngle));
end

%% 10. TGC Waveform structure.
TGC.CntrlPts = [139,535,650,710,770,932,992,1012];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% 11. Receive structures
%%%% ? Set RCV structure array for all acquisition
Receive = repmat(struct('Apod', ones(1,Resource.Parameters.numRcvChannels), ...
            'aperture', 33, ...
            'startDepth', P.startDepth, ...
            'endDepth', maxAcqLength, ...
            'TGC', 1, ...
            'bufnum', 1, ...
            'framenum', 1, ...
            'acqNum', 1, ...
            'sampleMode', 'NS200BW', ...
            'mode', 0, ...
            'callMediaFunc', 0),1,P.numAcqs*P.numFrames);
%%%% ? Set specific RCV event for all receive acquisitions
for iFrame = 1:Resource.RcvBuffer(1).numFrames
    Receive(P.numAcqs*(iFrame-1) + 1).callMediaFunc = 1;  % move points only once per super frame
    for iAcq = 1:P.numAcqs
        % -- Acquisitions for 'super' frame.
        rcvNum = P.numAcqs*(iFrame-1) + iAcq;
        Receive(rcvNum).Apod(:)=1;
%         Receive(rcvNum).callMediaFunc = 1;  % movepoints EVERY acquisition to illustrate superframe concept
        Receive(rcvNum).framenum = iFrame;
        Receive(rcvNum).acqNum = iAcq;
    end
end

%% 12. Recon structure arrays
%%%% Reconstuction for a super Frame
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:P.numAngles);
%%%% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, P.numAngles);
%%%% - Set specific ReconInfo attributes.
ReconInfo(1).mode = 'replaceIQ';
for iAngle = 1:P.numAngles
    ReconInfo(iAngle).txnum = iAngle;
    ReconInfo(iAngle).rcvnum = iAngle;
end
ReconInfo(P.numAngles).mode = 'accumIQ_replaceIntensity';

%% 13. Process structure array
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

%% 14. SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = TintPlane;  % 40 usec
SeqControl(3).command = 'triggerOut';
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'timeToNextAcq';  % time between super frames
SeqControl(5).argument = 1000000 - P.numAcqs*SeqControl(2).argument;  % 100000 usec = 100 msec
nsc = 6; % nsc is count of SeqControl objects

%% 15. Event objects to acquire all acquisitions
n = 1; % n is count of Events
for iSuperFrame = 1:Resource.RcvBuffer(1).numFrames
    for iCCframe = 1:P.numCCframes
        for iAngle=1:P.numAngles
            Event(n).info = 'Acquire RF';
            Event(n).tx = iAngle;
            Event(n).rcv = P.numAcqs*(iSuperFrame-1) + P.numAngles*(iCCframe-1)+iAngle;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = [2,3];
            n = n+1;
        end
    end
    % Set last acquisitions SeqControl for transferToHost.
    Event(n-1).seqControl = [5,3,nsc];
        SeqControl(nsc).command = 'transferToHost'; % transfer all acqs in one super frame
        nsc = nsc + 1;
    % Do reconstruction and processing for 1st sub frame
    Event(n).info = 'Reconstruct'; 
    Event(n).tx = 0;         
    Event(n).rcv = 0;        
    Event(n).recon = 1;      
    Event(n).process = 1;    
    Event(n).seqControl = 4;
    n = n+1;
end
% --- If this last event is not included, the sequence stops after one pass, and enters "freeze" state
%     Pressing the freeze button runs the "one-shot" sequence one more time
%     For live acquisition in mode 0, simply comment out the 'if/end' statements and manually freeze and exit when the data looks good.
if simulateMode==2 || simulateMode==0 %  In live acquisiton or playback mode, run continuously, but run only once for all frames in simulation
    Event(n).info = 'Jump back to first event';
    Event(n).tx = 0;        
    Event(n).rcv = 0;       
    Event(n).recon = 0;     
    Event(n).process = 0;
    Event(n).seqControl = 1; 
end


%$ 16. User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = P.numAcqs;

%% 17. save setup files and run VSX
save(['MatFiles/',filename]); 
disp([ mfilename ': NOTE -- Running VSX automatically!']), disp(' ')
VSX    
commandwindow  % just makes the Command window active to show printout

%% 18. after quitting VSX, save the RFData collected at high frame rate 
if simulateMode ~= 2    
    disp ('Info:  Saving the RF Data buffer -- please wait!'), disp(' ')
    save (['MatFiles/',RFdataFilename], '-v7.3','RcvData', 'P')
end
return

%% 19. **** Callback routines to be converted by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%RangeChangeCallback - Range change
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
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;    
    end
end
assignin('base','P',P);

evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Resource.Parameters.numRcvChannels-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback
