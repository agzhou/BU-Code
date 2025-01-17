% File name: RunSetUpL12_3VPAT_Acquire_US_MASTER.m
% Description: 
%   Sequence programming file for L12-3v Linear array, using 2-1 synthetic   
%   aperture plane wave transmits with multiple steering angles on a 128  
%   channel system. 128 transmit channels and 96 receive channels are  
%   active and positioned as follows (each char represents 4 elements) for  
%   each of the 2 synthetic apertures.
%
%   The receive data from each of these apertures are stored under  
%   different acqNums in the Receive buffer. The reconstruction sums the 
%   IQ data from the 2 aquisitions and computes intensity values to produce 
%   the full frame. Processing is asynchronous with respect to acquisition.
%
% Last update:
% 20190929 by Bingxue Liu

clear all

QSwitch_delay_us = 300; % Should be >= 300 for Laser External trigger + QSwitch external trigger

if QSwitch_delay_us < 300
    fprintf("QSwitch_delay_us = %d which is less than 300!\nThis would damage the laser. ABORTING!\n", QSwitch_delay_us);
    return
end

trigger = true;
transmit = false;

%%%%
cd C:\Users\BOAS-US\Documents\Vantage-3.4.0-1711281030
% addpath('D:\CODE\Functions');
load('D:\CODE\Mains\PAT\PATAQParameters_L12_3v.mat');
deftpath=PATAQInfo.savepath;
savepath=[uigetdir(deftpath),'\'];
% savepath=deftpath;
%% 1. PATAQ parameter 
prompt={'File name:',...
    'Average number', ...
    'Receive Start Depth, [Lamda]:',...
    'Receive End Depth, [Lamda]:',...
    'Sampling Mode (50, 67, 100, 200 [%])',...
    'Mode'};
name='PAT AQ parameters';
defaultvalue={PATAQInfo.Filename,...
    num2str(PATAQInfo.AverageNumber), ...
    num2str(PATAQInfo.StartDepth),...
    num2str(PATAQInfo.EndDepth),...
    num2str(PATAQInfo.SplMode) ...
    ,'0'};
inputValue=inputdlg(prompt,name, 1, defaultvalue);
P.Filename=inputValue{1};
P.AverageNumber=str2num(inputValue{2});
P.StartDepth=str2num(inputValue{3});
P.EndDepth=str2num(inputValue{4});
P.SplMode=str2num(inputValue{5});
P.simulateMode =str2num(inputValue{6});
if P.SplMode==200
    P.SampleMode='NS200BW';
    P.nSperWave=4; % sampling points per wavelength
elseif P.SplMode==50
    P.SampleMode='BS50BW';
    P.nSperWave=1; % sampling points per wavelength
elseif P.SplMode==67
    P.SampleMode='BS67BW';
    P.nSperWave=1.3333; % sampling points per wavelength
else
    P.SampleMode='BS100BW';
    P.nSperWave=2; % sampling points per wavelength
end
%% Save parameters

PATAQInfo=P;
PATAQInfo.savepath = savepath;
save('D:\CODE\Mains\PAT\PATAQParameters.mat_L12_3v','PATAQInfo');

%%%% RF data save name
RFdataFilename = ['RF-',P.Filename,'-1'];
Filename = P.Filename;
%% 2. System parameters
filename = mfilename; % used to launch VSX automatically
Resource.Parameters.numTransmit = 128;     % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;    % number of receive channels. CHANGE HERE
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.connector = 2;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = P.simulateMode;

%Resource.VDAS.halDebugLevel = 1;

P.nCh=Resource.Parameters.numRcvChannels ;
P.vSound=Resource.Parameters.speedOfSound;
% RcvProfile.condition = 'immediate';
%% 3. Transducer parameters
Trans.name = 'L12-3v';
Trans.units = 'wavelengths';
% Trans.frequency = 15;   % The center frequency for the A/D 4xFc sampling.
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 50; % mfr data sheet lists 30 Volt limit
P.LensDelay=Trans.lensCorrection;
P.frequency=Trans.frequency;
P.pitch=Trans.spacingMm;
%% 4. PData structure array

PData.PDelta = [Trans.spacing/4, 0, 0.25];
PData(1).Size(1) = ceil((P.EndDepth-P.StartDepth+1)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Resource.Parameters.numRcvChannels*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Resource.Parameters.numRcvChannels-1)/2,0,P.StartDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

%% 	5. Media object (for simulation)
pt1; % Specify Media object. 'pt1.m' script defines array of point targets.
Media.attenuation = -0.5;
Media.function = 'movePoints';
% Media.MP(1,:) = [0,0,70,1]; % Single point, [x, y, z, reflectivity]


%% 	6. Maximum acquisition length and number of z samples
maxAcqLength = ceil(sqrt(P.EndDepth^2 + ((Resource.Parameters.numRcvChannels-1)*Trans.spacing)^2));
actAcqLength =(maxAcqLength-P.StartDepth/2); % P.StarDepth/2
P.actZsamples=ceil(actAcqLength * 2 * P.nSperWave/128)*128; % P.nSperWave: sampling points per wavelength, rounded up to the next 128 VSX sample boundary


%% 	7. Resource parameters
% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.actZsamples*P.AverageNumber; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 1;       % 1 frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 1;

% Display
Resource.DisplayWindow(1).Title = 'L12-3v PAT Acquisition';    %%%%%, 4X sampling at 62.5 MHz';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 1;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

%% 	8. Transmit waveform definition
TW.type = 'parametric';
TW.Parameters = [6, 0.67, 2, 1];   % 6 MHz center frequency, 67% pulsewidth two cycle burst

%% 	9. TX structures
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(1).focus = 0;
TX(1).aperture = 33; 
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
if transmit
    TX(1).Apod = ones(1,Resource.Parameters.numTransmit);
else
    TX(1).Apod = zeros(1,Resource.Parameters.numTransmit);
end
TX(1).Delay = computeTXDelays(TX(1));

%%  10. Analog front end gain settings.

RcvProfile.LnaGain = 24;
RcvProfile.condition = 'immediate';
%%  11. TGC Waveform structure. --> UNIFORM??
TGC(1).CntrlPts = [500,590,650,710,770,830,890,950];    % Standard
% TGC(1).CntrlPts = round(950.^linspace(0,1,8));    % Exponential
TGC(1).rangeMax = P.EndDepth;
TGC(1).Waveform = computeTGCWaveform(TGC);

%% 12. Receive structures
%%%% ? Set RCV structure array for all acquisition
% - We need one Receive for every frame.
% sampling center frequency is 15.625, but we want the bandpass filter
% centered on the actual transducer center frequency of 18 MHz with 67%
% bandwidth, or 12 to 24 MHz.  Coefficients below were set using
% "filterTool" with normalized cf=1.15 (18 MHz), bw=0.85, 
% xsn wdth=0.41 resulting in -3 dB 0.71 to 1.6 (11.1 to 25 MHz), and
% -20 dB 0.57 to 1.74 (8.9 to 27.2 MHz)
% 
% BPF1 = [ -0.00009 -0.00128 +0.00104 +0.00085 +0.00159 +0.00244 -0.00955 ...
%          +0.00079 -0.00476 +0.01108 +0.02103 -0.01892 +0.00281 -0.05206 ...
%          +0.01358 +0.06165 +0.00735 +0.09698 -0.27612 -0.10144 +0.48608 ];
 partial = [zeros(1,64),ones(1,64)];    
Receive = repmat(struct(   'Apod', partial, ... %   CHANGE HERE ones(1,Resource.Parameters.numRcvChannels), ...
                    'aperture', 1, ...  % CHANGE HERE 33
                    'startDepth', P.StartDepth/2, ...  %P.StartDepth/2, ... 
                    'endDepth',   P.StartDepth/2 + actAcqLength, ...  %maxAcqLength,... %  %
                    'TGC', 1, ...
                    'bufnum', 1, ...
                    'framenum', 1, ...
                    'acqNum', 1, ...
                    'sampleMode', P.SampleMode, ...
...%                     'InputFilter', BPF1, ... 
                    'mode', 0, ...
                    'callMediaFunc', 0), 1, P.AverageNumber);
                
Receive(1).callMediaFunc = 1;

for iAcq = 1:P.AverageNumber
    Receive(iAcq).acqNum = iAcq;
end

%% 13. Recon structure arrays
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
...%                'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:P.AverageNumber ,...
               'newFrameTimeout', -1);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, P.AverageNumber);

ReconInfo(1).mode = 'replaceIQ';

for iAcq=1:P.AverageNumber
    ReconInfo(iAcq).rcvnum = iAcq;
end

if P.AverageNumber == 1
    ReconInfo(P.AverageNumber).mode = 'replaceIntensity';
else
    ReconInfo(P.AverageNumber).mode = 'accumIQ_replaceIntensity';
end

%% 14. Process structure array

pers = 0;%P.AverageNumber;

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

%% 15. Event objects to acquire all acquisitions
num_evt = 1;
start_evt = num_evt;

for iAcq=1:P.AverageNumber
    Event(num_evt).info = 'Start Acquisition';
    Event(num_evt).tx = 1;              % no TX structure.
    Event(num_evt).rcv = iAcq;          % use 1st Rcv structure.
    Event(num_evt).recon = 0;           % no reconstruction.
    Event(num_evt).process = 0;         % no processing
    Event(num_evt).seqControl = [2,3]; %  wait for trigger in to acquire + send trigger Out %+ transfer RF data to host
    num_evt = num_evt + 1;
end

Event(num_evt).info = 'Custom Sw & Hw sync timeout.';   % sync timeout is .5s by default
Event(num_evt).tx = 0;            % no TX structure.
Event(num_evt).rcv = 0;           % no Rcv structure.
Event(num_evt).recon = 0;         % no reconstruction.
Event(num_evt).process = 0;       % no processing
Event(num_evt).seqControl = 6;    % make the Sw wait at most 60s for the Hw to sync
num_evt = num_evt + 1;

Event(num_evt).info = 'Transfer RF Data to Host.';
Event(num_evt).tx = 0;            % no TX structure.
Event(num_evt).rcv = 0;           % no Rcv structure.
Event(num_evt).recon = 0;         % no reconstruction.
Event(num_evt).process = 0;       % no processing
Event(num_evt).seqControl = 4;    % transfer to host
num_evt = num_evt + 1;

% Perform reconstruction
Event(num_evt).info = 'Reconstruct RF Data.';
Event(num_evt).tx = 0;            % no TX structure.
Event(num_evt).rcv = 0;           % no Rcv structure.
Event(num_evt).recon = 1;         % reconstruction.
Event(num_evt).process = 1;       % reconstruction processing
Event(num_evt).seqControl = 5;    % return to matlab
num_evt = num_evt + 1;

% Jump back to Acquisition
Event(num_evt).info = 'Jump back to acquisition.';
Event(num_evt).tx = 0;            % no TX structure.
Event(num_evt).rcv = 0;           % no Rcv structure.
Event(num_evt).recon = 0;         % reconstruction.
Event(num_evt).process = 0;       % no processing
Event(num_evt).seqControl = 1;


%% 16. SeqControl structure arrays.

SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = start_evt;

% SeqControl(2).command = 'triggerIn';  % time between frames
% SeqControl(2).condition = 'Trigger_2_Rising';  % Wait for trigger 2 rising edge
if trigger
    SeqControl(2).command = 'triggerIn';  % time between frames
    SeqControl(2).condition = 'Trigger_2_Rising';  % Wait for trigger 1 rising edge
%     SeqControl(2).argument = 0;
else
    SeqControl(2).command = 'timeToNextAcq';  % time between frames
    SeqControl(2).argument = 50000;  % 50 msec
end

SeqControl(3).command = 'triggerOut';
SeqControl(3).argument = 0; %565;

SeqControl(4).command = 'transferToHost';

SeqControl(5).command = 'returnToMatlab';

SeqControl(6).command = 'sync';
SeqControl(6).argument = 60e6; % Wait 60 seconds

%% 17. User specified UI Control Elements
% % - Sensitivity Cutoff
% UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
%                   'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
%                   'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
% UI(1).Callback = text2cell('%SensCutoffCallback');
% 
% % - Range Change
% % - Range Change
% MinMaxVal = [64,300,EndDepth]; % default unit is wavelength
% AxesUnit = 'wls';
% if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
%     if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
%         AxesUnit = 'mm';
%         MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
%     end
% end
% UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
%                  'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
% UI(2).Callback = text2cell('%RangeChangeCallback');
% 
% % Specify factor for converting sequenceRate to frameRate.
% frameRateFactor = P.numAcqs;

%% 18. save setup files and run VSX
save(['MatFiles/',filename]); 
disp([ mfilename ': NOTE -- Running VSX automatically!']), disp(' ')
VSX    
%% 19. after quitting VSX, save the RFData collected at high frame rate 
if P.simulateMode ~= 2    
    disp ('Info:  Saving the RF Data buffer -- please wait!'), disp(' ')
    if exist([savepath,RFdataFilename,'.mat'])
        NewPlane=1;
        while (exist([savepath,RFdataFilename,'.mat'])==2)
            NewPlane=NewPlane+1;
            RFdataFilename = ['RF-IQ-',P.Filename,'-',num2str(NewPlane)];
        end
        prompt={'File exist! New Plane Sequence: '};
        name='File info';
        defaultvalue={num2str(NewPlane)};
        inputValue=inputdlg(prompt,name, 1, defaultvalue);
        P.Plane=str2num(inputValue{1});
        RFdataFilename = ['RF-IQ-',P.Filename,'-',num2str(P.Plane)];
    end
    tic
    RFRAW=RcvData{1};Img=ImgData{1};
    savefast ([savepath,RFdataFilename], 'RFRAW','P','Img')
    toc
    disp ('RF DATA SAVED!')
end
return
