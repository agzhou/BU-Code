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
%   Sequence programming file for L22-14v Linear array, using 2-1 synthetic   
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
% 2-28-2018 by Jianbo Tang

clear all

flash2Qdelay = 300; % Should be >= 300 microsecs for Laser External trigger + QSwitch external trigger

if flash2Qdelay < 300
    fprintf("QSwitch_delay_us = %d which is less than 300!\nThis would damage the laser. ABORTING!\n", flash2Qdelay);
    return
end

PA_PRF = 20; % PA PRF in Hz. To be set in accordance with laser rep rate, when not using the input trigger mode.
                % When using the input trigger mode, remove the TTNA (SeqControl(4)) from the PA events to avoid
                % "missed TTNA" messages.

%%%%
cd C:\Users\BOAS-US\Desktop\Vantage-3.4.3-1807021300  %C:\Users\BOAS-US\Documents\Vantage-3.4.0-1711281030
addpath('D:\CODE\Functions');

transmit = 0;     % (logical) transmit=0 turns off the transmitters by setting TX.Apod to zero for all transmitters
    if transmit==0
        disp(' *** PhotoAcoustic mode: Using one-way receive-only reconstruction ***')
    else
        disp(' *** Ultrasound Transmit mode: Using conventional T/R reconstruction ***')
    end
    
load('D:\CODE\Mains\DAQParameters.mat');
deftpath=DAQInfo.savepath;
savepath=[uigetdir(deftpath),'\'];


%% 1.2 US DAQ parameter
prompt={'File name:',...
    'SuperFrame/PDI rate(PDIFR), [Hz]',...
    'Number of SuperFrame (nSupFrame)', ...
    'Coherence Compounding frame rate(CCFR), [Hz]:',...
    'Number of Coherence Compounding frames(nCC):',...
    'Compounding Angle (Angle), [degree]:',...
    'Number of planes for CC (nAngle):',...
    'Receive Start Depth, [Lamda]:',...
    'Receive End Depth, [Lamda]:',...
    'TW frequency [MHz]',...
    'TW number of HalfCycle',...
    'Sampling Mode (50, 67, 100, 200 [%])',...
    'Mode'};
name='Ultrafast US DAQ parameters';
defaultvalue={DAQInfo.Filename,...
    num2str(DAQInfo.PDIFR),...
    num2str(DAQInfo.numSupFrames),...
    num2str(DAQInfo.CCFR),...
    num2str(DAQInfo.numCCframes),...
    num2str(DAQInfo.CCangle),...
    num2str(DAQInfo.numAngles),...
    num2str(DAQInfo.StartDepth),...
    num2str(DAQInfo.EndDepth),...
    '18', '3','100','0'};
inputValue=inputdlg(prompt,name, 1, defaultvalue);
Filename=inputValue{1};
P.PDIFR=str2num(inputValue{2});
P.numSupFrames=str2num(inputValue{3});
P.CCFR=str2num(inputValue{4});
P.numCCframes=str2num(inputValue{5});
P.CCangle=str2num(inputValue{6});
P.numAngles=str2num(inputValue{7});
P.StartDepth=str2num(inputValue{8});
P.EndDepth=str2num(inputValue{9});
P.TWfrequency=str2num(inputValue{10});
P.TWnHC=str2num(inputValue{11});
P.SplMode=str2num(inputValue{12});
simulateMode =str2num(inputValue{13});
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
%% time intervale between plane wave emission
P.t2NextPlaneDAQ=33;    % us, time interval between successive plane wave, travel time + system overhead time, the default is 33 us for 150 wavelengths acquisition
tCCUS=round(1e6/P.CCFR); % Time intervale between successive CCUS frame,us 
P.tIntPDI=1e6/P.PDIFR;

if P.tIntPDI-tCCUS*P.numCCframes<10 || (P.tIntPDI-tCCUS*P.numCCframes)>4e6
    P.tIntPDI=tCCUS*(P.numCCframes)+0.1e6;
    P.PDIFR=1e6/P.tIntPDI;
    prompt={['Incorrect superFrame rate, Approriate Range: ', num2str(1e6/(tCCUS*P.numCCframes+4e6)) ,'-',num2str(1e6/(tCCUS*(P.numCCframes)+0.1e6))]};
    name='Wrong superFrame rate';
    defaultvalue={num2str(P.PDIFR)};
    inputValue=inputdlg(prompt,name, 1, defaultvalue);
    P.PDIFR=str2num(inputValue{1});
    P.tIntPDI=1e6/P.PDIFR;
end
DAQInfo=P;
DAQInfo.StartDepth=P.StartDepth;
DAQInfo.EndDepth=P.EndDepth;
DAQInfo.Filename=Filename;
DAQInfo.savepath=savepath;
P.Filename=Filename;
save('D:\CODE\Mains\DAQParameters.mat','DAQInfo');

%%%% RF data save name
RFdataFilename = ['RF-',num2str(P.CCangle),'-',num2str(P.numAngles),'-',num2str(P.CCFR),'-',num2str(P.numCCframes),'-',num2str(P.numSupFrames),'-',Filename,'-', '1'];
%% time intervale between plane wave emission
P.t2NextPlaneDAQ=33;    % us, time interval between successive plane wave, travel time + system overhead time, the default is 33 us for 150 wavelengths acquisition
tCCUS=round(1e6/P.CCFR); % Time intervale between successive CCUS frame,us 
% P.Plane=1;      % current imaging plane
% P.numSupFrames = 1;      % no. of super frames (real-time images are produced 1 per frame)
% P.PDIFR=2;      % Power Doppler frame rate
% P.numCCframes = 200;      % no. of CC frames in a super frame
% P.CCFR=1000;    % Coherence compounding frame rate, hz
% P.CCangle=20;   % imaging compounding angle, in degree
% P.numAngles=11; % no. of angles for a CC US image
% StartDepth = 30;   % Acquisition depth in wavelengths
% EndDepth = 150;    % This should preferrably be a multiple of 128 samples.
%%%%%%%%%%%% select working mode %%%%%%%%%%%%%%%%%%%%
%  simulateMode = 0 acquire data with hardware
%  simulateMode = 1 forces simulate mode, even if hardware is present.
%  simulateMode = 2 stops sequence and processes RcvData continuously.
%%%% Define the super frame
P.numAcqs =P.numAngles* P.numCCframes ;      % no. of Acquisitions in a "superframe"
if (P.numAngles > 1)
    dtheta = (P.CCangle*pi/180)/(P.numAngles-1);
    startAngle = -P.CCangle*pi/180/2;
else
    dtheta = 0;
    startAngle=0;
end % set dtheta to range over +/- P.CCangle/2 degrees.


%% 1.1 PA DAQ parameter 

load('D:\CODE\Mains\PAT\PATAQParameters_L22_14v.mat');
deftpath_pa=PATAQInfo.savepath;
savepath_pa=[uigetdir(deftpath_pa),'\'];
% savepath=deftpath;

prompt={'File name:',...
    'Average number', ...
    'Number of wavelength(nSuperFrame)',...
    'Receive Start Depth, [Lamda]:',...
    'Receive End Depth, [Lamda]:',...
    'Sampling Mode (50, 67, 100, 200 [%])',...
    'Mode'};
name='PAT AQ parameters';
defaultvalue={PATAQInfo.Filename,...
    num2str(PATAQInfo.AverageNumber), ...
    '2', ...
    num2str(PATAQInfo.StartDepth),...
    num2str(PATAQInfo.EndDepth),...
    num2str(PATAQInfo.SplMode) ...
    ,'0'};
inputValue=inputdlg(prompt,name, 1, defaultvalue);
P(2).Filename=inputValue{1};
P(2).AverageNumber=str2num(inputValue{2});
P(2).numSupFrames=str2num(inputValue{3});
P(2).StartDepth=str2num(inputValue{4});
P(2).EndDepth=str2num(inputValue{5});
P(2).SplMode=str2num(inputValue{6});
P(2).simulateMode =str2num(inputValue{7});
if P(2).SplMode==200
    P(2).SampleMode='NS200BW';
    P(2).nSperWave=4; % sampling points per wavelength
elseif P(2).SplMode==50
    P(2).SampleMode='BS50BW';
    P(2).nSperWave=1; % sampling points per wavelength
elseif P(2).SplMode==67
    P(2).SampleMode='BS67BW';
    P(2).nSperWave=1.3333; % sampling points per wavelength
else
    P(2).SampleMode='BS100BW';
    P(2).nSperWave=2; % sampling points per wavelength
end

P(2).numAcqs = P(2).AverageNumber*P(2).numSupFrames; % no. of all PA acquisitions for all supframes

%%%%%% Save parameters
PATAQInfo=P(2);
PATAQInfo.savepath = savepath(2);
save('D:\CODE\Mains\PAT\PATAQParameters.mat','PATAQInfo');

%%%% RF data save name
RFdataFilename_pa = ['RF-',P(2).Filename,'-1'];
Filename_pa = P(2).Filename;

%% 2. System parameters
filename = mfilename; % used to launch VSX automatically
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.connector = 1;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = simulateMode;

Resource.VDAS.halDebugLevel = 1;

P(1).nCh=Resource.Parameters.numRcvChannels ;
P(2).nCh=Resource.Parameters.numRcvChannels ;
P(1).vSound=Resource.Parameters.speedOfSound;
P(2).vSound=Resource.Parameters.speedOfSound;
%% 3. Transducer parameters
Trans.name = 'L22-14v';
Trans.units = 'wavelengths';
% Trans.frequency = 15;   % The center frequency for the A/D 4xFc sampling.
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 30; % mfr data sheet lists 30 Volt limit
P(1).LensDelay=Trans.lensCorrection;
P(1).frequency=Trans.frequency;
P(1).pitch=Trans.spacingMm;
P(2).LensDelay=Trans.lensCorrection;
P(2).frequency=Trans.frequency;
P(2).pitch=Trans.spacingMm;
%% 4. PData structure array
PData(1).PDelta = [Trans.spacing/2, 0, 0.5];% [0.4, 0, 0.25];
PData(1).Size(1) = ceil((P(1).EndDepth-P(1).StartDepth+1)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Resource.Parameters.numTransmit*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Resource.Parameters.numTransmit-1)/2,0,P(1).StartDepth]; % x,y,z of upper lft crnr.

PData(2).PDelta = [Trans.spacing/2, 0, 0.5];%[0.4, 0, 0.25];
PData(2).Size(1) = ceil((P(2).EndDepth-P(2).StartDepth+1)/PData(2).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(2).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1));
PData(2).Size(3) = 1;      % single image page
PData(2).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P(2).StartDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

%% 	5. Media object (for simulation)
pt1; % Specify Media object. 'pt1.m' script defines array of point targets.
Media.attenuation = -0.5;
Media.function = 'movePoints';
% Media.MP(1,:) = [0,0,70,1]; % Single point, [x, y, z, reflectivity]
%% 	6. Maximum acquisition length and number of z samples
maxAcqLength = ceil(sqrt(P(1).EndDepth^2 + ((Resource.Parameters.numRcvChannels-1)*Trans.spacing)^2));
actAcqLength =(maxAcqLength-P(1).StartDepth); 
P(1).actZsamples=ceil(actAcqLength*2*P(1).nSperWave/128)*128; % *2 for round trip; P.nSperWave: sampling points per wavelength, rounded up to the next 128 VSX sample boundary

maxAcqLength_pa = ceil(sqrt(P(2).EndDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
actAcqLength_pa =(maxAcqLength_pa-P(2).StartDepth/2);
P(2).actZsamples=ceil(actAcqLength_pa * 2 * P(2).nSperWave/128)*128; % P.nSperWave: sampling points per wavelength, rounded up to the next 128 VSX sample boundary

DataSize=P(1).actZsamples*Resource.Parameters.numRcvChannels*P(1).numAcqs*2/1024^2; % MBytes  ?????
allDataSize = (P(1).actZsamples*Resource.Parameters.numRcvChannels*P(1).numAcqs*2+P(2).actZsamples*Resource.Parameters.numRcvChannels*P(2).numAcqs*2)/1024^2;
tDMA=DataSize/(6.6*1024)*1e6+2e3;% in us; DMA rate: 6.6 Gbytes/s; 2 ms overhead time
P(1).tMiniIntSupFrame=ceil(tDMA/power(10,floor(log10(tDMA))))*power(10,floor(log10(tDMA)))+10e3; % in us; time interval between successive super frame
P(1).tIntPDI=1e6/P(1).PDIFR;
if P(1).tIntPDI<P(1).tMiniIntSupFrame
    P(1).tIntPDI=P(1).tMiniIntSupFrame;
    P(1).PDIFR=1E6/P(1).tIntPDI;
    prompt={'Greater than maximum frame rate! Maximu FR: '};
    name='File info';
    defaultvalue={num2str(P(1).PDIFR)};
    inputValue=inputdlg(prompt,name, 1, defaultvalue);
    P(1).PDIFR=str2num(inputValue{1});
    P(1).tIntPDI=1e6/P.PDIFR;
end
    
%% 	7. Resource parameters
%%%% ? Receive buffer
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P(1).actZsamples*P(1).numAcqs + P(2).actZsamples*P(2).numAcqs;   % this size allows for maximum range in wavelength
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P(1).numSupFrames;       % number of 'super frames'
%%%% ? IQ buffer and Image buffer
Resource.InterBuffer(1).numSupFrames = 1;  % only one intermediate buffer needed.
Resource.ImageBuffer(1).numSupFrames = P(1).numSupFrames; % for reduced online visualization
%%%%%% Resource.ImageBuffer(1).numSupFrames = P.numSupFrames*P.numAcqs; % for all acquired data reconstruction
for iwl = 1: P(2).numSupFrames
Resource.InterBuffer(iwl+1).numSupFrames = 1;  % only one intermediate buffer needed.
Resource.ImageBuffer(iwl+1).numSupFrames = 1; % for reduced online visualization
end

%%%% ? Display window parameters
Resource.DisplayWindow(1).Title = 'L22-14v-CCHFR';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
    DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numSupFrames = 1;% ? 20
%%%%%% Resource.DisplayWindow(1).numSupFrames = P.numSupFrames*P.numAcqs; % for displaying all acquired images
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).splitPalette = 1;

%% 	8. Transmit waveform definition
TW(1).type = 'parametric';
TW(1).Parameters = [P.TWfrequency, 0.67, P.TWnHC, 1];   % 18 MHz center frequency, 67% pulsewidth 1.5 cycle burst
[Wvfm1Wy64Fc, peak, numsamples, rc, TWout] = computeTWWaveform(TW);
P(1).PeakDelay=peak;

TW(2).type = 'parametric';
TW(2).Parameters = [18, 0.67, 2, 1];   % 18 MHz center frequency, 67% pulsewidth two cycle burst
%% 	9. TX structures
%%%% ? Set TX structure array for P.numAngles transmits
% emitElem=ones(1,Trans.numelements);
emitElem=kaiser(Resource.Parameters.numTransmit,1)';
% emitElem=zeros(1,Trans.numelements);
% emitElem(15:114)=1;

TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', emitElem, ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, P(1).numAngles+1);    
%%%% ? Set specific TX event for the P.numAngles transmits
if fix(P(1).numAngles/2) == P(1).numAngles/2       % if P.numAngles even
    P(1).startAngle = (-(fix(P(1).numAngles/2) - 1) - 0.5)*dtheta;
else
    P(1).startAngle = -fix(P(1).numAngles/2)*dtheta;
end
for iAngle = 1:P(1).numAngles   % P.numAngles transmit events
    TX(iAngle).Steer = [(P(1).startAngle+(iAngle-1)*dtheta),0.0];
    TX(iAngle).Delay = computeTXDelays(TX(iAngle));
    P(1).TXDelay(iAngle,:)=TX(iAngle).Delay;
end
P(1).dAngle=dtheta;

%-- only one TX struct needed for PA
TX(P(1).numAngles+1).waveform = 2;            % TW(2)
TX(P(1).numAngles+1).Steer = [0.0,0.0];       % theta, alpha = 0.
if transmit
    TX(P(1).numAngles+1).Apod = ones(1,Trans.numelements);    % This is the conventional T/R condition and invokes the default beamformer
else
  TX(P(1).numAngles+1).Apod = zeros(1,Trans.numelements);   % THIS COMMAND TURNS OFF ALL TRANSMITTERS AND INVOKES THE RECEIVE-ONLY BEAMFORMER 
end
TX(P(1).numAngles+1).Delay = computeTXDelays(TX(P(1).numAngles+1));

%% 10. Analog front end gain settings.
RcvProfile(1).LnaGain = 18;     % 12, 18, or 24 dB  (18=default)
RcvProfile(1).condition = 'immediate';

RcvProfile(2).LnaGain = 24;
RcvProfile(2).condition = 'immediate';

%% 11. TGC Waveform structure.
% TGC.CntrlPts = [650,720,800,860,910,980,1000,1000];
TGC(1).CntrlPts = [800,880,940,980,1000,1000,1000,1000];
TGC(1).rangeMax = P(1).EndDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));

TGC(2).CntrlPts = [500,590,650,710,770,830,890,950];    % Standard
% TGC(1).CntrlPts = round(950.^linspace(0,1,8));    % Exponential
TGC(2).rangeMax = P(2).EndDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));


%% 11. Receive structures
%%%% ? Set RCV structure array for all acquisition
BPF1 = [ -0.00009 -0.00128 +0.00104 +0.00085 +0.00159 +0.00244 -0.00955 ...
         +0.00079 -0.00476 +0.01108 +0.02103 -0.01892 +0.00281 -0.05206 ...
         +0.01358 +0.06165 +0.00735 +0.09698 -0.27612 -0.10144 +0.48608 ];
          
rcvElem=ones(1,Resource.Parameters.numTransmit);
% rcvElem=zeros(1,Trans.numelements);
% rcvElem(45:84)=1;
Receive = repmat(struct('Apod',  rcvElem, ...
                 'startDepth', P(1).StartDepth, ...
                 'endDepth', P(1).StartDepth+actAcqLength,...
                 'TGC', 1, ...
                 'bufnum', 1, ...
                 'framenum', 1, ...
                 'acqNum', 1, ...
                 'sampleMode', P(1).SampleMode, ...
                 'InputFilter', BPF1, ...
                 'mode', 0, ...
                 'callMediaFunc', 0), 1, P(1).numAcqs*P(1).numSupFrames+P(2).numAcqs);
%%%% ? Set specific RCV event for all receive acquisitions
for iFrame = 1:Resource.RcvBuffer(1).numFrames
    Receive(P(1).numAcqs*(iFrame-1) + 1).callMediaFunc = 1;  % move points only once per super frame
    for iAcq = 1:P(1).numAcqs
        % -- Acquisitions for 'super' frame.
        rcvNum = P(1).numAcqs*(iFrame-1) + iAcq;
        Receive(rcvNum).Apod(:)=1;
%         Receive(rcvNum).callMediaFunc = 1;  % movepoints EVERY acquisition to illustrate superframe concept
        Receive(rcvNum).framenum = iFrame;
        Receive(rcvNum).acqNum = iAcq;
    end
end             
P(1).maxDepth=Receive(1).endDepth;    % the maximum data acquisition depth

for iFrame_pa = 1: P(2).numSupFrames
    Receive(rcvNum+(iFrame_pa-1)*P(2).AverageNumber + 1).callMediaFunc = 1; % move points every wavelegth in PA in simulation;
    for iAcq_pa  = 1: P(2).AverageNumber
    rcvNum_pa = rcvNum + (iFrame_pa-1)*P(2).AverageNumber + iAcq_pa;
    Receive(rcvNum_pa).startDepth = P(2).StartDepth/2;
    Receive(rcvNum_pa).endDepth = maxAcqLength_pa;
    Receive(rcvNum_pa).TGC = 2;
    Receive(rcvNum_pa).sampleMode = P(2).SampleMode;
    Receive(rcvNum_pa).acqNum = rcvNum_pa;
    Receive(rcvNum_pa).InputFilter = []; % no filter for PA
    end
end   
P(2).maxDepth = Receive(rcvNum_pa).endDepth;

%% 12. Recon structure arrays
%%%% Reconstuction for a super Frame
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'newFrameTimeout', -1,...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', zeros(1,1)),1,P(1).numSupFrames + P(2).numSupFrames);
for iFrame = 1: P(1).numSupFrames
    Recon(iFrame).newFrameTimeout = P(1).tIntPDI*1.2/1000;  % ????
    Recon(iFrame).RINums = 1:P(1).numAngles;
end
for iFrame_pa = 1: P(2).numSupFrames
    k = iFrame_pa+iFrame;
    Recon(k).IntBufDest = [k,1];
    Recon(k).ImgBufDest = [k,-1];
    Recon(k).pdatanum = 2;
    Recon(k).RINums = (P(1).numAngles+1+P(2).AverageNumber*(iFrame_pa-1)) : (P(1).numAngles+P(2).AverageNumber*iFrame_pa);
end
%%%% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, P(1).numAngles+P(2).numAcqs);
%%%% - Set specific ReconInfo attributes.
ReconInfo(1).mode = 'replaceIQ';
for iAngle = 1:P(1).numAngles
    ReconInfo(iAngle).txnum = iAngle;
    ReconInfo(iAngle).rcvnum = iAngle;
end
ReconInfo(P(1).numAngles).mode = 'accumIQ_replaceIntensity';

for iFrame_pa = 1: P(2).numAcqs
    ReconInfo(iFrame_pa+iAngle).txnum = P(1).numAngles + 1;
    ReconInfo(iFrame_pa+iAngle).rcvnum = P(1).numAcqs + iFrame_pa;
end
for iwl = 1: P(2).numSupFrames
    ReconInfo(P(1).numAngles+1+(iwl-1)*P(2).AverageNumber).mode = 'replaceIQ';
if P(2).AverageNumber == 1
    ReconInfo(P(1).numAngles+P(2).AverageNumber*iwl).mode = 'replaceIntensity';
else
    ReconInfo(P(1).numAngles+P(2).AverageNumber*iwl).mode = 'accumIQ_replaceIntensity';
end
end

%% 13. Process structure array
cpt = 22;       % define here so we can use in UIControl below
cpers = 80;     % define here so we can use in UIControl below
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

for k = 1: P(2).numSupFrames
    Process(k+1) = Process(1);
    Process(k+1).Parameters(2) = {k+1};
    Process(k+1).Parameters(6) = {2};
    Process(k+1).Parameters(14) = {cpers};
    %Process(k+1).Parameters.threshold = cpt;
end
    
%% 14. SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = P(1).t2NextPlaneDAQ;  % in usec, time between successive angles 
SeqControl(3).command = 'timeToNextAcq';  % time between successive CCUS images,
SeqControl(3).argument = tCCUS-SeqControl(2).argument*(P(1).numAngles-1);  % time between successive CCUS images, CCUS frame rate=1/tCCUS 

SeqControl(4).command = 'triggerOut';
SeqControl(5).command = 'returnToMatlab';
SeqControl(6).command = 'timeToNextAcq';  % time between super frames
SeqControl(6).argument = P(1).tIntPDI-tCCUS*P(1).numCCframes;  % in us; time for data trasfer to host, amd to next superframe/PDI acquisition
% set receive profile
SeqControl(7).command = 'setRcvProfile';
SeqControl(7).argument = 1;
SeqControl(8).command = 'setRcvProfile';
SeqControl(8).argument = 2;
% input trigger
SeqControl(9).command = 'triggerIn';
SeqControl(9).condition = 'Trigger_2_Rising'; % Trigger input 2, enable with rising edge
SeqControl(9).argument = 100; % 25 sec timeout delay
% (Timeout range is 1:255 in 250 msec steps; 0 means timeout disabled)
% noop delay between trigger in and start of acquisition
SeqControl(10).command = 'noop';
SeqControl(10).argument = round(flash2Qdelay*5); % noop counts are in 0.2 microsec increments
% sync command
SeqControl(11).command = 'sync';
SeqControl(11).argument = 10e6; % Wait 30 seconds
% -- PRF for PA ensemble
SeqControl(12).command = 'timeToNextAcq';
SeqControl(12).argument = round(1/(PA_PRF*1e-06)); % (50 msecs for PA_PRF=20 Hz)



nsc = 13; % nsc is count of SeqControl objects

%% 15. Event objects to acquire all acquisitions
n = 1; % n is count of Events
for iSuperFrame = 1:Resource.RcvBuffer(1).numFrames
         % Acquire PA ensemble.
    for iFrame_pa = (P(1).numAcqs+1):(P(1).numAcqs+P(2).numAcqs)
        % Wait for input trigger from flash lamp firing
        Event(n).info = 'Wait for Trigger IN';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 9; %wait for trigger in, timeout
        n = n+1;

        % Pause for optical buildup
        Event(n).info = 'noop';
        Event(n).tx = 0;
        Event(n).rcv = 0;       
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [10,11];
        n = n+1;

        % send trigger output at start of every PA acquisition to fire Q-switch
        Event(n).info = 'Acquire PA event';
        Event(n).tx = P(1).numAngles+1;
        Event(n).rcv = (P(1).numAcqs+P(2).numAcqs)*(iSuperFrame-1)+iFrame_pa;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 4; % send trigger to laser Q switch
        n = n+1;
        
%         Event(n).info = 'Custom Sw & Hw sync timeout.';   % sync timeout is .5s by default
%         Event(n).tx = 0;            % no TX structure.
%         Event(n).rcv = 0;           % no Rcv structure.
%         Event(n).recon = 0;         % no reconstruction.
%         Event(n).process = 0;       % no processing
%         Event(n).seqControl = 11;    % make the Sw wait at most 60s for the Hw to sync
%         n = n + 1;
    end
     Event(n-1).seqControl = [4,7]; % 12% replace last PA acquisition Event's seqControl with longer TTNA and RCV profile change
     
    for iCCframe = 1:P(1).numCCframes
        for iAngle=1:P(1).numAngles
            Event(n).info = 'Acquire RF';
            Event(n).tx = iAngle;
            Event(n).rcv = P(1).numAcqs*(iSuperFrame-1) + P(1).numAngles*(iCCframe-1)+iAngle;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = [2];
            n = n+1;
        end
        Event(n-1).seqControl = [3,4];
    end
    Event(n-1).seqControl = [6,4,8]; % time to superframe  Set last US acquisition Event's seqControl (longer TTNA and new RCV profile)
    
    Event(n).info = 'Transfer Data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc; %%%%???????? additional time to transfer data?
    n = n+1;
    SeqControl(nsc).command = 'transferToHost'; % tDMA transfer frame to host buffer (each needs a different value of nsc)
      nsc = nsc+1;

%     Event(n).info = 'recons and 2D process';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = [1,2:P(2).numSupFrames+1]; % recostruct 1 US frame + n(#of wls) PA frames 
%     Event(n).process = 1;
%     Event(n).seqControl = 0;
%     n = n+1;
% 
%     Event(n).info = 'PA image display';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 0;
%     Event(n).process = 2:P(2).numSupFrames+1; 
%     Event(n).seqControl = 0;
% %     if floor(i/frameRateFactor) == i/frameRateFactor     % Exit to Matlab only every 3rd frame to prevent slowdown
% %         Event(n).seqControl = nsc;
% %         SeqControl(nsc).command = 'returnToMatlab';
% %         nsc = nsc+1;
% %     end
%     n = n+1;
end
Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;


%$ 16. User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
% - Range Change
MinMaxVal = [64,300,P(1).EndDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
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
return
% commandwindow  % just makes the Command window active to show printout
% P(1).nSmplPerWvlnth=Receive(1).samplesPerWave;
%% 18. after quitting VSX, save the RFData collected at high frame rate 
% if simulateMode ~= 2    
%     disp ('Info:  Saving the RF Data buffer -- please wait!'), disp(' ')
%     if exist([savepath,RFdataFilename,'.mat'])
%         NewPlane=1;
%         while (exist([savepath,RFdataFilename,'.mat'])==2)
%             NewPlane=NewPlane+1;
%             RFdataFilename = ['RF-',num2str(P.CCangle),'-',num2str(P.numAngles),'-',num2str(P.CCFR),'-',num2str(P.numCCframes),'-',num2str(P.numSupFrames),'-',Filename,'-',num2str(NewPlane)];
%         end
%         prompt={'File exist! New Plane Sequence: '};
%         name='File info';
%         defaultvalue={num2str(NewPlane)};
%         inputValue=inputdlg(prompt,name, 1, defaultvalue);
%         P.Plane=str2num(inputValue{1});
%         RFdataFilename = ['RF-',num2str(P.CCangle),'-',num2str(P.numAngles),'-',num2str(P.CCFR),'-',num2str(P.numCCframes),'-',num2str(P.numSupFrames),'-',Filename,'-',num2str(P.Plane)];
% %         IQdataFilename = ['IQ-',num2str(P.CCangle),'-',num2str(P.numAngles),'-',num2str(P.CCFR),'-',num2str(P.numCCframes),'-',num2str(P.numSupFrames),'-',Filename,'-',num2str(P.Plane)]; 
%     end
%     tic
%     RFRAW=RcvData{1};IQ=IQData{1};
%     savefast ([savepath,RFdataFilename], 'RFRAW','P')
% %     savefast ([savepath,IQdataFilename], 'IQ','P')
%     toc
%     disp ('RF DATA SAVED!')
% end
% return

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
    set(hObject,'Value',evalin('base','EndDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
EndDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        EndDepth = UIValue*scaleToWvl;    
    end
end
assignin('base','P',P);

evalin('base','PData(1).Size(1) = ceil((P(1).EndDepth-P(1).StartDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(EndDepth^2 + ((Resource.Parameters.numRcvChannels-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = EndDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback
