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
% 11-1-2019 by Bingxue Liu

clear all
cd C:\Users\BOAS-US\Documents\Vantage-3.4.0-1711281030
addpath('D:\CODE\Functions');
load('D:\CODE\Mains\DAQParameters_L12_3v.mat');
deftpath=DAQInfo.savepath;
savepath=[uigetdir(deftpath),'\'];

%% 1. DAQ parameter 
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
    'Sampling Mode (100, 200 %)'};
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
    '6', '3','100'};
inputValue=inputdlg(prompt,name, 1, defaultvalue);
Filename=inputValue{1};
P.PDIFR=str2num(inputValue{2});
P.numSupFrames=str2num(inputValue{3});
P.CCFR=str2num(inputValue{4});
P.numCCframes=str2num(inputValue{5});
P.CCangle=str2num(inputValue{6});
P.numAngles=str2num(inputValue{7});
StartDepth=str2num(inputValue{8});
EndDepth=str2num(inputValue{9});
P.TWfrequency=str2num(inputValue{10});
P.TWnHC=str2num(inputValue{11});
P.SplMode=str2num(inputValue{12});
if P.SplMode==200
    P.SampleMode='NS200BW';
    P.nSperWave=4; % sampling points per wavelength
else
    P.SampleMode='BS100BW';
    P.nSperWave=2; % sampling points per wavelength
end
%% time-to-next acquisition
P.t2NextPlaneDAQ=33;    % us, time-to-next acquisition between successive plane wave, travel time + system overhead time, the default is 33 us for 150 wavelengths acquisition
t2NextCCPlaneDAQ=round(1e6/P.CCFR); %us, time-to-next acquisition between successive CCUS frame,us 
P.tIntPDI=round(1e6/P.PDIFR);
P.t2NextSupPlaneDAQ=1e6/P.PDIFR; %us, time-to-next acquisition between super frames
if P.t2NextSupPlaneDAQ-t2NextCCPlaneDAQ*P.numCCframes<0 
    P.t2NextSupPlaneDAQ=t2NextCCPlaneDAQ*(P.numCCframes)+2.5e6;
    P.PDIFR=1e6/P.t2NextSupPlaneDAQ;
    prompt={['Incorrect superFrame rate, Approriate Range: ', num2str(1e6/(t2NextCCPlaneDAQ*P.numCCframes+4e6)) ,'-',num2str(1e6/(t2NextCCPlaneDAQ*(P.numCCframes)+2.5e6))]};
    name='Wrong superFrame rate';
    defaultvalue={num2str(P.PDIFR)};
    inputValue=inputdlg(prompt,name, 1, defaultvalue);
    P.PDIFR=str2num(inputValue{1});
    P.t2NextSupPlaneDAQ=1e6/P.PDIFR; 
elseif (P.t2NextSupPlaneDAQ-t2NextCCPlaneDAQ*P.numCCframes)>4e6
    P.t2NextSupPlaneDAQ=4e6;
end

DAQInfo=P;
DAQInfo.StartDepth=StartDepth;
DAQInfo.EndDepth=EndDepth;
DAQInfo.Filename=Filename;
DAQInfo.savepath=savepath;
save('D:\CODE\Mains\DAQParameters.mat_L12_3v','DAQInfo');

%%%% RF data save name
RFdataFilename = ['RF-',num2str(P.CCangle),'-',num2str(P.numAngles),'-',num2str(P.CCFR),'-',num2str(P.numCCframes),'-',num2str(P.numSupFrames),'-',Filename,'-', '1'];
fileID = fopen([savepath,'DAQ-INFO.txt'],'at');
fmt = ['Data acqusiition information.\n\n'];
fprintf(fileID,fmt);
fclose(fileID);
%% %%%%%%%%%% select working mode %%%%%%%%%%%%%%%%%%%%
simulateMode = 0;  % acquire data with hardware
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
%% 2. System parameters
filename = mfilename; % used to launch VSX automatically
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
% Resource.Parameters.numTransmit = 80;      % number of transmit channels.
% Resource.Parameters.numRcvChannels = 80;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.connector = 2;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = simulateMode;

P.nCh=Resource.Parameters.numRcvChannels ;
P.vSound=Resource.Parameters.speedOfSound;
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
PData(1).PDelta = [0.4, 0, 0.25];
PData(1).Size(1) = ceil((EndDepth-StartDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Resource.Parameters.numTransmit*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Resource.Parameters.numTransmit-1)/2,0,StartDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

%% 	5. Media object (for simulation)
pt1; % Specify Media object. 'pt1.m' script defines array of point targets.
Media.attenuation = -0.5;
Media.function = 'movePoints';
% Media.MP(1,:) = [0,0,70,1]; % Single point, [x, y, z, reflectivity]
%% 	6. Maximum acquisition length and number of z samples
maxAcqLength = ceil(sqrt(EndDepth^2 + ((Resource.Parameters.numRcvChannels-1)*Trans.spacing)^2));
actAcqLength =(maxAcqLength-StartDepth); 
P.actZsamples=ceil(actAcqLength *2*P.nSperWave/128)*128; % *2 for round trip; P.nSperWave: sampling points per wavelength, rounded up to the next 128 VSX sample boundary

DataSize=P.actZsamples*Resource.Parameters.numRcvChannels*P.numAngles*P.numCCframes*2/1024^2; % MBytes
tDMA=DataSize/(6.6*1024)*1e6+2e3;% in us; DMA rate: 6.6 Gbytes/s; 2 ms overhead time
P.tMiniIntSupFrame=ceil(tDMA/power(10,floor(log10(tDMA))))*power(10,floor(log10(tDMA)))+10e3; % in us; time interval between successive super frame
P.t2NextSupPlaneDAQ=1e6/P.PDIFR;
if P.t2NextSupPlaneDAQ<P.tMiniIntSupFrame
    P.t2NextSupPlaneDAQ=P.tMiniIntSupFrame;
    P.PDIFR=1E6/P.t2NextSupPlaneDAQ;
    prompt={'Greater than maximum frame rate! Maximu FR: '};
    name='File info';
    defaultvalue={num2str(P.PDIFR)};
    inputValue=inputdlg(prompt,name, 1, defaultvalue);
    P.PDIFR=str2num(inputValue{1});
    P.t2NextSupPlaneDAQ=1e6/P.PDIFR;
end
    
%% 	7. Resource parameters
%%%% ? Receive buffer
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.actZsamples*P.numAcqs;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numSupFrames;       % number of 'super frames'
%%%% ? IQ buffer and Image buffer
Resource.InterBuffer(1).numSupFrames = 1;  % only one intermediate buffer needed.
Resource.ImageBuffer(1).numSupFrames = P.numSupFrames; % for reduced online visualization
%%%%%% Resource.ImageBuffer(1).numSupFrames = P.numSupFrames*P.numAcqs; % for all acquired data reconstruction
%%%% ? Display window parameters
Resource.DisplayWindow(1).Title = 'L12-3v-CCHFR';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
    DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numSupFrames = 20;
%%%%%% Resource.DisplayWindow(1).numSupFrames = P.numSupFrames*P.numAcqs; % for displaying all acquired images
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Resource.VDAS.halDebugLevel = 1;

%% 	8. Transmit waveform definition
TW(1).type = 'parametric';
TW.Parameters = [P.TWfrequency, 0.67, P.TWnHC, 1];   % 18 MHz center frequency, 67% pulsewidth 1.5 cycle burst
[Wvfm1Wy64Fc, peak, numsamples, rc, TWout] = computeTWWaveform(TW);
P.PeakDelay=peak;
%% 	9. TX structures
%%%% ? Set TX structure array for P.numAngles transmits
% emitElem=ones(1,Trans.numelements);
emitElem=kaiser(Resource.Parameters.numTransmit,1)';
% emitElem=zeros(1,Trans.numelements);
% emitElem(15:114)=1;

TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'aperture', 33, ...
                   'Apod', emitElem, ...
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
    P.TXDelay(iAngle,:)=TX(iAngle).Delay;
end
P.dAngle=dtheta;
%% 10. TGC Waveform structure.
% TGC.CntrlPts = [650,720,800,860,910,980,1000,1000];
TGC.CntrlPts = [750,820,880,910,970,980,1000,1000];
TGC.rangeMax = EndDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% 11. Receive structures
%%%% ? Set RCV structure array for all acquisition
BPF1 = [ -0.00009 -0.00128 +0.00104 +0.00085 +0.00159 +0.00244 -0.00955 ...
         +0.00079 -0.00476 +0.01108 +0.02103 -0.01892 +0.00281 -0.05206 ...
         +0.01358 +0.06165 +0.00735 +0.09698 -0.27612 -0.10144 +0.48608 ];

Receive = repmat(struct('Apod', ones(1,Resource.Parameters.numTransmit), ...
                        'aperture', 33, ...
                        'startDepth', StartDepth, ...
                        'endDepth', StartDepth+actAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', P.SampleMode, ...
                        'InputFilter', BPF1, ... 
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, P.numAcqs*P.numSupFrames);


P.nSmplPerWvlnth=P.nSperWave;
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
               
P.startDepth=Receive(1).startDepth;% the defined image start depth
P.endDepth=EndDepth;               % the defined image end depth
P.maxDepth=Receive(1).endDepth;    % the maximum data acquisition depth

%% 12. Recon structure arrays
%%%% Reconstuction for a super Frame
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'newFrameTimeout', P.t2NextSupPlaneDAQ*1.2/1000,...
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
% Define the external processing event for RF data saving. 
Process(2).classname = 'External';
Process(2).method = 'SaveRFonline';
Process(2).Parameters = {'srcbuffer','receive',... % name of buffer to process.
    'srcbufnum',1,... 
    'srcframenum',1,... 
    'dstbuffer','none'};

%% 14. SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = P.t2NextPlaneDAQ;  % in usec, time between successive angles 
SeqControl(3).command = 'timeToNextAcq';  % time between successive CCUS images,
SeqControl(3).argument = t2NextCCPlaneDAQ-SeqControl(2).argument*(P.numAngles-1);  % time between successive CCUS images, CCUS frame rate=1/t2NextCCPlaneDAQ 

SeqControl(4).command = 'triggerOut';
SeqControl(5).command = 'returnToMatlab';
SeqControl(6).command = 'timeToNextAcq';  % time between super frames
SeqControl(6).argument = min(P.t2NextSupPlaneDAQ-t2NextCCPlaneDAQ*P.numCCframes,4e6);  % in us; time for data trasfer to host, amd to next superframe/PDI acquisition
nsc = 7; % nsc is count of SeqControl objects

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
            Event(n).seqControl = [2,4];
            n = n+1;
        end
            Event(n-1).seqControl = [3,4];
    end
    % Set last acquisitions SeqControl for transferToHost.
    Event(n-1).seqControl = [6,4,nsc];
        SeqControl(nsc).command = 'transferToHost'; % transfer all acqs in one super frame
        nsc = nsc + 1;
    %
    Event(n).info = 'Call external Processing function.'; 
    Event(n).tx = 0;       % no TX structure.
    Event(n).rcv = 0;      % no Rcv structure.
    Event(n).recon = 0;    % no reconstruction.
    Event(n).process = 2;  % call SaveRFonline
%     Event(n).seqControl = [5];
    n=n+1;
%     % Do reconstruction and processing for 1st sub frame
%     Event(n).info = 'Reconstruct'; 
%     Event(n).tx = 0;         
%     Event(n).rcv = 0;        
%     Event(n).recon = 1;      
%     Event(n).process = 1;    
%     Event(n).seqControl = 5;
%     n = n+1;
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
MinMaxVal = [64,300,EndDepth]; % default unit is wavelength
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
save('D:\CODE\Mains\P.mat','P');
disp([ mfilename ': NOTE -- Running VSX automatically!']), disp(' ')
VSX    
% commandwindow  % just makes the Command window active to show printout
P.nSmplPerWvlnth=Receive(1).samplesPerWave;
%% 18. after quitting VSX, save the RFData collected at high frame rate 
% if simulateMode ~= 2    
%     disp ('Info:  Saving the RF Data buffer -- please wait!'), disp(' ')
%     if exist([savepath,RFdataFilename,'.mat'])
%         NewPlane=1;
%         while (exist([savepath,RFdataFilename,'.mat'])==2)
%             NewPlane=NewPlane+1;
%             RFdataFilename = [num2str(P.CCangle),'-',num2str(P.numAngles),'-',num2str(P.CCFR),'-',num2str(P.numCCframes),'-',num2str(P.numSupFrames),'-',Filename,'-',num2str(NewPlane),'-RF'];
%         end
%         prompt={'File exist! New Plane Sequence: '};
%         name='File info';
%         defaultvalue={num2str(NewPlane)};
%         inputValue=inputdlg(prompt,name, 1, defaultvalue);
%         P.Plane=str2num(inputValue{1});
%         RFdataFilename = [num2str(P.CCangle),'-',num2str(P.numAngles),'-',num2str(P.CCFR),'-',num2str(P.numCCframes),'-',num2str(P.numSupFrames),'-',Filename,'-',num2str(P.Plane),'-RF'];
%     end
%     tic
%     RFRAW=RcvData{1};
%     savefast ([savepath,RFdataFilename], 'RFRAW','P')
%     toc
%     disp ('RF DATA SAVED!')
% end
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

evalin('base','PData(1).Size(1) = ceil((EndDepth-StartDepth)/PData(1).PDelta(3));');
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
