% File name: RunSetUpL122_14VPAT_Acquire_US_MASTER.m
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
% 20201016 by Bingxue Liu

clear all

QSwitch_delay_us = 300; % Should be >= 300 for Laser External trigger + QSwitch external trigger

if QSwitch_delay_us < 300
    fprintf("QSwitch_delay_us = %d which is less than 300!\nThis would damage the laser. ABORTING!\n", QSwitch_delay_us);
    return
end

trigger = true;
transmit = false;

%%%%
cd C:\Users\BOAS-US\Desktop\Vantage-3.4.3-1807021300 %C:\Users\BOAS-US\Documents\Vantage-3.4.0-1711281030
% addpath('D:\CODE\Functions');
load('D:\CODE\Mains\PAT\PATAQParameters_L22_14v.mat');
deftpath=PATAQInfo.savepath;
simulateMode = 2;
addpath('D:\CODE\Functions');
%% 1. Load RFdata and parameter
if simulateMode ==2
    [FileName,FilePath]=uigetfile(deftpath);
    disp('Loading Data...')
    load ([FilePath, FileName]);  % this file includes RcvData and the numAcqs and numFrames parameters for properly setting Resources
    PATAQInfo.savepath=FilePath;
    FileInfo=strsplit(FileName,'-');
    PATAQInfo.Filename=FileInfo{4};
    PATAQInfo.FileNameFull=FileName;
    save('D:\CODE\Mains\PATAQParameters.mat','PATAQInfo','P');
    RcvData{1}=RFRAW;
    clear RFRAW RFRAW0;
    clear Img Img0;
    P.simulateMode = simulateMode;
        %% 1.1. IQ beamforming parameter
    prompt={'PDelta-X',...
        'PDelta-Z'};
    name='Beamforming parameters';
    defaultvalue={'0.5', '0.5'};
    inputValue=inputdlg(prompt,name, 1, defaultvalue);
    PDeltaX=str2num(inputValue{1});
    PDeltaZ=str2num(inputValue{2});
%else
end   
disp('Data Loaded!')

IQdataFilename = ['IQ-',FileName(4:end)];
IMGdataFilename = ['Img-',FileName(4:end)];
%% 2. System parameters
filename = mfilename; % used to launch VSX automatically
Resource.Parameters.numTransmit = 128;     % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels. CHANGE HERE
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.connector = 1;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = simulateMode;

%Resource.VDAS.halDebugLevel = 1;

P.nCh=Resource.Parameters.numRcvChannels ;
P.vSound=Resource.Parameters.speedOfSound;
% RcvProfile.condition = 'immediate';
%% 3. Transducer parameters
Trans.name = 'L22-14v';
Trans.units = 'wavelengths';
Trans.connType = 0;
% Trans.frequency = 15;   % The center frequency for the A/D 4xFc sampling.
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 30; % mfr data sheet lists 30 Volt limit
P.LensDelay=Trans.lensCorrection;
P.frequency=Trans.frequency;
P.pitch=Trans.spacingMm;
P.wavelength=P.vSound/P.frequency*1e-3;% trans wavelength, mm
% P.EndDepth = P.endDepth;
% P.StartDepth = P.startDepth;
%% 4. PData structure array

PData.PDelta = [Trans.spacing*PDeltaX, 0, PDeltaZ];
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
Resource.InterBuffer(1).numFrames = P.AverageNumber;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = P.AverageNumber;

% Display
Resource.DisplayWindow(1).Title = 'L22-14v PAT Reconstruction';    %%%%%, 4X sampling at 62.5 MHz';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Matlab';
Resource.DisplayWindow(1).numFrames = P.AverageNumber;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

%% 	8. Transmit waveform definition
TW.type = 'parametric';
TW.Parameters = [18, 0.67, 2, 1];   % 18 MHz center frequency, 67% pulsewidth two cycle burst

%% 	9. TX structures
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(1).focus = 0;% TX(1).aperture = 0; 
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
if transmit
    TX(1).Apod = ones(1,Trans.numelements);
else
    TX(1).Apod = zeros(1,Trans.numelements);
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
%  partial = [zeros(1,64),ones(1,64)];    
Receive = repmat(struct(   'Apod', ones(1,Resource.Parameters.numRcvChannels), ... %   CHANGE HERE ones(1,Resource.Parameters.numRcvChannels), ...%'aperture', 1, ...  % CHANGE HERE 33
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
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
...%                'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1 ), 1, P.AverageNumber);           
               %'newFrameTimeout', -1)
              

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...  
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, P.AverageNumber);

% ReconInfo(1).mode = 'replaceIQ';
% 
for iAcq=1:P.AverageNumber
    Recon(iAcq).IntBufDest = [1, iAcq];
    Recon(iAcq).RINums = iAcq;
    ReconInfo(iAcq).rcvnum = iAcq;
end
% 
% if P.AverageNumber == 1
%     ReconInfo(P.AverageNumber).mode = 'replaceIntensity';
% else
%     ReconInfo(P.AverageNumber).mode = 'accumIQ_replaceIntensity';
% end

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


                     
%% 15. SeqControl structure arrays.
num_evt = 1;
start_evt = num_evt;

SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;

% SeqControl(2).command = 'triggerIn';  % time between frames
% SeqControl(2).condition = 'Trigger_2_Rising';  % Wait for trigger 2 rising edge
% if trigger
%     SeqControl(2).command = 'triggerIn';  % time between frames
%     SeqControl(2).condition = 'Trigger_2_Rising';  % Wait for trigger 1 rising edge
% %     SeqControl(2).argument = 0;
% else
%     SeqControl(2).command = 'timeToNextAcq';  % time between frames
%     SeqControl(2).argument = 50000;  % 50 msec
% end

SeqControl(2).command = 'triggerOut';

% SeqControl(3).command = 'transferToHost';
% 
% SeqControl(4).command = 'returnToMatlab';

% SeqControl(6).command = 'sync';
% SeqControl(6).argument = 60e6; % Wait 60 seconds
%% 16. Event objects to acquire all acquisitions

% for iAcq=1:P.AverageNumber
%     Event(num_evt).info = 'Start Acquisition';
%     Event(num_evt).tx = 1;              % no TX structure.
%     Event(num_evt).rcv = iAcq;          % use 1st Rcv structure.
%     Event(num_evt).recon = 0;           % no reconstruction.
%     Event(num_evt).process = 0;         % no processing
%     Event(num_evt).seqControl = [2,3]; %  wait for trigger in to acquire + send trigger Out %+ transfer RF data to host
%     num_evt = num_evt + 1;
% end

% Event(num_evt).info = 'Custom Sw & Hw sync timeout.';   % sync timeout is .5s by default
% Event(num_evt).tx = 0;            % no TX structure.
% Event(num_evt).rcv = 0;           % no Rcv structure.
% Event(num_evt).recon = 0;         % no reconstruction.
% Event(num_evt).process = 0;       % no processing
% Event(num_evt).seqControl = 6;    % make the Sw wait at most 60s for the Hw to sync
% num_evt = num_evt + 1;

% Event(num_evt).info = 'Transfer RF Data to Host.';
% Event(num_evt).tx = 0;            % no TX structure.
% Event(num_evt).rcv = 0;           % no Rcv structure.
% Event(num_evt).recon = 0;         % no reconstruction.
% Event(num_evt).process = 0;       % no processing
% Event(num_evt).seqControl = 4;    % transfer to host
% num_evt = num_evt + 1;

% Perform reconstruction
for irecon = 1: P.AverageNumber 
Event(num_evt).info = 'Reconstruct';
Event(num_evt).tx = 0;            % no TX structure.
Event(num_evt).rcv = 0;           % no Rcv structure.
Event(num_evt).recon = irecon;         % reconstruction.
Event(num_evt).process = 1;       % reconstruction processing
Event(num_evt).seqControl = 2;    % trigger out
num_evt = num_evt + 1;
end

% % Jump back to Acquisition
% Event(num_evt).info = 'Jump back to acquisition.';
% Event(num_evt).tx = 0;            % no TX structure.
% Event(num_evt).rcv = 0;           % no Rcv structure.
% Event(num_evt).recon = 0;         % reconstruction.
% Event(num_evt).process = 0;       % no processing
% Event(num_evt).seqControl = 1;




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
VSXBF

%% 19. after quitting VSX, save the reconstructed ImgData 
if simulateMode == 2    
    disp ('Info:  Saving the Img Data and IQ Data buffers -- please wait!'), disp(' ')
    IQ=squeeze(IQData{1});
    Img = squeeze(ImgData{1});
    save ([FilePath,IQdataFilename],'Img', 'IQ', 'P')
end
FileName
delete(findobj('tag','UI'));
close all
% CmbnID
return
%% 20. **** Callback routines to be converted by text2cell function. ****
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

%% post-processing: remove jitter 
[nr, nc] = size(Img(:,:,1));
Imgk = zeros(nr+1,nc,P.AverageNumber);
Imgp = Img;
[~,MaxP] = find(Img(:,:,1)==max(max(Img(:,:,1)))); % column position of maximum value of the first image
for i = 1: P.AverageNumber
    indMax(i) = find(Img(:,MaxP,i)==max(Img(:,MaxP,i)));
    MaxMost = mode(indMax);
    if indMax(i)> MaxMost
%        k = indMax(i)-MaxMost;
%        Imgk(:,:,i) = padarray(Img(:,:,i), k, 0, 'post');
        Imgk(:,:,i) = vertcat(Img(:,:,i),Img(nr, :, 1));
        Imgp(:,:,i) = Imgk(2:end,:,i);
    elseif indMax(i) < MaxMost
        Imgk(:,:,i) = vertcat(Img(1, :, 1),Img(:,:,i));
        Imgp(:,:,i) = Imgk(1:end-1,:,i);
    end
end

for i = 1: P.AverageNumber
    figure(1);
    plot(Img(:,MaxP,i));
    hold on;
    figure(2);
    plot(Imgp(:,MaxP,i));
    hold on;
end

Imgbar = sum(Imgp, 3)/P.AverageNumber;
Imgvar = std(Imgp, 0, 3);
Imgsnr = 20*log10(Imgbar./Imgvar);
figure; imagesc(Imgbar); axis image; colorbar;% corrected averaged image
figure; imagesc(Img(:,:,1)); axis image; colorbar;% first image
figure; imagesc(Imgvar); axis image; colorbar;  % variance image
figure; imagesc(Imgsnr); axis image; colorbar;% snr image