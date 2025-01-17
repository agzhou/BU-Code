clear all

thisFileName = 'RunSetUpL12_3VPAT_Acquire_LASER_MASTER_test';
cd C:\Users\BOAS-US\Documents\Vantage-3.4.0-1711281030

P.startDepth = 1;   % Acquisition depth in wavelengths
P.endDepth = 300;   % This should preferrably be a multiple of 128 samples.

% Define system parameters.
filename = mfilename; % used to launch VSX automatically
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.connector = 2;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

trigger = true;
transmit = false;

% Resource.Parameters.verbose = 2;
% Specify media points
Media.MP(1,:) = [0,0,100,1]; % [x, y, z, reflectivity]
Media.MP(2,:) = [20,0,80,.5]; % [x, y, z, reflectivity]
Media.MP(3,:) = [-30,0,150,.2]; % [x, y, z, reflectivity]
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Media.function = 'movePoints';

%%
% Specify Trans structure array.
Trans.name = 'L12-3v';
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 50; % mfr data sheet lists 30 Volt limit
P.LensDelay=Trans.lensCorrection;
P.frequency=Trans.frequency;
P.pitch=Trans.spacingMm;

% Specify PData structure array.
PData(1).PDelta = [0.4, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Resource.Parameters.numTransmit*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Resource.Parameters.numTransmit-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

%%
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Resource.Parameters.numRcvChannels-1)*Trans.spacing)^2));
actAcqLength =(maxAcqLength-P.startDepth/2);P.nSperWave=4;
P.actZsamples=ceil(actAcqLength * 2 * P.nSperWave/128)*128; % P.nSperWave: sampling points per wavelength, rounded up to the next 128 VSX sample boundary


% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.actZsamples;%6144; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 1;       % 100 frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer neededm
Resource.ImageBuffer(1).numFrames = 1;

% Display
Resource.DisplayWindow(1).Title = 'L22-14v PAT Acquisition, 4X sampling at 62.5 MHz';
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

%%
% Specify Transmit waveform structure.  
TW.type = 'parametric';
TW.Parameters = [18, 0.67, 2, 1];   % 18 MHz center frequency, 67% pulsewidth two cycle burst

% Specify TX structure array.
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



% Specify TGC Waveform structure.
TGC(1).CntrlPts = [500,590,650,710,770,830,890,950];
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays. 
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
     
     
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Resource.Parameters.numTransmit-1)*Trans.spacing)^2));
Receive = struct(   'Apod', ones(1,Resource.Parameters.numTransmit), ...
                    'aperture',33, ...
                    'startDepth', P.startDepth, ...
                    'endDepth', maxAcqLength,...
                    'TGC', 1, ...
                    'bufnum', 1, ...
                    'framenum', 1, ...
                    'acqNum', 1, ...
                    'sampleMode', 'NS200BW', ...
...%                     'InputFilter', BPF1, ... 
                    'mode', 0, ...
                    'callMediaFunc', 1);

%%
% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', 1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1 ,...
               'newFrameTimeout', -1);

% Define ReconInfo structures.
ReconInfo = struct('mode', 'replaceIntensity', ...  
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1);

Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level 
                         'interpMethod','4pt',...  
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

%%
% Specify sequence events.

num_evt = 1;

Event(num_evt).info = 'Acquire RF Data.';
Event(num_evt).tx = 1;            % use 1st TX structure.
Event(num_evt).rcv = 1;           % use 1st Rcv structure.
Event(num_evt).recon = 0;         % no reconstruction.
Event(num_evt).process = 0;       % no processing
Event(num_evt).seqControl = [2,6,4]; %  wait for trigger in to acquire + transfer RF data to host


acq_evt = num_evt;
num_evt = num_evt + 1;

Event(num_evt).info = 'Custom Sw & Hw sync timeout.';   % sync timeout is .5s by default
Event(num_evt).tx = 0;            % use 1st TX structure.
Event(num_evt).rcv = 0;           % use 1st Rcv structure.
Event(num_evt).recon = 0;         % no reconstruction.
Event(num_evt).process = 0;       % no processing
Event(num_evt).seqControl = 5;    % make the Sw wait at most 60s for the Hw to sync

num_evt = num_evt + 1;

% Event(num_evt).info = 'Transfer RF Data to Host.';
% Event(num_evt).tx = 0;            % use 1st TX structure.
% Event(num_evt).rcv = 0;           % use 1st Rcv structure.
% Event(num_evt).recon = 0;         % no reconstruction.
% Event(num_evt).process = 0;       % no processing
% Event(num_evt).seqControl = 4;    % transfer data to host
% 
% num_evt = num_evt + 1;
% 
% Perform reconstruction
Event(num_evt).info = 'Reconstruct RF Data.';
Event(num_evt).tx = 0;            % no TX structure.
Event(num_evt).rcv = 0;           % no Rcv structure.
Event(num_evt).recon = 1;         % reconstruction.
Event(num_evt).process = 1;       % processing
Event(num_evt).seqControl = 3;

num_evt = num_evt + 1;

% Jump back to Acquisition
Event(num_evt).info = 'Jump to acquisition.';
Event(num_evt).tx = 0;            % no TX structure.
Event(num_evt).rcv = 0;           % no Rcv structure.
Event(num_evt).recon = 0;         % reconstruction.
Event(num_evt).process = 0;       % no processing
Event(num_evt).seqControl = 1;

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = acq_evt;
if trigger
    SeqControl(2).command = 'triggerIn';  % time between frames
    SeqControl(2).condition = 'Trigger_2_Rising';  % Wait for trigger 1 rising edge
%     SeqControl(2).argument = 0;
else
    SeqControl(2).command = 'timeToNextAcq';  % time between frames
    SeqControl(2).argument = 50000;  % 50 msec
end
SeqControl(3).command = 'returnToMatlab';

SeqControl(4).command = 'transferToHost';

SeqControl(5).command = 'sync';
SeqControl(5).argument = 60e6; % Wait 1 minute % microsecond

SeqControl(6).command = 'triggerOut';
SeqControl(6).argument = 0;

% Save all the structures to a .mat file.

save(['MatFiles/',thisFileName]); 
disp([ mfilename ': NOTE -- Running VSX automatically!']), disp(' ')
VSX
%%
%     disp ('Info:  Saving the RF Data buffer -- please wait!'), disp(' ')
%     if exist([savepath,RFdataFilename,'.mat'])
%         NewPlane=1;
%         while (exist([savepath,RFdataFilename,'.mat'])==2)
%             NewPlane=NewPlane+1;
%             RFdataFilename = ['PATRF-',num2str(P.CCangle),'-',num2str(P.numAngles),'-',num2str(P.CCFR),'-',num2str(P.numCCframes),'-',num2str(P.numSupFrames),'-',Filename,'-',num2str(NewPlane)];
%         end
%         prompt={'File exist! New Plane Sequence: '};
%         name='File info';
%         defaultvalue={num2str(NewPlane)};
%         inputValue=inputdlg(prompt,name, 1, defaultvalue);
%         P.Plane=str2num(inputValue{1});
%         RFdataFilename = ['RF-',num2str(P.CCangle),'-',num2str(P.numAngles),'-',num2str(P.CCFR),'-',num2str(P.numCCframes),'-',num2str(P.numSupFrames),'-',Filename,'-',num2str(P.Plane)];
%     tic
%     RFRAW=RcvData{1};
%     savefast ([savepath,RFdataFilename], 'RFRAW','P')
%     toc
%     disp ('RF DATA SAVED!')
% end
% return
