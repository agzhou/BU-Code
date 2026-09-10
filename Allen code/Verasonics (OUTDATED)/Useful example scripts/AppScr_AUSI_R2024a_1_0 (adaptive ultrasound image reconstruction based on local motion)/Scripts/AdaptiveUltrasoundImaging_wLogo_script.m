%
% File name: AdaptiveUltrasoundImaging_wLogo_script.m - Example of reconstruction
%            using two different interleaved transmits. This script is
%            based on the AdaptiveUltrasoundImaging_script.m script.
%            This script is provided for users who wish to define a fixed, 
%            arbitrary region instead of using the region calculated by the
%            Doppler processing.
%
% Description:
%  This method uses pixel-oriented processing to adjust image reconstruction 
%  at individual spatial points. Two acquisition types were interleaved to 
%  concurrently acquire plane wave data for motion detection and high frame rate, 
%  and focused beam data for high image quality. Real-time spatial motion
%  sequence is the first part of this code. However, slow and fast
%  motion media regions are defined by the createVerasonicsV function, which
%  takes PData(1).Size(1) and PData(1).Size(2) as inputs and outputs a vector 
%  to replace the last PData(1).Region.PixelsLA variable. This parameter 
%  then defines the 'V' shape for the fast-motion region. The fast-motion
%  region (inside the 'V') uses plane wave transmits, resulting in lower
%  contrast and artifacts from side lobes, but the slow-motion region 
%  (outside the 'V') maintains higher image quality, contrast, and spatial 
%  resolution
%  The slow-motion region was reconstructed using accumulation of IQ data 
%  from 48 wide focused beams acquisitions, while fast motion region was 
%  reconstructed using flash plane wave acquisitions. 
%  The final image is a composite, with different areas reconstructed using 
%  varying amounts of data depending on the motion. 
%  ** Bmode and Motion detection sequence: Real-time spatial motion was
%            detected using speckle tracking. Based on a defined threshold, 
%            pixels were classified into "fast and slow-motion" regions.  
%     * Plane wave compounding for initial Bmode
%         + 5-angle transmit beams
%         + PRFtotal: 1.25kHz
%     * Fast media tracking
%         + 20 Tx/Rx per ensemble
%         + PRF: 3kHz
%     * Full 128-channel apodization for Tx/Rx
%  ** Adaptive ultrasound image (AUSI) sequence: Wide beam (WB) transmit scanning 
%            with spatial compounding. 
%     * 48 scan lines with 5 different steering directions (Based on 
%       SetUpL11_5gHWideBeamSC_opt.m).
%     * Frames are averaged using a running average that updates each frame. 
%         + Tx aperture size set to 32 elements. 
%         + Rx size set to 128 elements.
%     * An interleaved Plane Wave (PW) transmit and receive acquisition is 
%       performed too.
%         + All 128 Tx/Rx elements are active.
%         + Only use for region where motion was previously detected
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

% Parameters needed to determine kernel size for morphological operations
% to smooth the calculated slow-motion region
% radius1 = 3; dim1 = -radius1+1:radius1-1;
% SE1 = my_strel(radius1,dim1);

% P(1) is used for initial Bmode, P(2) is used for Doppler, local motion
P(1).startDepthMm = 1; % startDepth in mm
P(1).endDepthMm = 30;  % endDepth in mm
P(1).maxDepthMm = 70;  % maxDepth for RangeChange and RcvBuffer

% P(3) for High Image quality
P(3).numTx = 32;   % no. of elements in TX aperture.
P(3).numRays = 48; % no. of rays in frame
numAngles = 5;
P(3).dtheta = 12*(pi/180);  % angle delta for beams steered left or right
P(3).txFocusMm = 300;   % transmit focal pt in wavelengths


% Set 2D parameters for Media motion
na = 5;      % Set na = number of flash angles for 2D.
if (na > 1)
    dtheta2D = (30*pi/180)/(na-1);
    startAngle = -30*pi/180/2;
else
    dtheta2D = 0;
    startAngle=0;
end % set dtheta2D to range over +/- 15 degrees.

% Set Media motion parameters
ne = 16;     % Set ne = number of acquisitions in Media motion ensemble.
dopAngle = 12 * pi/180;
dopPRF = 3.0e+03; % Media motion PRF in Hz.
pwrThres = 0.27;  % Default threshold for Power Doppler

% Define Resources
AUSIFrames = 16; % number of frames in Receive buffer for Adaptive US Imaging (AUSI)

% Specify system parameters.
Resource.Parameters.numTransmit = 128;          % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;        % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.VDAS.watchdogTimeout = 120 * 1000;

% Specify Trans structure array.
Trans.name = 'L11-5gH';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % computeTrans is used for known transducers.
Trans.maxHighVoltage = 50;    % set a reasonable high voltage limit.

% Convert mm to wavelength
demodFreq = Trans.frequency; % demodulation frequency
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

% Calculate units in wavelengths
P(1).startDepth = P(1).startDepthMm*scaleToWvl;  % startDepth in wavelength
P(1).endDepth = P(1).endDepthMm*scaleToWvl;
P(1).maxDepth = P(1).maxDepthMm*scaleToWvl;
maxBufLength = ceil(sqrt(P(1).maxDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
maxBufSizePerAcq = 128*ceil(maxBufLength*8*(demodFreq/Trans.frequency)/128);

P(2).startDepth = P(1).startDepth;   % Acquisition depth in wavelengths
P(2).endDepth = P(1).endDepth ;   % This should preferrably be a multiple of 128 samples.

P(3).startDepth = P(1).startDepth;
P(3).endDepth = P(1).endDepth ;
P(3).txFocus = P(3).txFocusMm*scaleToWvl;

% Specify PData structure arrays.
% - PData(1) is defined for all BMode images including the initial Bmode
% generated by compounding five plane waved, 48 Wide Beams for high image
% quality and one Flash plane wave. Here we define the regions for each
% transmit
PData(1).PDelta = [0.7, 0, 0.5];
PData(1).Size(1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3)); % rows
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1)); % cols
PData(1).Size(3) = 1;             % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P(1).startDepth]; % x,y,z of uppr lft crnr.
PData(1).Region = repmat(struct('Shape',struct('Name','Rectangle',...
                                               'Position',[0,0,P(1).startDepth],...
                                               'width',5*Trans.spacing*Trans.numelements/(P(3).numRays),...
                                               'height',P(1).endDepth-P(1).startDepth)),1,1+numAngles*P(3).numRays+1); %Regions for for Flash, HIQ, AUSI
PData(1).Region(1).Shape.width = Trans.spacing*127;

%PData region for AUSI_WideBeams
% Compute the x coords of the TX beam centers
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P(3).numRays-1)):(63.5*Trans.spacing);
% Specify P.numRays rectangular regions centered on TX beam origins (use default angle of 0.0).
for n = 2:P(3).numRays+1, PData(1).Region(n).Shape.Position(1) = TxOrgX(n-1); end
m = P(3).numRays+1;
% Define numRays steered left parallelogram regions, centered on TX beam origins. Adjust the angle
%   so that the steering goes to zero over 8 beams at the left and right edge.
for n = 1:P(3).numRays
    if n<=8
        angle = -((n-1)/8)*P(3).dtheta;
    elseif n>(P(3).numRays-8)
        angle = -((P(3).numRays-n)/8)*P(3).dtheta;
    else
        angle = -P(3).dtheta;
    end
    PData(1).Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData(1).Region(n+m).Shape.height = (P(3).endDepth-P(3).startDepth)/cos(angle);
    PData(1).Region(n+m).Shape.angle = angle;
end
m = m + P(3).numRays;
% Define numRays steered right parallelogram regions, centered on TX beam origins. Adjust the angle
%   so that the steering goes to zero over 8 beams at the left and right edge.
for n = 1:P(3).numRays
    if n<=8
        angle = -((n-1)/8)*P(3).dtheta/2;
    elseif n>(P(3).numRays-8)
        angle = -((P(3).numRays-n)/8)*P(3).dtheta/2;
    else
        angle = -P(3).dtheta/2;
    end
    PData(1).Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData(1).Region(n+m).Shape.height = (P(3).endDepth-P(3).startDepth)/cos(angle);
    PData(1).Region(n+m).Shape.angle = angle;
end
m = m + P(3).numRays;
% Define numRays steered right parallelogram regions, centered on TX beam origins. Adjust the angle
%   so that the steering goes to zero over 8 beams at the left and right edge.
for n = 1:P(3).numRays
    if n<=8
        angle = ((n-1)/8)*P(3).dtheta/2;
    elseif n>(P(3).numRays-8)
        angle = ((P(3).numRays-n)/8)*P(3).dtheta/2;
    else
        angle = P(3).dtheta/2;
    end
    PData(1).Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData(1).Region(n+m).Shape.height = (P(3).endDepth-P(3).startDepth)/cos(angle);
    PData(1).Region(n+m).Shape.angle = angle;
end
m = m + P(3).numRays;
% Define numRays steered right parallelogram regions, centered on TX beam origins. Adjust the angle
%   so that the steering goes to zero over 8 beams at the left and right edge.
for n = 1:P(3).numRays
    if n<=8
        angle = ((n-1)/8)*P(3).dtheta;
    elseif n>(P(3).numRays-8)
        angle = ((P(3).numRays-n)/8)*P(3).dtheta;
    else
        angle = P(3).dtheta;
    end
    PData(1).Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData(1).Region(n+m).Shape.height = (P(3).endDepth-P(3).startDepth)/cos(angle);
    PData(1).Region(n+m).Shape.angle = angle;
end

% PData region for Flash plane wave transmits used for the Adaptive
% reconstruction
PData(1).Region(n+m+1).Shape.width = Trans.spacing*127;
PData(1).Region = computeRegions(PData(1));
PdataREF = PData; % saving the initial PData since it will change at each iteration
% Here user can use their own function instead of createVerasonicsV to
% create the slow and fast motion media regions
rows = PData(1).Size(1);
cols = PData(1).Size(2);
[~,pixels] = createVerasonicsV(rows, cols); 
PData(1).Region(1+numAngles*P(3).numRays+1).PixelsLA = int32(pixels);
PData(1).Region(1+numAngles*P(3).numRays+1).numPixels = length(PData(1).Region(1+numAngles*P(3).numRays+1).PixelsLA);
% Region1 = Rregion_total-Region2
R2PixelsLA = PData(1).Region(1+numAngles*P(3).numRays+1).PixelsLA;
for jj = 2:length(PData(1).Region)-1
    newR1PixelsLA = PData(1).Region(jj).PixelsLA;
    PData(1).Region(jj).PixelsLA  = setdiff(newR1PixelsLA,R2PixelsLA);
    PData(1).Region(jj).numPixels = length(PData(1).Region(jj).PixelsLA);
end

% PData for Media motion detection (Doppler)
PData(2).PDelta = [0.7, 0, 0.5];
PData(2).Size(1) = ceil((P(2).endDepth-P(2).startDepth)/PData(2).PDelta(3)); % rows
PData(2).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1)); % cols
PData(2).Size(3) = 1;             % single image page
PData(2).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P(2).startDepth]; % x,y,z of uppr lft crnr.

% Specify Media object for simulateMode
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
% - RcvBuffer(1) is for both intial BMode and media motion acquisitions.
Resource.RcvBuffer(1).rowsPerFrame = 2048*(na + ne);
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 4;           % 20 frames allocated for RF acqusitions.
% - RcvBuffer(2-6) is for AUSI acquisitions.
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = 2*P(3).numRays*maxBufSizePerAcq; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = AUSIFrames;
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = 2*P(3).numRays*maxBufSizePerAcq; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(3).numFrames = AUSIFrames;  % all RcvBuffers should have same no. of frames
Resource.RcvBuffer(4).datatype = 'int16';
Resource.RcvBuffer(4).rowsPerFrame = 2*P(3).numRays*maxBufSizePerAcq; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(4).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(4).numFrames = AUSIFrames;
Resource.RcvBuffer(5).datatype = 'int16';
Resource.RcvBuffer(5).rowsPerFrame = 2*P(3).numRays*maxBufSizePerAcq; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(5).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(5).numFrames = AUSIFrames;
Resource.RcvBuffer(6).datatype = 'int16';
Resource.RcvBuffer(6).rowsPerFrame = 2*P(3).numRays*maxBufSizePerAcq; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(6).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(6).numFrames = AUSIFrames;

% InterBuffer(1) is for Flash reconstructions.
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for 2D.
% InterBuffer(2) is for Media motion reconstructions.
Resource.InterBuffer(2).numFrames = 1;          % one intermediate frame needed for Doppler.
Resource.InterBuffer(2).pagesPerFrame = ne;     % ne pages per ensemble
% InterBuffer(3) is for AUSI_WB reconstructions.
Resource.InterBuffer(3).numFrames = 1;   % one intermediate buffer needed.
% InterBuffer(4) is for AUSI_Flash reconstructions.
Resource.InterBuffer(4).numFrames = 1;   % one intermediate buffer needed.

% ImageBuffer(1) is for 5-Plane Wave compounding image.
Resource.ImageBuffer(1).numFrames = 100;
% ImageBuffer(2) is for Media motion image.
Resource.ImageBuffer(2).numFrames = 30;
% ImageBuffer(3) is for AUSI_WB image.
Resource.ImageBuffer(3).numFrames = 300;
% ImageBuffer(4) for AUSI_Flash buffer
Resource.ImageBuffer(4).numFrames = 50*150;
% ImageBuffer(5) for final AUSI image
Resource.ImageBuffer(5).numFrames = 50*150;

% DisplayWindow 
Resource.DisplayWindow(1).Title = [Trans.name,'_Adaptive_US_image'];
Resource.DisplayWindow(1).pdelta = 0.3;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [100,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 50;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);
% Resource.DisplayWindow(1).splitPalette = 1;

% ------Specify structures used in Events------
% Specify Transmit waveform structures.
% - 2D Flash transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];
% - Media motion transmit waveform, transmit frequency should be equivalent to
% supported demodFrequency
TW(2).type = 'parametric';
TW(2).Parameters = [6.25,0.67,6,1];
% - AUSI wide beam transmit waveform
TW(3).type = 'parametric';
TW(3).Parameters = [Trans.frequency,0.67,1,1];
TW(3).equalize = 0;

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
end

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'TXPD', [], ...
                   'peakCutOff', 0.25, ...
                   'peakBLMax', 15.0, ...
                   'Delay', zeros(1,Trans.numelements)), 1, na+1+numAngles*P(3).numRays+1); % na TXs for Flash + 1 for Media Motion + WB + Flash
% - Set event specific TX attributes.
for n = 1:na   % na transmit events for Flash
    TX(n).Steer = [(startAngle+(n-1)*dtheta2D),0.0];
    TX(n).Delay = computeTXDelays(TX(n));  
end
% - only one TX struct needed for Media motion
TX(na+1).waveform = 2;
TX(na+1).Steer = [dopAngle,0.0];
TX(na+1).Delay = computeTXDelays(TX(na+1));
tx_pos = na+1;

% - transmits for AUSI Wide Beams
for n = 1:P(3).numRays  % specify P.numRays transmit events
    TX(n+tx_pos).Origin(1) = TxOrgX(n);
    TX(n+tx_pos).focus = P(3).txFocus;
    TX(n+tx_pos).waveform = 3;
    % Compute transmit aperture apodization
    TX(n+tx_pos).Apod = +(((scaleToWvl*Trans.ElementPos(:,1))>(TxOrgX(n)-Trans.spacing*P(3).numTx/2))& ...
                 ((scaleToWvl*Trans.ElementPos(:,1))<(TxOrgX(n)+Trans.spacing*P(3).numTx/2)))';
    [RIndices,CIndices,V] = find(TX(n+tx_pos).Apod);
    V = kaiser(size(V,2),2.5);
    % V = tukeywin(size(V,2),0.45);
    TX(n+tx_pos).Apod(CIndices) = V;
    % Compute transmit delays
    TX(n+tx_pos).Delay = computeTXDelays(TX(n+tx_pos));
end
m = P(3).numRays+tx_pos;
for n = 1:P(3).numRays
    TX(n+m).Origin(1) = TX(n+tx_pos).Origin(1);
    TX(n+m).focus = P(3).txFocus;
    TX(n+m).waveform = 3;
    TX(n+m).Apod = TX(n+tx_pos).Apod;
    if n<=8
        TX(n+m).Steer = [-((n-1)/8)*P(3).dtheta,0.0];
    elseif n>(P(3).numRays-8)
        TX(n+m).Steer = [-((P(3).numRays-n)/8)*P(3).dtheta,0.0];
    else
        TX(n+m).Steer = [-P(3).dtheta,0.0];
    end
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
end
m = m + P(3).numRays;
for n = 1:P(3).numRays
    TX(n+m).Origin(1) = TX(n+tx_pos).Origin(1);
    TX(n+m).focus = P(3).txFocus;
    TX(n+m).waveform = 3;
    TX(n+m).Apod = TX(n+tx_pos).Apod;
    if n<=8
        TX(n+m).Steer = [-((n-1)/8)*P(3).dtheta/2,0.0];
    elseif n>(P(3).numRays-8)
        TX(n+m).Steer = [-((P(3).numRays-n)/8)*P(3).dtheta/2,0.0];
    else
        TX(n+m).Steer = [-P(3).dtheta/2,0.0];
    end
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
end
m = m + P(3).numRays;
for n = 1:P(3).numRays
    TX(n+m).Origin(1) = TX(n+tx_pos).Origin(1);
    TX(n+m).focus = P(3).txFocus;
    TX(n+m).waveform = 3;
    TX(n+m).Apod = TX(n+tx_pos).Apod;
    if n<=8
        TX(n+m).Steer = [((n-1)/8)*P(3).dtheta/2,0.0];
    elseif n>(P(3).numRays-8)
        TX(n+m).Steer = [((P(3).numRays-n)/8)*P(3).dtheta/2,0.0];
    else
        TX(n+m).Steer = [P(3).dtheta/2,0.0];
    end
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
end
m = m + P(3).numRays;
for n = 1:P(3).numRays
    TX(n+m).Origin(1) = TX(n+tx_pos).Origin(1);
    TX(n+m).focus = P(3).txFocus;
    TX(n+m).waveform = 3;
    TX(n+m).Apod = TX(n+tx_pos).Apod;
    if n<=8
        TX(n+m).Steer = [((n-1)/8)*P(3).dtheta,0.0];
    elseif n>(P(3).numRays-8)
        TX(n+m).Steer = [((P(3).numRays-n)/8)*P(3).dtheta,0.0];
    else
        TX(n+m).Steer = [P(3).dtheta,0.0];
    end
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
end

% calculate TXPD
h = waitbar(0,'Program TX parameters, please wait!');
steps = numAngles*P(3).numRays;
for i = 1:steps
    TX(i+tx_pos).TXPD = computeTXPD(TX(i+tx_pos),PData(1));
    waitbar(i/steps)
end
close(h)

% Specify TX for AUSI Flash structure array.
n = numAngles*P(3).numRays + tx_pos + 1;
TX(n).waveform = 1;            % use 1st TW structure.
TX(n).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(n).focus = 0;
TX(n).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(n).Apod = ones(1,Trans.numelements);
TX(n).Delay = computeTXDelays(TX(n));
TX(n).TXPD = computeTXPD(TX(n),PData(1));

% Specify TPC structures.
TPC(1).name = 'Flash';
TPC(1).maxHighVoltage = 50;
TPC(2).name = 'MediaMotion';
TPC(2).maxHighVoltage = 35;

RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(2).LnaGain = 21; % Profile used for Media motion

% Specify TGC Waveform structures.
% - Flash TGC
TGC.CntrlPts = [200 750 800 1023 1023 1023 1023 1023];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - Media Motion TGC
TGC(2).CntrlPts = [10 950 950 1023 1023 1023 1023 1023];
TGC(2).rangeMax = P(2).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays.
%   We need to acquire all the 2D Flash and Media motion data within a single RcvBuffer frame.  This allows
%   the transfer-to-host DMA after each frame to transfer a large amount of data, improving throughput.
% - We need na Receives for a 2D frame and ne Receives for a Media motion frame.
maxAcqLngth2D = sqrt(P(1).endDepth^2 + (Trans.numelements*Trans.spacing)^2) - P(1).startDepth;
maxAcqLngthDop =  sqrt(P(2).endDepth^2 + (96*Trans.spacing)^2) - P(2).startDepth;
wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
wl2sPer128 = 128/(2*2);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
Receive1 = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P(1).startDepth, ...
                        'endDepth', maxAcqLngth2D, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ... % 200% Bandwidth for 2D
                        'mode', 0, ...
                        'demodFrequency',Trans.frequency,...
                        'callMediaFunc', 0), 1, (na+ne)*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (na + ne)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive1(k+1).callMediaFunc = 1;
    for j = 1:na  % acquisitions for 2D
        Receive1(j+k).framenum = i;
        Receive1(j+k).acqNum = j;
    end
    for j = (na+1):(na+ne)
        % Doppler acquisition
        Receive1(j+k).startDepth = P(2).startDepth;
        Receive1(j+k).endDepth = P(2).startDepth + wl2sPer128*ceil(maxAcqLngthDop/wl2sPer128);
        Receive1(j+k).sampleMode = 'BS100BW';
        Receive1(j+k).demodFrequency = TW(2).Parameters(1);
        Receive1(j+k).TGC = 2;
        Receive1(j+k).framenum = i;
        Receive1(j+k).acqNum = j;        % Media motion acqNums continue after 2D
    end
end

% Received buffer for AUSI, Wide Beams and Flash acquisitions
Receive2 = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P(3).startDepth, ...
                        'endDepth', maxAcqLngth2D, ...
                        'TGC', 1, ...
                        'bufnum', 2, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'demodFrequency',Trans.frequency,...
                        'callMediaFunc', 0), 1, (2*numAngles*P(3).numRays)*AUSIFrames);
% - Set event specific Receive attributes for each Bmode frame.
for ii = 1:AUSIFrames
%     Receive(2*P.numRays*(ii-1)+1).callMediaFunc = 1;
    for j = 1:2:2*P(3).numRays
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j).framenum = ii;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j).acqNum = j;

        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+2*P(3).numRays).bufnum = 3;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+2*P(3).numRays).framenum = ii;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+2*P(3).numRays).acqNum = j;

        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+4*P(3).numRays).bufnum = 4;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+4*P(3).numRays).framenum = ii;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+4*P(3).numRays).acqNum = j;

        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+6*P(3).numRays).bufnum = 5;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+6*P(3).numRays).framenum = ii;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+6*P(3).numRays).acqNum = j;

        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+8*P(3).numRays).bufnum = 6;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+8*P(3).numRays).framenum = ii;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+8*P(3).numRays).acqNum = j;
    end
    % - Set event specific Receive attributes for each AUSI frame.
    for j = 2:2:2*P(3).numRays
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j).framenum = ii;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j).acqNum = j;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j).callMediaFunc = 1;

        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+2*P(3).numRays).framenum = ii;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+2*P(3).numRays).acqNum = j;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+2*P(3).numRays).callMediaFunc = 1;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+2*P(3).numRays).bufnum = 3;

        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+4*P(3).numRays).framenum = ii;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+4*P(3).numRays).acqNum = j;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+4*P(3).numRays).callMediaFunc = 1;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+4*P(3).numRays).bufnum = 4;

        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+6*P(3).numRays).framenum = ii;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+6*P(3).numRays).acqNum = j;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+6*P(3).numRays).callMediaFunc = 1;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+6*P(3).numRays).bufnum = 5;

        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+8*P(3).numRays).framenum = ii;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+8*P(3).numRays).acqNum = j;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+8*P(3).numRays).callMediaFunc = 1;
        Receive2(numAngles*2*P(3).numRays*(ii-1)+j+8*P(3).numRays).bufnum = 6;
    end

end

Receive = [Receive1,Receive2];

% Specify Recon structure arrays.
% - We need two Recon structures, one for 2D, one for Media motion. These will be referenced in the same
%   event, so that they will use the same (most recent) acquisition frame.
Recon1 = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', zeros(1,1)), 1, 2);
% - Set Recon values for 2D frame.
Recon1(1).RINums(1,1:na) = (1:na);  % na ReconInfos needed for na angles
k = na + 1;
% - Set Recon values for Media motion ensemble.
Recon1(2).pdatanum = 2;
Recon1(2).IntBufDest = [2,1];
Recon1(2).ImgBufDest = [0,0];
Recon1(2).RINums(1,1:ne) = (k:(k+ne-1));   % ne ReconInfos needed for Media motion ensemble.
k = k+ne-1;

% Recon for AUSI Wide Beams
% - We need three Recon structures, one for each steering direction.
Recon2 = repmat(struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [3,1], ...
               'ImgBufDest', [3,-1], ...
               'RINums',1+k:2:2*P(3).numRays+k), 1, numAngles);
% - Set specific Recon attributes.
Recon2(2).RINums = (2*P(3).numRays+k+1):2:(2*P(3).numRays+2*P(3).numRays+k);
Recon2(3).RINums = (4*P(3).numRays+k+1):2:(2*P(3).numRays+4*P(3).numRays+k);
Recon2(4).RINums = (6*P(3).numRays+k+1):2:(2*P(3).numRays+6*P(3).numRays+k);
Recon2(5).RINums = (8*P(3).numRays+k+1):2:(2*P(3).numRays+8*P(3).numRays+k);

% Recon for AUSI Flash
Recon3 = repmat(struct('senscutoff', 0.8, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [4,1], ...
               'ImgBufDest', [4,-1], ...
               'RINums',2),1,numAngles*length(2:2:2*P(3).numRays));
for ii = 1:numAngles*length(2:2:2*P(3).numRays)
    Recon3(ii).RINums = ii*2+k;
end

Recon2 = [Recon2,Recon3];
Recon = [Recon1,Recon2];

% Define ReconInfo structures.
% - For 2D, we need na ReconInfo structures for na steering angles.
% - For Media motion, we need ne ReconInfo structures.
ReconInfo1 = repmat(struct('mode', 'accumIQ', ...    % accumulate IQ data.
                   'Pre',[],...
                   'Post',[],...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 0.2, ...
                   'regionnum', 1, ...
                   'pagenum',1, ...
                   'threadSync', 1), 1, na + ne);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
if na>1
    ReconInfo1(1).mode = 'replaceIQ'; % replace IQ data
    for j = 1:na
        ReconInfo1(j).txnum = j;
        ReconInfo1(j).rcvnum = j;
    end
    ReconInfo1(na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo1(1).mode = 'replaceIntensity';
end

%  - ReconInfos for Media motion ensemble.
k = na;
for j = 1:ne
    ReconInfo1(k+j).mode = 'replaceIQ';
    ReconInfo1(k+j).txnum = na + 1;
    ReconInfo1(k+j).rcvnum = na + j;
    ReconInfo1(k+j).pagenum = j;
end

% ReconInfo for AUSI Wide Beams
cont=1;
ReconInfo2 = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'Pre',[],...
                   'Post',[],...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 0.5, ...
                   'regionnum', 0, ...
                   'pagenum',1, ...
                   'threadSync', 1), 1, 2*P(3).numRays*numAngles);
% - Set specific ReconInfo attributes.
ReconInfo2(1).Pre = 'clearInterBuf';
for j = 1:2:2*P(3).numRays*numAngles
    ReconInfo2(j).txnum = 0.5*(j+1)+na+1;
    ReconInfo2(j).rcvnum = j+size(Receive1,2);
    ReconInfo2(j).regionnum = 0.5*(j+1)+1; %last +1 is for first region used on 2D initial image
    if sum(cont==(P(3).numRays:P(3).numRays:P(3).numRays*numAngles))==1
        ReconInfo2(j).Post = 'IQ2IntensityImageBuf';
        if (j+2)<2*P(3).numRays*numAngles
            ReconInfo2(j+2).Pre = 'clearInterBuf';
        end
    end
    cont = cont+1;
end

for j = 2:2:2*P(3).numRays*numAngles
    ReconInfo2(j).Pre = 'clearInterBuf';
    ReconInfo2(j).mode = 'replaceIntensity';
    ReconInfo2(j).Post = 'IQ2IntensityImageBuf';
    ReconInfo2(j).scaleFactor = 0.5;
    ReconInfo2(j).txnum = numAngles*P(3).numRays+1+na+1;
    ReconInfo2(j).rcvnum = j+size(Receive1,2);
    ReconInfo2(j).regionnum = numAngles*P(3).numRays+1+1;%last +1 is for first region used on 2D initial image
end

ReconInfo = [ReconInfo1,ReconInfo2];

% Specify Process structure arrays.
cpt = 150;  % define here so we can use in UIControl below
persf = 80;
persp = 90;
DopState = 'freq';

% Define all Process structures
ImgFlash = 1;
Process(ImgFlash).classname = 'Image';
Process(ImgFlash).method = 'imageDisplay';
Process(ImgFlash).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',2,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',20,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',60,...
                         'mappingMethod','full',...
                         'display',1,...      % don't display image after processing
                         'displayWindow',1};

MedMotion = ImgFlash+1;
Process(MedMotion).classname = 'Doppler';                   % process structure for 1st Doppler ensemble
Process(MedMotion).method = 'computeCFIPowerEst';
Process(MedMotion).Parameters = {'IntBufSrc',[2,1],...      % buffer and frame num of interbuffer
                         'SrcPages',[3,ne-2],...    % start and last pagenum
                         'ImgBufDest',[2,-1],...    % buffer and frame num of image desitation
                         'pdatanum',2,...           % number of PData structure
                         'prf',dopPRF,...           % Doppler PRF in Hz
                         'wallFilter','regression',...  % 1 -> quadratic regression
                         'pwrThreshold',pwrThres,...
                         'maxPower',50,...
                         'postFilter',1};

% Function GetMotion is commented since script uses another function to get
% the slow and fast media motion regions
% GetROI = MedMotion+1;
% Process(GetROI).classname = 'External';
% Process(GetROI).method = 'GetRegion';
% Process(GetROI).Parameters = {'srcbuffer','image',... % name of buffer to process.
%     'srcbufnum',2,...
%     'srcframenum',0};

Buf2Base = MedMotion+1;
Process(Buf2Base).classname = 'External';
Process(Buf2Base).method = 'ImgBuff1toBase';
Process(Buf2Base).Parameters = {'srcbuffer','image',... % name of buffer to process.
    'srcbufnum',3,...
    'srcframenum',-1,...
    'dstbuffer','none'};

MergBuff = Buf2Base+1;
Process(MergBuff).classname = 'External';
Process(MergBuff).method = 'mergeImgBuffers';
Process(MergBuff).Parameters = {'srcbuffer','image',... % name of buffer to process.
    'srcbufnum',4,...
    'srcframenum',-1,...
    'dstbuffer','image',...
    'dstbufnum',5,...
    'dstframenum',-1 };

ImgAUSI1 = MergBuff+1;
Process(ImgAUSI1).classname = 'Image';
Process(ImgAUSI1).method = 'imageDisplay';
Process(ImgAUSI1).Parameters = {'imgbufnum',3,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',45,...            % pgain is image processing gain
                         'reject',0,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',20,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','runAverage5',...
                         'compressMethod','log',...
                         'compressFactor',30,...
                         'mappingMethod','lowerHalf',...
                         'display',0,...      % don't display image after processing
                         'displayWindow',1};

ImgAUSI2 = ImgAUSI1+1;
Process(ImgAUSI2).classname = 'Image';
Process(ImgAUSI2).method = 'imageDisplay';
Process(ImgAUSI2).Parameters = {'imgbufnum',4,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',0.15,...            % pgain is image processing gain
                         'reject',0,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',20,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','log',...
                         'compressFactor',30,...
                         'mappingMethod','upperHalf',...
                         'display',1,...      % don't display image after processing
                         'displayWindow',1};

% external function for UI control
ExtUI = ImgAUSI2+1;
Process(ExtUI).classname = 'External';
Process(ExtUI).method = 'UIControl';
Process(ExtUI).Parameters = {'srcbuffer','none'}; 

% external function move to AUSI
Move2ausi = ExtUI+1;
Process(Move2ausi).classname = 'External';
Process(Move2ausi).method = 'Move2AUSI';
Process(Move2ausi).Parameters = {'srcbuffer','none'};

% Specify SeqControl structure arrays.
% -- Time between 2D flash angle acquisitions
Seq_TTNA_2D = 1;
SeqControl(Seq_TTNA_2D).command = 'timeToNextAcq';
SeqControl(Seq_TTNA_2D).argument = 160; % time in usec
% -- Change to Profile 2 (Media motion)
Seq_TPC2 = Seq_TTNA_2D+1;
SeqControl(Seq_TPC2).command = 'setTPCProfile';
SeqControl(Seq_TPC2).condition = 'next';
SeqControl(Seq_TPC2).argument = 2;
% -- Time between 2D acquisition and Media motion ensemble. Set to allow time for profile change.
Seq_TTNA2MM = Seq_TPC2+1;
SeqControl(Seq_TTNA2MM).command = 'timeToNextAcq';
SeqControl(Seq_TTNA2MM).argument = 8000; % time in usec
% -- PRF for Media motion ensemble
Seq_TTNA_MM = Seq_TTNA2MM+1;
SeqControl(Seq_TTNA_MM).command = 'timeToNextAcq';
SeqControl(Seq_TTNA_MM).argument = round(1/(dopPRF*1e-06)); % (for 3KHz dopPRF & if 14 ensemble = 4.7 msecs)
% -- Change to Profile 1 (2D)
Seq_TPC1 = Seq_TTNA_MM+1;
SeqControl(Seq_TPC1).command = 'setTPCProfile';
SeqControl(Seq_TPC1).condition = 'next';
SeqControl(Seq_TPC1).argument = 1;
% -- Time between Media motion and next 2D acquisition. Set to allow time for profile change.
Seq_TTNA2N2D = Seq_TPC1+1;
SeqControl(Seq_TTNA2N2D).command = 'timeToNextAcq';
SeqControl(Seq_TTNA2N2D).argument = 8000; % time in usec
% -- Jump back to start.
Seq_jump1 = Seq_TTNA2N2D+1;
SeqControl(Seq_jump1).command = 'jump';
SeqControl(Seq_jump1).argument = 1;
% set receive profile
Seq_RCVP1 = Seq_jump1+1;
SeqControl(Seq_RCVP1).command = 'setRcvProfile';
SeqControl(Seq_RCVP1).argument = 1;
Seq_RCVP2 = Seq_RCVP1+1;
SeqControl(Seq_RCVP2).command = 'setRcvProfile';
SeqControl(Seq_RCVP2).argument = 2;
% -- Time between AUSI acquisitions (Wide Beam and Flash)
Seq_TTNA_AUSI = Seq_RCVP2+1;
SeqControl(Seq_TTNA_AUSI).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(Seq_TTNA_AUSI).argument = 160;  % in usec
% -- Time bewteen AUSI frames
Seq_TTNA_AUSIf = Seq_TTNA_AUSI+1;
SeqControl(Seq_TTNA_AUSIf).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(Seq_TTNA_AUSIf).argument = 20000 - (2*P(3).numRays-1)*SeqControl(Seq_TTNA_AUSI).argument;  % usec time between frames
% -- Return to Matlab
Seq_RTMat = Seq_TTNA_AUSIf+1;
SeqControl(Seq_RTMat).command = 'returnToMatlab';
% - The remainder of the SeqControl structures are defined dynamically in the sequence events.
%   The variable nsc keeps track of the next SeqControl structure index.
nsc = length(SeqControl)+1; % Count of SeqControl objects.

% Specify Event structure arrays.
frameRateFactor = 2;
n = 1;
Event(n).info = 'ext func for UI control';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = ExtUI;
Event(n).seqControl = 0;
n = n+1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire 2D frame
    for j = 1:na
        Event(n).info = 'Acquire 2D flash angle';
        Event(n).tx = j;
        Event(n).rcv = (na+ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = Seq_TTNA_2D;
        if j == 1
            Event(n).seqControl = [Seq_TTNA_2D,Seq_RCVP1];
        end
        n = n+1;
    end
    Event(n-1).seqControl = [Seq_TPC2,Seq_TTNA2MM];   % replace last 2D acquisition Event's seqControl
    % Acquire Doppler ensemble.
    for j = (na+1):(na+ne)
        Event(n).info = 'Acquire Media motion ensemble';
        Event(n).tx = na+1;
        Event(n).rcv = (na+ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = Seq_TTNA_MM;
        if j == na+1
            Event(n).seqControl = [Seq_TTNA_MM,Seq_RCVP2];
        end
        n = n+1;
    end
    Event(n-1).seqControl = [Seq_TPC1,Seq_TTNA2N2D,nsc]; % replace last Media motion acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recons and 2D process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = [1,2];
    Event(n).process = ImgFlash;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Media motion processing';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = MedMotion;
    Event(n).seqControl = Seq_RTMat;
%     if floor(i/frameRateFactor) == i/frameRateFactor && i ~= Resource.RcvBuffer(1).numFrames
%         Event(n).seqControl = Seq_RTMat;
%     end
    n = n+1; 
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = Seq_jump1;
n = n+1;

% Once AUSI button is on, sequence start point moves here, now we are in
% the Adaptive US Image loop sequence. It starts calling the external
% function that calculates the two regions: slow and fast motion regions
% using data collected in the intial loop above.
AUSI_seq = n;
% This Event is commented since we are not using the GetRegion function
% Event(n).info = 'Get region';
% Event(n).tx = 0;
% Event(n).rcv = 0;
% Event(n).recon = 0;
% Event(n).process = GetROI;
% Event(n).seqControl = 0;
% n = n+1;

Event(n).info = 'Move to AUSI loop';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = Move2ausi;
Event(n).seqControl = Seq_RTMat;
n = n+1;

AUSI_loop = n;
for ii = 1:AUSIFrames
     for hh = 1:numAngles
        ww = 2*P(3).numRays*(hh-1);
        for j = 1:2:2*P(3).numRays       
            Event(n).info = 'AUSI_WideBeam';
            Event(n).tx = 0.5*(j+1)+(na+1)+ww*0.5;
            Event(n).rcv = (numAngles*2*P(3).numRays)*(ii-1)+j+size(Receive1,2)+ww;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = Seq_TTNA_AUSI;
            if j == 1
                Event(n).seqControl = [Seq_TTNA_AUSI,Seq_RCVP1,Seq_TPC1];
            end

            Event(n+1).info = 'AUSI_Flash';
            Event(n+1).tx = numAngles*P(3).numRays+1+na+1;
            Event(n+1).rcv = (numAngles*2*P(3).numRays)*(ii-1)+j+1+size(Receive1,2)+ww;
            Event(n+1).recon = 0;
            Event(n+1).process = 0;
            Event(n+1).seqControl = Seq_TTNA_AUSI;

            n = n+2;
        end
        Event(n-1).seqControl = [Seq_TTNA_AUSIf,nsc]; % modify last acquisition Event's seqControl
        SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
        nsc = nsc+1;

        n = length(Event)+1;
        Event(n).info = 'recon for slow FR';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 2+hh;
        Event(n).process = 0;
        Event(n).seqControl = 0;
        n = n+1;

        Event(n).info = 'Image Buff1 to Base';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = Buf2Base;
        Event(n).seqControl = 0;
        n = n+1;

        for mm = 1:length(2:2:2*P(3).numRays)
            Event(n).info = 'recon high FR';
            Event(n).tx = 0;
            Event(n).rcv = 0;
            Event(n).recon = mm+3+numAngles+ww*0.5-1;
            Event(n).process = 0;
            Event(n).seqControl = 0;
            n = n+1;

            Event(n).info = 'merge ImgBuffers';
            Event(n).tx = 0;
            Event(n).rcv = 0;
            Event(n).recon = 0;
            Event(n).process = MergBuff;
            Event(n).seqControl = 0;
            n = n+1;

            Event(n).info = 'Process AUSI image from Buff3';
            Event(n).tx = 0;
            Event(n).rcv = 0;
            Event(n).recon = 0;
            Event(n).process = ImgAUSI1;
            Event(n).seqControl = 0;
            n = n+1;
         
            Event(n).info = 'Process AUSI image from Buff4';
            Event(n).tx = 0;
            Event(n).rcv = 0;
            Event(n).recon = 0;
            Event(n).process = ImgAUSI2;
            Event(n).seqControl = 0;
            n = n+1;
        end

        n = n-1;
        if floor(ii/2) == ii/2     % Exit to Matlab every 3rd frame reconstructed
            Event(n).seqControl = Seq_RTMat;
        end
        n = n+1;
    end

end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc;
  SeqControl(nsc).command = 'jump';
  SeqControl(nsc).argument = AUSI_loop;
  nsc = nsc+1;

% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Doppler Power Threshold Slider
UI(2).Control = {'UserB3','Style','VsSlider','Label','DopPwrThres','SliderMinMaxVal',[0.0,1.0,pwrThres],...
                 'SliderStep',[0.02,0.1],'ValueFormat','%3.2f'};
UI(2).Callback = text2cell('%-UI#3Callback');

% - Peak CutOff
UI(3).Control = {'UserB2','Style','VsSlider','Label','Peak Cutoff',...
                  'SliderMinMaxVal',[0,20.0,TX(1).peakCutOff],...
                  'SliderStep',[0.005,0.020],'ValueFormat','%1.3f'};
UI(3).Callback = text2cell('%PeakCutOffCallback');

% - Controls to move between Media Motion detection and AUSI method
UI(4).Control = vsv.seq.uicontrol.VsButtonControl('LocationCode','UserB4','Label','Get_media' );
UI(5).Control = vsv.seq.uicontrol.VsButtonControl('LocationCode','UserB4','Label','AUSI');

% Define external functions
EF(1).Function = vsv.seq.function.ExFunctionDef('UIControl', @UIControl);
EF(2).Function = vsv.seq.function.ExFunctionDef('Move2AUSI', @Move2AUSI);
EF(3).Function = vsv.seq.function.ExFunctionDef('ImgBuff1toBase', @ImgBuff1toBase);
EF(4).Function = vsv.seq.function.ExFunctionDef('mergeImgBuffers', @mergeImgBuffers);
% EF(5).Function = vsv.seq.function.ExFunctionDef('GetRegion', @GetRegion);

% Save all the structures to a .mat file.
clear Receive1 Receive2 Recon1 Recon2 Recon3 ReconInfo1 ReconInfo2
save('MatFiles/L11-5gH_AUSI_sequence_logo');
return

% **** Callback routines to be encoded by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
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
return
%SensCutoffCallback

%-UI#3Callback - Doppler Power change
Process = evalin('base','Process');
for k = 1:2:length(Process(2).Parameters)
    if strcmp(Process(2).Parameters{k},'pwrThreshold'), Process(2).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Doppler threshold.
Control = evalin('base','Control');
if isempty(Control(1).Command), n=1; else, n=length(Control)+1; end
Control(n).Command = 'set&Run';
Control(n).Parameters = {'Process',2,'pwrThreshold',UIValue};
assignin('base','Control', Control);
%-UI#3Callback

%PeakCutOffCallback
TX = evalin('base', 'TX');
for ii=1:size(TX,2)
    TX(ii).peakCutOff = UIValue;
end
assignin('base','TX',TX);
% Set Control.Command to set TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
%PeakCutOffCallback

function UIControl(varargin)
    UI = evalin('base','UI');
    set(UI(4).handle,'Callback',@Get_media);
    set(UI(5).handle,'Callback',@AUSI);
end

function AUSI(varargin)
    % make related UI controls visible
    UI = evalin('base','UI');
    for i = 4:5
        set(UI(i).handle,'Visible','on','Interruptible','off');
    end
    set(UI(5).handle,'Visible','off');
    set(UI(2).handle.List{1},'Enable', 'off'); % disable playing with Doppler threshold again
    set(UI(2).handle.List{2},'Enable', 'off');
    set(UI(2).handle.List{3},'Enable', 'off');
    
    % change the start event
    nStart = evalin('base','AUSI_seq');
    Control = evalin('base','Control');
    if isempty(Control(1).Command), n=1; else n =length(Control)+1; end
    Control(n).Command = 'set&Run';
    Control(n).Parameters = {'Parameters',1,'startEvent',nStart};
    evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
    assignin('base','Control',Control);
end

function Get_media(varargin)
    % make related UI controls unvisible
    UI = evalin('base','UI');
    for i = 4:5
        set(UI(i).handle,'Visible','off');
    end
    set(UI(5).handle,'Visible','on');
    set(UI(2).handle.List{1},'Enable', 'on'); % enable playing with Doppler threshold again
    set(UI(2).handle.List{2},'Enable', 'on');
    set(UI(2).handle.List{3},'Enable', 'on');
    
    % change the start event
    nStart = 1;
    
    Control = evalin('base','Control');
    if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
    Control(n).Command = 'set&Run';
    Control(n).Parameters = {'Parameters',1,'startEvent',nStart};
    evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
    assignin('base','Control',Control);
end

function Move2AUSI(varargin)
    % change the start event
    nStart = evalin('base','AUSI_loop');
    
    Control = evalin('base','Control');
    if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
    Control(n).Command = 'set&Run';
    Control(n).Parameters = {'Parameters',1,'startEvent',nStart};
    evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
    assignin('base','Control',Control);
end

function ImgBuff1toBase(imageBuff1)
  assignin('base','ImgBkg',imageBuff1);
end

function ImgFinal = mergeImgBuffers(imageBuff2)
   ImgBkg = evalin('base','ImgBkg');
   ImgFinal = ImgBkg + imageBuff2;
end
