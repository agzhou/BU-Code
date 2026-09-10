% Take note of spatial and temporal resolution. Note that the time lagged
% stuff needs more subframes so the frame rate this displays is not
% necessarily accurate

% for straight up and down, t_tx_to_rcv = 2*h/c = 2*endDepth_wl*(c/f_c)/c =
% 2*endDepth/(c^2 * f_c)

% for endDepth = 100 and f_c = 13.6 MHz, t = 6.200827438413382e-12 s

%% Specify system parameters
clear

% cd 'C:\Users\BOAS-US\Desktop\Vantage-4.9.5-2409181500'
cd 'G:\My Drive\Verasonics files\Vantage-4.9.2-2308102000'
activate

% savepath = 'G:\Allen\Data\10-03-2024 phantom\L22-14v\';

movePointsOrNot = 0; % Media movePoints on or off
runVSX = 1;
simOrNot = 1;

fps_target = 1;
supFrameBurstRate = 1;

startDepthMM = 0; % start depth in mm
endDepthMM = 10;

initialVoltage = 20;


numSupFrames = 1;
numSubFrames = 1;
na = 5; % # of acquisitions per frame (# angles)
maxAngle = 5; % degrees
angleRange = [-maxAngle, maxAngle].*pi/180; % Angle range in radians

% savepath = strcat("G:\Allen\Data\10-04-2024 phantom\L22-14v\", "verex2 TGC ", num2str(na), " angles ", num2str(numFrames), " frames ", num2str(endDepthMM), " endDepthMM ", num2str(maxAngle), " maxAngle test\");
% savepath = 'D:\Allen\Simulation Results\L22-14v\datatest\'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% savepath = "G:\Allen\Data\12-16-2024 Phantom\L22-14v\run 3 42 angles -10 to 10 deg\";
% savepath = char(savepath);
% mkdir(savepath)

if na >= 2 % Need at least 2 acquisitions to use multiple angles
    angles = linspace(angleRange(1), angleRange(2), na);
else
    angles = 0;
end

% savepath = strcat("G:\Allen\Data\01-09-2025 phantom anechoic\L22-14v\run 9 ", num2str(na), " angles -", num2str(maxAngle), " to ", num2str(maxAngle), " deg\");
% savepath = char(savepath);
savepath = 'D:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\PSF sim\';
mkdir(savepath)

numChannels = 128;
% Resource is a structure, define system parameters
Resource.Parameters.numTransmit = numChannels; % number of transmit channels
Resource.Parameters.numRcvChannels = numChannels; % number of receive channels
% Resource.Parameters.connector = 1; % transducer connector to use since the current plate for the 256 bit system is split into two 128 bit connectors. 1 is left and 2 is right
Resource.Parameters.speedOfSound = 1540; % speed of sound in m/s, the 1540 is for average human tissue

%% Define Transducer

Trans.name = 'L22-14v';
Trans.frequency = 13.6; % Not needed if using the default center frequency
Trans.units = 'wavelengths'; % or mm
% Trans.units = 'mm';

Trans = computeTrans(Trans); % Generate required attributes for the probe into the Trans structure; e.g., the transducer element positions
% Trans.maxHighVoltage = ; % set maximum high voltage that is allowed to the transducer

% Bae and Song 2018 test for uniform angle condition (runs out of memory)
L = Trans.spacingMm*Trans.numelements/1e3; % m
wl = Resource.Parameters.speedOfSound / Trans.frequency / 1e6; % m
% xf = 0;
% % dt_max = wl / (L/2 + xd * tand(maxAngle) - xf);
% dt_max = wl / (L/2 - xf);
% angles = [-maxAngle:dt_max:maxAngle].*pi/180;
% numAngles = length(angles);
% na = length(angles);
% 
% GL_appear_interval_mm = wl./(maxAngle*pi/180*2/na)*1e3
% NN = 6;
% N1 = NN;
% aR1 = [-maxAngle, maxAngle].*pi/180; % angle Range 1
% aR2 = [-maxAngle, maxAngle].*pi/180.*2./3;
% N2 = NN;
% da_1 = (aR1(end)-aR1(1))./N1;
% da_2 = (aR2(end)-aR2(1))./N2;
% a1 = aR1(1):da_1:aR1(end);
% a2 = aR2(1):da_2:aR2(end);
% angles = [a1 a2]
% na = length(angles);
% numAngles = na;

%% Simulation things - Media structure

Resource.Parameters.simulateMode = simOrNot; % run script in simulate mode. Set to 0 if not

xl = 40;
zl = 50;
zlmult = 1.5;
% Set up Media model for the simulation, which generates the scattering points with 3D location
% and reflectivity. For 1D transducer arrays, they are aligned on the
% x-axis with the center at x = 0, and scan depth is in z.
% Media.MP(1,:) = [-10, -10, zl, 1.0]; % [x, y, z, reflectivity]. x, y, z are defined as # of wavelengths.
% Media.MP(2,:) = [10, -10, zl, 1.0];
% Media.MP(3,:) = [10, 10, zl, 1.0];
% Media.MP(4,:) = [-10, 10, zl, 1.0];

Media.MP(1,:) = [0, 0, zl, 1.0];
% Media.MP(2,:) = [xl, 0, zl, 1.0];
% Media.MP(3,:) = [-xl, 0, zl, 1.0];
% Media.MP(4,:) = [0, 0, zl*zlmult, 1.0];
% Media.MP(5,:) = [xl, 0, zl*zlmult, 1.0];
% Media.MP(6,:) = [-xl, 0, zl*zlmult, 1.0];
% Media.MP(7,:) = [0, 0, zl*2, 1.0];
% Media.MP(8,:) = [xl, 0, zl*2, 1.0];
% Media.MP(9,:) = [-xl, 0, zl*2, 1.0];

% Media.MP(1,:) = [-xl, 0, zl, 1.0];

% 
% testmp = 4/1e3;
% Media.MP(4,:) = [testmp/wl, 0, zl*2, 1.0];
% Media.MP(5,:) = [-testmp/wl, 0, zl*2, 1.0];

% Media.attenuation = 0;
Media.attenuation = -0.7; %  media attenuation in dB/cm/MHz
% Media.FlowObj = struct();

% randomPts;

% Media.function = 'movePointsX'; % move points in x after each frame
Media.function = 'movePointsZ';


startDepth = ceil(startDepthMM/1e3/wl); % in wavelengths
endDepth = ceil(endDepthMM/1e3/wl); % in wavelengths

numElements = Trans.numelements; % the structure gives # row elements + # column elements

%% Transmission Waveform (TW)
tw.A = Trans.frequency; % frequency of transmission pulse, sets half cycle period of the waveform...
tw.B = 0.67; % amount of time (0.1 - 1.0) that the transmission drivers are active in the half cycle period. Controsl how much power is delivered.
             % Apparently using B = 0.67 approximates a sine wave.
tw.C = 2; % number of half cycles in the transmission waveform. 2 half cycles = 1 full cycle burst
tw.D = 1; % initial polarity of the first half cycle (1 = +, 0 = -)
TW(1).type = 'parametric';
TW(1).Parameters = [tw.A, tw.B, tw.C, tw.D];

% Note: can modify B for transmission apodization......... See tutorial, "The modification of the B parameters for Transmit apodization can be
% performed automatically by VSX by setting the weighting values in the TX.Apod attribute of a TX structure..."

% In the case where one really would like independent parametric waveforms on individual
% channels, the TW.Parameters array can be expanded to a 128 (or 256) row, two
% dimensional array, where each row specifies the waveform for the transmitter of the same
% number as the row index.
% The Vantage system supports other methods for defining transmit waveforms, including the
% TW.types of ‘envelope’, ‘pulseCode’, and ‘function’. These additional methods are
% described in the Sequence Programming Manual.

TPC.hv = initialVoltage;
%% Transmit action - TX structure

% Need a TX structure for each unique transmit action in the imaging
% sequence

% na*2 transmissions of a plane wave in pairs, one by all row elements and then one by all
% column elements
TX = repmat(struct('waveform', 1, ...
                   'focus', 0, ... % plane wave
                   'Steer', [0.0, 0.0], ... % theta, alpha (beam angle projected in xz from +z axis, beam angle wrt xz)
                   'Apod', ones(1, Trans.numelements)), 1, na);
for n = 1:na
    TX(n).Steer = [angles(n), 0];
    TX(n).Delay = computeTXDelays(TX(n));
end

%% Define Time Gain Control waveform (TGC)
% Accounts for decrease in amplitude of echoes for longer distance traveled

% TGC curve definition
% TGC(1).CntrlPts = [500,590,650,710,770,830,890,950]; % 0 to 1023, minimum to maximum gain
                                                     % Values represent the
                                                     % gain at increasing
                                                     % depth in the
                                                     % acquisition period.
                                                     % They are equally
                                                     % distributed over the
                                                     % 0 to rangeMax depth
                                                     % (in wavelengths)
% From SetUpL22_14vFlashAngles.m
% TGC.CntrlPts = [330 560 780 1010 1023 1023 1023 1023]; % [0,511,716,920,1023,1023,1023,1023];
% TGC.CntrlPts = [0,511,716,920,1023,1023,1023,1023];
TGC.CntrlPts = [1023 1023 1023 1023 1023 1023 1023 1023];

TGC(1).rangeMax = endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC); % Parameters can be adjusted later with GUI sliders


%% Receiver array object

% Receive(1).sampleMode = 'NS200BW'; % RF signal is Nyquist Sampled at x% BandWidth of the transducer's center frequency.
                                   %  ‘BS100BW’, ‘BS67BW’ and ‘BS50BW’.
                                   %  Recommended to use the lowest option,
                                   %  along with the corresponding
                                   %  inputFilter to lower DMA transfer
                                   %  bandwidth needs and improve SNR.
% Receive(1).sampleMode = 'BS100BW'; % idk but the L22_14v BW is 8 MHz, which is 50% the center freq. So I doubled that.                                   

maxAcqLength = ceil(sqrt(endDepth^2 + (numElements*Trans.spacing)^2)); % account for the longest distance an echo could travel
% nspa = 4*(2*maxAcqLength); % see p104. Receive does 2x of maxAcqLength

% % change 9/12/24
% spw = 3.6; % with 13.8889 MHz
% nspa = spw*(2*(maxAcqLength - startDepth)); % see p104. Receive does 2x of (maxAcqLength - startDepth)
% nspa = 128 * ceil(nspa/128); % # samples per acquisition
% maxAcqLength_adjusted = nspa / spw / 2;

Receive = repmat(struct('Apod', ones(1, Trans.numelements), ...
                        'startDepth', startDepth, ...
                        'endDepth', maxAcqLength + startDepth, ... % change 9/12/24, was previously just maxAcqLength
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ... 
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0, ... % 1 to call Media func above %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        'LowPassCoef', [], ...
                        'InputFilter', []), 1, numSupFrames * numSubFrames *na);

k = 1;
for nsupf = 1:numSupFrames
    j = 0;
    for nsubf = 1:numSubFrames
        
        Receive(k).callMediaFunc = movePointsOrNot;
        for n = 1:na
            j = j + 1;
            
            Receive(k).framenum = nsupf;
            Receive(k).acqNum = j;
            k = k + 1;
        end
    end
end

% Lowpass and Bandpass digital filters (hardware only, no simulation).
% Empty values cause VSX to program a default set of values, according to
% transducer frequency. See decimation/downsampling.
% Receive(1).LowPassCoef = []; % cutoff is frequencies higher than 2*f_center. See tutorial p19 for possible inputs
% Receive(1).InputFilter = []; % See tutorial p20 for possible inputs. Good to eliminate any DC components

% Note: after running VSX, Receive will have startSample and endSample
% fields. Since a column can have data from multiple acquisitions
% (different acqNums), it helps to have this index

% showGeometry

%% PData structure (Pixel Data --> image reconstruction range)
% For 2D scans and slices of 3D scans, it's always a rectangular area at a
% fixed location in the transducer coord system


% PData.PDelta = [Trans.spacing, 0, 0.086647727272727]; % Spacing between pixels in x, y, z, in wavelengths
PData.PDelta = [Trans.spacing, 0, Trans.spacing/2]; % Spacing between pixels in x, y, z, in wavelengths

PData.Coord = 'rectangular'; % rectangular coords, could change to polar or spherical
% Set PData array dimensions --> # of rows, columns, sections (planes
% parallel to the xy plane)
% For a 3D scan, rows - y axis, columns - x axis, sections - z axis
% PData.Size(1) = ceil((endDepth - startDepth)./PData.PDelta(3)); % # rows
PData.Size(1) = floor((endDepth - startDepth)./PData.PDelta(3)); % # rows
% (old)
% PData.Size(1) = ceil((Receive(1).endDepth - Receive(1).startDepth)./PData.PDelta(3)); % # rows
PData.Size(2) = ceil(numElements.*Trans.spacing./PData.PDelta(1)); % # cols
PData.Size(3) = 1; % depth, is 1 unit deep for a 2D image

% Define the location (x, y, z) of the upper left corner of the array
% half_probe_dist = (numElements-1)./2.*Trans.spacing;
% PData.Origin = [-half_probe_dist, half_probe_dist, startDepth];
PData.Origin = [-(numElements-1)./2.*Trans.spacing, 0, startDepth];
% Upper left corner if you look aligned with positive z

% Set a local region to view/use for processing
% PData.Region(1) = struct('Shape',struct('Name','PData'));
% 
% PData.Region(2).Shape = struct('Name', 'Slice', 'Orientation', 'xz', ...
%                             'oPAIntersect', PData.Origin(2) - (numElements-1).*Trans.spacing./2); % out of Plane Axis Intersection
% PData.Region(3).Shape = struct('Name', 'Slice', 'Orientation', 'yz', ...
%                             'oPAIntersect', PData.Origin(1) + (numElements-1).*Trans.spacing./2);
% PData.Region(4).Shape = struct('Name', 'Slice', 'Orientation', 'xy', ...
%                             'oPAIntersect', Media.MP(3)); % currently set to the plane intersecting the only scatter point
% 
% %     'Position', [0, 0, 10], ...
% %                      'width', PData.Size(2), 'height', PData.Size(1)./2);
%                     % Position is relative to the global coords
% PData.Region = computeRegions(PData);

% Display window

% xz
Resource.DisplayWindow(1).Title = 'xz plane';
Resource.DisplayWindow(1).pdelta = 0.3; % pixel spacing (in wavelengths) on the display window, for all dimensions
llx = 100; % lower left corner x on screen
xmult = 200;
lly = 150; % lower left corner y
Resource.DisplayWindow(1).Position = [llx, lly, ...
                                      ceil(PData.Size(2).* PData.PDelta(1) ./ Resource.DisplayWindow(1).pdelta), ... % width (x)
                                      ceil(PData.Size(1).* PData.PDelta(3) ./ Resource.DisplayWindow(1).pdelta)]; % height (z)
Resource.DisplayWindow(1).ReferencePt = [PData.Origin(1), 0, PData.Origin(3)]; % Display Window location wrt transducer coords
Resource.DisplayWindow(1).AxesUnits = 'wavelengths'; % can change to mm
Resource.DisplayWindow(1).Colormap = gray(256);
% Resource.DisplayWindow(1).Orientation = 'xz';
Resource.DisplayWindow(1).numFrames = numSupFrames; % Define buffer size for a history of displayed frames


%% Allocate storage space for RF data acquisition
% RcvBuffer dimensions are (samples, channels, frames)
% RcvBuffer is accessible in Matlab as a cell array

% If not set by the user, VSX will set the effective ADC sample rate to 4x
% the transducer center frequency (system max seems to be 250 MHz). There
% are preset values for sampling rate, and the system will choose the
% closest one above whatever the user sets.

% RcvBuffer dimensions: (samples, channels, frames)!!!!!!!!!!!!!!!!!!!!!!!

Resource.RcvBuffer(1).datatype = 'int16'; % 16 bit signed integers are the only supported datatype

%%%% from Nikunj's SetUpCustomIntegratedRecon.m code
if strcmp(Receive(1).sampleMode,'custom')
    error('No handling of condition for custom Receive sampling. Refer to VsUpdate line 712 to implement');
else
    fs = 4*Trans.frequency;
    samplesPerWave = 4;
end

% if statement included to match verasonics automatic extension to
% multiples of 128 samples
nSmpls = 2*(maxAcqLength - startDepth) * samplesPerWave; % maxAcqLength is the Receive(1).endDepth
% nSmpls = 2*(Receive(1).endDepth - Receive(1).startDepth) * samplesPerWave;
if abs(round(nSmpls/128) - nSmpls/128) < .01
    numRcvSamples = 128*round(nSmpls/128);
else
    numRcvSamples = 128*ceil(nSmpls/128);
end

startSample = (0:(na-1))*numRcvSamples + 1;
endSample = startSample + numRcvSamples - 1;
%%%%

% spw = 3.6765; % samples per wave, it isn't always exactly 4... check p107
% nspa = spw*(2*(Receive(1).endDepth - Receive(1).startDepth));
% nspa = 128 * ceil(nspa/128); % # samples per acquisition
% maxAcqLength_adjusted = nspa / spw / 2;
Resource.RcvBuffer(1).rowsPerFrame = numRcvSamples * na * numSubFrames;
maxAcqLength_adjusted = numRcvSamples / samplesPerWave / 2;

for lss = 1:length(startSample)
    Receive(lss).startSample = startSample(lss);
    Receive(lss).endSample = endSample(lss);    
%     Receive(lss).decimSampleRate = samplesPerWave * Trans.frequency;
%     Receive(lss).decimSampleRate = 62.5;

end

Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels; % Usually 1:1 to # of receive channels available in the system. Can change to 256 with the 2D probe and new connector plate.
Resource.RcvBuffer(1).numFrames = numSupFrames; % minimum # frames of RF data to acquire; RcvBuffer contains all the data needed for a whole frame, including multiple acquisition passes needed for reconstruction. Software can re-process RcvBuffer frames
Resource.Parameters.verbose = 2; % Describe errors in varying levels

%% Reconstruction
% numRegions = 1;
% 
% Resource.ImageBuffer(1).numFrames = numFrames; % Define an ImageBuffer with a # of frames
% Resource.InterBuffer(1).numFrames = numFrames;
% Resource.InterBuffer(1).pagesPerFrame = na; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Recon = struct('senscutoff', 0.6, ... % Threshold for which the reconstruction doesn't consider an element's contribution due to directivity of the element, for a certain pixel (whose echoes are at an angle to the element). Should be in radians.
% %                'pdatanum', 1, ... % Which PData structure to use
% %                'rcvBufFrame', -1, ... % Use the most recently transferred frame
% %                'IntBufDest', [1, 1], ... % idk but it's for the IQ (complex) data
% %                'ImgBufDest', [1, -1], ... % [buffer #, frame #] Auto-increment ImageBuffer for each reconstruction???? % something is [first/oldest frame, last/newest frame]
% %                'RINums', [1:2*na]); % The ReconInfo structure #(s). Each Recon must have its own unique set of ReconInfo #s
% 
% sco = 0.6;
% % sco = 0.4;
% Recon = struct('senscutoff', sco, ... % Threshold for which the reconstruction doesn't consider an element's contribution due to directivity of the element, for a certain pixel (whose echoes are at an angle to the element). Should be in radians.
%                'pdatanum', 1, ... % Which PData structure to use
%                'rcvBufFrame', -1, ... % Use the most recently transferred frame
%                'IntBufDest', [1, -1], ... % idk but it's for the IQ (complex) data
%                'ImgBufDest', [1, -1], ... % [buffer #, frame #] Auto-increment ImageBuffer for each reconstruction???? % something is [first/oldest frame, last/newest frame]
%                'RINums', [1:na]); % The ReconInfo structure #(s). Each Recon must have its own unique set of ReconInfo #s
% 
% % Recon = repmat(Recon, 1, numFrames);
% % for nf = 1:numFrames
% %     Recon(nf).IntBufDest = [1, nf];
% %     Recon(nf).ImgBufDest = [1, nf];
% % end
% 
% % The final image it shows only shows the image buffer data for the last
% % page (angle)
% ReconInfo = repmat(struct('mode', 'accumIQ', ... % reconstruct, and replace intensity data in ImageBuffer and IQ data in InterBuffer (see Table 12.4 in Tutorial)
%                    'txnum', 1, ...                 % TX structure to use
%                    'rcvnum', 1, ...                % RX structure to use
%                    'regionnum', 1), ...
%                    1, na);                % PData Region to process in
% 
% 
% for n = 1:na % need to change this and above for more than 1 frame
%     % - Set specific ReconInfo attributes.
%     % ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
%     ReconInfo(n).txnum = n;
%     ReconInfo(n).rcvnum = n;
%     ReconInfo(n).pagenum = n; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     ReconInfo(1).regionnum = 1; %1 for the whole volume, 5 for the slices
% 
% end
% ReconInfo(end).mode = 'accumIQ_replaceIntensity';
% %% Process the Reconstructed data
% % e.g., scaling and compression to make an image look good on the screen
% 
% pgainValue = 1.0; % Image processing gain
% persValue = 0;
% rejectLevel = 2;
% % rejectLevel = 300;
% compFac = 40;
% 
% % First image (xz)
% Process(2).classname = 'Image';
% Process(2).method = 'imageDisplay';             % To not overwrite orig data while processing, system uses another ImageP buffer as output.
% Process(2).Parameters = {'imgbufnum', 1, ...             % which ImageBuffer to process
%                          'framenum', -1, ...              % -1 means use last frame in ImageBuffer
%                          'pdatanum', 1, ...              % PData structure which was used in Reconstruction
%                          'pgain', pgainValue, ...
%                          'reject', rejectLevel, ...                % Make intensity values below this threshold appear as black (reduce low intensity noise)
%                          'persistMethod', 'simple', ...   % simple: Add a fraction of the previous weighted average frames' invensity values to the current one. See manual for 'dynamic' option, which is good when there is a lot of motion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                          'persistLevel', persValue, ...   % NewAvg = PL * PrevAvg + (1 - PL) * NewFrame; PL = persistLevel/100
%                          'interpMethod', '4pt', ...
%                          'grainRemoval', 'medium', ...    % low, medium, high. Remove pixels that differ significantly from their neighbors
%                          'processMethod', 'none', ...     % see manual, reduces variation in line structures detected within the filter kernel???
%                          'averageMethod', 'none', ...     % None or can do Running averages (2 or 3), can do things like spatial compounding...
%                          'compressMethod', 'log', ...     % log or power (x^a fraction) compression
%                          'compressFactor', compFac, ...        % Higher compressFactor means smaller powers for the power option (more compression), or a more rapid rise to the log curve (raise brightness of low intensity values). Not a real log bc intensities of 0 need to be mapped to 0, not -inf
%                          'mappingMethod', 'full', ...     % Portion of the colormap to use. lowerHalf and upperHalf would be used to do the combined B-mode and Doppler imaging, for example.
%                          'display', 1, ...                % 1: show processed image on screen, 0: don't but still tto useransfer the processed data to the DisplayData buffer
%                          'displayWindow', 1};             % which displayWindow 

%% External processing function call

Process(1).classname = 'External';
Process(1).method = 'saveRcvData'; % Function name
% Process(1).Parameters = {'srcbuffer', 'bufferName', ...
%                          'srcbufnum', 1, ... % # of buffer to process
%                          'srcframenum', 1, ... % starting frame #
%                          'srcsectionnum', 1, ...
%                          ' srcpagenum', 1, ...
%                          'dstbuffer', 'bufferName', ... % destination buf
%                          'dstbufnum', 1, ...
%                          'dstframenum', 1, ...
%                          'dstsectionnum', 1, ...
%                          'dstpagenum', 1};

Process(1).Parameters = {'srcbuffer', 'receive', ...
                         'srcbufnum', 1, ... % # of buffer to process
                         'dstbuffer', 'none'};

% Process(2).classname = 'External';
% Process(2).method = 'procSVD2DRegionSaveIQAcrossAngles'; % Function name
% 
% Process(2).Parameters = {'srcbuffer', 'receive', ...
%                          'srcbufnum', 1, ... % # of buffer to process
%                          'srcframenum', -1, ... % starting frame #
%                          'dstbuffer', 'none'};

%%
makeParameterStructureSmall;

%% Define sequence of Events for data acquisition

% Flow:
% 1. Transmit (TX)
% 2. Receive (Receive)
% 3. Reconstruction (Recon)
% 4. Processing (Process)
% 5. Control (SeqControl)


% SeqControl(1).argument = 150;  
scInd = 1; % sequence control index
SeqControl(scInd).command = 'timeToNextAcq'; % In us
                                         % Very useful if you are switching
                                         % the TPC (voltage) across acqs,
                                         % since it takes 800 us - 8 ms to
                                         % switch

% SeqControl(scInd).argument = 1./fps_target.*1e6; % convert fps to s, then us
timePerAcq = 1 / fps_target / (na * 2) * 1e6; % frame rate converted to acq time step (us)

timePerAcqLimits = [10, 4190000];
if timePerAcq < timePerAcqLimits(1)
    warning('Acquisition time too short, setting to minimum of 10 us')
    SeqControl(scInd).argument = timePerAcqLimits(1); 
elseif timePerAcq > timePerAcqLimits(2)
    warning('Acquisition time too long, setting to maximum of 4190000 us')
    SeqControl(scInd).argument = timePerAcqLimits(2);
else
    SeqControl(scInd).argument = timePerAcq;
end

scInd = scInd + 1;

SeqControl(scInd).command = 'returnToMatlab';

scInd = scInd + 1;

SeqControl(scInd).command = 'jump'; % jump to
SeqControl(scInd).argument = 1;     % first event


% SeqControl(3).command = 'transferToHost'; % probably define this before Event stuff for longer loops
% Transfer data to host, needed for hardware but not simulation, which writes directly to RcvBuffer
% NEED A UNIQUE SEQCONTROL FOR EACH TRANSFERTOHOST COMMAND!!!!!!!!!!!!!!!!!!

n = 0;

for nsupf = 1:numSupFrames
    for nsubf = 1:numSubFrames
        for a = 1:na % go through all the angles for each frame
            n = n + 1;
            Event(n).info = 'Transmit and receive';
            Event(n).tx = a; % Use ath TX structure
            Event(n).rcv = (nsupf - 1) * na * numSubFrames + (nsubf - 1)*na + a; % Receive structure
            Event(n).recon = 0; % 0 means no reconstruction
            Event(n).process = 0; % 0 means no processing
            Event(n).seqControl = 1;
    
        end
    end
    
    % Transfer all the acquisitions for one superframe 
    n = n + 1;

    Event(n).info = 'Transfer data'; % want a transferToHost after all the acquisitions for one frame
    Event(n).tx = 0; 
    Event(n).rcv = 0; 
    Event(n).recon = 0;
    Event(n).process = 0; 
    scInd = scInd + 1;

    SeqControl(scInd).command = 'transferToHost';
    
    scInd = scInd + 1;
    
    SeqControl(scInd).command = 'waitForTransferComplete';
    SeqControl(scInd).argument = scInd - 1;

    Event(n).seqControl = [scInd, scInd-1];

%     Event(n).info = ['Frame ' num2str(nf) ': Doing IQ data thing................'];
%     Event(n).tx = 0; 
%     Event(n).rcv = 0; 
%     Event(n).recon = 0; 
%     Event(n).process = 2; 
%     Event(n).seqControl = 0; 
% 
%     i = i + 1;

%     n = n + 1;
% 
%     Event(n).info = 'Save data - ext proc func';
%     Event(n).tx = 0; 
%     Event(n).rcv = 0; 
%     Event(n).recon = 0;
%     Event(n).process = 1; 
%     Event(n).seqControl = 0; 


%     i = i + 1;
% 
%     Event(n).info = ['Frame ' num2str(nf) ': Reconstruction'];
%     Event(n).tx = 0; 
%     Event(n).rcv = 0; 
%     Event(n).recon = 1;  %%
%     Event(n).process = 0; 
%     Event(n).seqControl = 0;
% 
%     i = i + 1;
% 
%     Event(n).info = ['Frame ' num2str(nf) ': Processing xz'];
%     Event(n).tx = 0; 
%     Event(n).rcv = 0; 
%     Event(n).recon = 0; 
%     Event(n).process = 1; 
%     Event(n).seqControl = 3; 
    

%     i = i + 1;
% 
%     Event(n).info = 'Jump';
%     Event(n).tx = 0; 
%     Event(n).rcv = 0; 
%     Event(n).recon = 0;
%     Event(n).process = 0; 
%     Event(n).seqControl = 3; 

end

%% Save all the data/structures to a .mat file.
currentDir = cd; currentDir = regexp(currentDir, filesep, 'split');
filename = 'L22-14v_Allen.mat';

save(fullfile(currentDir{1:find(contains(currentDir,"Vantage"),1)})+"\MatFiles\"+filename);


%% 
if runVSX
    disp("running VSX")
    VSX
%     evalin('base', 'setup_L22_14v_ctsacq1')
end

% save([savepath 'params.mat'], 'angles', 'startDepth', 'startDepthMM', 'endDepth', 'endDepthMM', 'Event', 'fps_target', 'maxAcqLength_adjusted', 'maxAngle', 'Media', 'na', 'nf', 'Receive', 'Resource', 'SeqControl', 'TGC', 'Trans', 'TW', 'TX', 'wl', 'numElements', '-v7.3')

%% Save parameters in a structure P
makeParameterStructure;
% save([savepath, 'params.mat'], 'angles', 'startDepth', 'startDepthMM', 'endDepth', 'endDepthMM', 'Event', 'fps_target', 'maxAcqLength_adjusted', 'maxAngle', 'Media', 'na', 'nf', 'Receive', 'Resource', 'SeqControl', 'TGC', 'Trans', 'TW', 'TX', 'wl', 'numElements', '-v7.3')
savefast([savepath, 'params.mat'], 'P')

saveRcvData(RcvData{1})
