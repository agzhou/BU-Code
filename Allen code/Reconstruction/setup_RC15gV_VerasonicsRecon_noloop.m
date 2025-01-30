
%% Description
% Verasonics reconstruction of RF data with the RC15gV probe that loops
% over all the files in a folder

% Output of IQ data for all subframes
% C-R and R-C pairs of TX-RX

% This version uses the same Recon structure (assumes the same Receive
% parameters) for each subframe

% Last updated on 01/29/2024 and working for an input of one superframe

%% TO DO
% make into a function
% add adjustable pixel spacing

%% Specify system parameters
clear

cd 'C:\Users\BOAS-US\Desktop\Vantage-4.9.5-2409181500'
% cd 'G:\My Drive\Verasonics files\Vantage-4.9.2-2308102000'

activate

% savepath = "G:\Allen\Data\12-04-2024 RC15gV saving tests\trial 1 exp\";
% savepath = char(savepath);
% mkdir(savepath)

% Uncomment/initialize this variable to some arbitrary value if you want
% the script to autoquit runAcq

% Mcr_AutoScriptTest = 1;

%% Load parameters and RcvData, then reshape the RcvData
% [paramsFilename, paramsPath] = uigetfile;
% [RcvDataFilename, RcvDataPath] = uigetfile;
% 
% load([paramsPath, paramsFilename])
% load([RcvDataPath, RcvDataFilename])

%% Data loading test
clearvars

datapath = 'G:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\';

load([datapath, 'params.mat']) % load acquisition parameters

filenameStructure = ['RF-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

filenum = 25;                                           % TEMPORARY FOR TESTING !!!!!!!!!!!
saveAllAngles = 0; % choose if you want to save the matrix with pages for each angle or not
load([datapath, filenameStructure, num2str(filenum)]);

%% Put RcvData into a cell array for VSX
r = RcvData;
clear RcvData;

RcvData{1} = r;
%% Assign variables and structures

pair = 2; % The R-C and C-R pair of acquisitions per angle
numChannels = P.Resource.Parameters.numTransmit;

% assignFromParameterStructure;
assignStructVars(P);                                                        % Assign parameters from P structure into the workspace

maxAcqLength = maxAcqLength_adjusted;
Trans_acq = Trans;
TX_acq = TX;
TW_acq = TW;
Resource_acq = Resource;
Receive_acq = Receive;
Receive = Receive(1:numFramesPerBuffer * na * pair);
clear Event Process Recon ReconInfo SeqControl Trans TW TX Resource



%% Resource, define system parameters
Resource.Parameters.numTransmit = numChannels;                              % number of transmit channels
Resource.Parameters.numRcvChannels = numChannels;                           % number of receive channels
% Resource.Parameters.connector = 1;                                        % transducer connector to use since the current plate for the 256 bit system is split into two 128 bit connectors. 1 is left and 2 is right
Resource.Parameters.speedOfSound = Resource_acq.Parameters.speedOfSound;    % speed of sound in m/s

Resource.RcvBuffer = Resource_acq.RcvBuffer(1);
Resource.RcvBuffer.rowsPerFrame = size(RcvData{1}, 1);
Resource.RcvBuffer.numFrames = numFramesPerBuffer;
Resource.RcvBuffer.lastFrame = 1; % reset the counter
Resource.Parameters.simulateMode = 2; % Enable mode 2, which processes data in the buffers
%% Define Transducer

Trans.name = Trans_acq.name; 
Trans.frequency = Trans_acq.frequency;
Trans.units = Trans_acq.units;

Trans = computeTrans(Trans); % Generate required attributes for the probe into the Trans structure; e.g., the transducer element positions

%% TW
TW(1).type = TW_acq.type;
TW(1).Parameters = TW_acq.Parameters;

%% TX
TX_fn = fieldnames(TX_acq);
TX = rmfield(TX_acq, TX_fn(6:14));

%% Receive (need to be careful about the superframe and subframe definition since ReconInfo depends on it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
% Receive = repmat(struct('Apod', zeros(1, Trans.numelements), ... 
%                         'startDepth', startDepth, ...
%                         'endDepth', maxAcqLength, ...
%                         'TGC', 1, ...
%                         'bufnum', 1, ...
%                         'framenum', 1, ...
%                         'acqNum', 1, ...
%                         'sampleMode', 'NS200BW', ...
%                         'mode', 0, ...
%                         'callMediaFunc', 0, ...
%                         'LowPassCoef', [], ...
%                         'InputFilter', []), 1, pair*P_new.numSupFrames*P_new.numSubFrames*na);
% j = 1;
% % an = 0;
% for nsupf = 1:P_new.numSupFrames
%     for nsubf = 1:P_new.numSubFrames
%         % Move points after all the acquisitions for one frame
% %         Receive(j).callMediaFunc = movePointsOrNot;
%     %     Receive(j).mode = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         an = 0;
% 
%         for n = 1:na
%             an = an + 1;
%             Receive(j).framenum = nsubf;
%             Receive(j).acqNum = an;
%             Receive(j).Apod(Trans.numelements/2 + 1 : end) = ones(1, Trans.numelements/2);
%             j = j + 1;
%         end
%     
%         for n = 1:na
%             an = an + 1;
%             Receive(j).framenum = nsubf;
%             Receive(j).acqNum = an;
%             Receive(j).Apod(1:Trans.numelements/2) = ones(1, Trans.numelements/2);
%             j = j + 1;
%         end
%         
%     end
% end
%% PData structure (Pixel Data --> image reconstruction range)
% For 2D scans and slices of 3D scans, it's always a rectangular area at a
% fixed location in the transducer coord system

numElements = Trans.numelements./2; % the structure gives # row elements + # column elements

PData.PDelta = [Trans.spacing, Trans.spacing, 0.5]; % Spacing between pixels in x, y, z, in wavelengths

PData.Coord = 'rectangular'; % rectangular coords, could change to polar or spherical
% Set PData array dimensions --> # of rows, columns, sections (planes
% parallel to the xy plane)
% For a 3D scan, rows - y axis, columns - x axis, sections - z axis
PData.Size(1) = ceil(numElements.*Trans.spacing./PData.PDelta(2)); % # rows
PData.Size(2) = ceil(numElements.*Trans.spacing./PData.PDelta(1)); % # cols
PData.Size(3) = ceil((endDepth - startDepth)./PData.PDelta(3)); % sections

% Define the location (x, y, z) of the upper left corner of the array
half_probe_dist = (numElements-1)./2.*Trans.spacing;
PData.Origin = [-half_probe_dist, half_probe_dist, startDepth];
% PData.Origin = [-half_probe_dist, -half_probe_dist, startDepth];

% Upper left corner if you look aligned with positive z

% Set a local region to view/use for processing
PData.Region(1) = struct('Shape',struct('Name','PData'));

PData.Region(2).Shape = struct('Name', 'Slice', 'Orientation', 'xz', ...
                            'oPAIntersect', PData.Origin(2) - (numElements-1).*Trans.spacing./2); % out of Plane Axis Intersection
PData.Region(3).Shape = struct('Name', 'Slice', 'Orientation', 'yz', ...
                            'oPAIntersect', PData.Origin(1) + (numElements-1).*Trans.spacing./2);
PData.Region(4).Shape = struct('Name', 'Slice', 'Orientation', 'xy', ...
                            'oPAIntersect', Media.MP(3)); % currently set to the plane intersecting the only scatter point

PData.Region = computeRegions(PData);

%% Reconstruction
numRegions = 3;

Resource.ImageBuffer(1).numFrames = numFramesPerBuffer; % Define an ImageBuffer with a # of frames
Resource.InterBuffer(1).numFrames = numFramesPerBuffer;

if saveAllAngles
    Resource.InterBuffer(1).pagesPerFrame = na * pair;
end
% Recon = struct('senscutoff', 0.6, ... % Threshold for which the reconstruction doesn't consider an element's contribution due to directivity of the element, for a certain pixel (whose echoes are at an angle to the element). Should be in radians.
%                'pdatanum', 1, ... % Which PData structure to use
%                'rcvBufFrame', -1, ... % Use the most recently transferred frame
%                'IntBufDest', [1, 1], ... % idk but it's for the IQ (complex) data
%                'ImgBufDest', [1, -1], ... % [buffer #, frame #] Auto-increment ImageBuffer for each reconstruction???? % something is [first/oldest frame, last/newest frame]
%                'RINums', [1:2*na]); % The ReconInfo structure #(s). Each Recon must have its own unique set of ReconInfo #s

sco = 0.6; %%%%

Recon = repmat(struct('senscutoff', sco, ... % Threshold for which the reconstruction doesn't consider an element's contribution due to directivity of the element, for a certain pixel (whose echoes are at an angle to the element). Should be in radians.
               'pdatanum', 1, ... % Which PData structure to use
               'rcvBufFrame', -1, ... % Use the most recently transferred frame
               'IntBufDest', [1, -1], ... % IQ (complex) data, [buffer, frame], -1 means use the next available frame as output
               'ImgBufDest', [1, -1], ... % [buffer #, frame #]
               'RINums', [1:pair*na]), 1, 1); % The ReconInfo structure #(s). Each Recon must have its own unique set of ReconInfo #s

ReconInfo = repmat(struct('mode', 'accumIQ', ... % reconstruct, and replace intensity data in ImageBuffer and IQ data in InterBuffer (see Table 12.4 in Tutorial)
                   'txnum', 1, ...                 % TX structure to use
                   'rcvnum', 1, ...                % RX structure to use
                   'regionnum', 1), 1, pair*na);                % PData Region to process in

rii = 0; % recon info index

% Modify ReconInfo
for n = 1:pair*na
    rii = rii + 1;
    
    ReconInfo(rii).txnum = n;
    ReconInfo(rii).rcvnum = rii;
    if saveAllAngles
        ReconInfo(rii).pagenum = n;
        ReconInfo(rii).mode = 'replaceIQ'; % replace IQ data
    end
%     ReconInfo(1).regionnum = 1; %1 for the whole volume, 5 for the slices

end

if ~saveAllAngles
    ReconInfo(1).mode = 'replaceIQ'; % replace IQ in the buffer for each new frame processed
    ReconInfo(end).mode = 'accumIQ_replaceIntensity'; % at the last acquisition, update the ImgData
end
%% Process the Reconstructed data
% e.g., scaling and compression to make an image look good on the screen

% pgainValue = 1.0; % Image processing gain
% persValue = 0;
% rejectLevel = 2;
% % rejectLevel = 300;
% compFac = 40;

% First image (xz)
% Process(1).classname = 'Image';
% Process(1).method = 'imageDisplay';             % To not overwrite orig data while processing, system uses another ImageP buffer as output.
% Process(1).Parameters = {'imgbufnum', 1, ...             % which ImageBuffer to process
%                          'framenum', -1, ...              % -1 means use last frame in ImageBuffer
%                          'pdatanum', 1, ...              % PData structure which was used in Reconstruction
%                          'srcData', 'intensity3D', ...
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
% 
% % Second image (xy)
% Process(2).classname = 'Image';
% Process(2).method = 'imageDisplay';             % To not overwrite orig data while processing, system uses another ImageP buffer as output.
% Process(2).Parameters = {'imgbufnum', 1, ...             % which ImageBuffer to process
%                          'framenum', -1, ...              % -1 means use last frame in ImageBuffer
%                          'pdatanum', 1, ...              % PData structure which was used in Reconstruction
%                          'srcData', 'intensity3D', ...
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
%                          'displayWindow', 2};             % which displayWindow 
% 
% % Third image (yz)
% Process(3).classname = 'Image';
% Process(3).method = 'imageDisplay';             % To not overwrite orig data while processing, system uses another ImageP buffer as output.
% Process(3).Parameters = {'imgbufnum', 1, ...             % which ImageBuffer to process
%                          'framenum', -1, ...              % -1 means use last frame in ImageBuffer
%                          'pdatanum', 1, ...              % PData structure which was used in Reconstruction
%                          'srcData', 'intensity3D', ...
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
%                          'displayWindow', 3};             % which displayWindow 

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
%                          'srcframenum', -1, ... % last frame transferred

% Process(2).classname = 'External';
% Process(2).method = 'ShowTimeTag'; % Function name
% Process(2).Parameters = {'srcbuffer','receive',... % name of buffer to process.
%                          'srcbufnum',1,...
%                          'srcframenum',-1, ...
%                          'dstbuffer','none'};

%% New Event structure

% Flow:
% 1. Transmit (TX)
% 2. Receive (Receive)
% 3. Reconstruction (Recon)
% 4. Processing (Process)
% 5. Control (SeqControl)


SeqControl(1).command = 'noop'; % VSX errors if there is no SeqControl structure

n = 0;

for nf = 1:numFramesPerBuffer
    n = n + 1;

    Event(n).info = ['Frame ' num2str(nf) ': Reconstruction'];
    Event(n).tx = 0; 
    Event(n).rcv = 0; 
    Event(n).recon = 1;
    Event(n).process = 0; 
    Event(n).seqControl = 0; 
end

% 
% n = n + 1;
% 
% Event(n).info = 'Save data - ext proc func';
% Event(n).tx = 0; 
% Event(n).rcv = 0; 
% Event(n).recon = 0;
% Event(n).process = 1; 
% Event(n).seqControl = 0; 


%% Save all the data/structures to a .mat file.
currentDir = cd; currentDir = regexp(currentDir, filesep, 'split');
filename = 'RC15gV_Allen_recon.mat';

save(fullfile(currentDir{1:find(contains(currentDir,"Vantage"),1)})+"\MatFiles\"+filename);


%% Run VSX automatically and make parameter structure for RF file naming
tic

    disp("running VSX")
    VSX

toc
%% 
IQ = IData{1} + 1i .* QData{1};
figure; imagesc(abs(squeeze(IQ(40, :, :, 1, 1)))')
figure; imagesc(abs(squeeze(IQ(:, 40, :, 1, 1)))')
% saveRcvData(RcvData{1})
% ImgData_temp = ImgData{1};
% figure; imagesc(squeeze(ImgData_temp(40, :, :, 2))')
IQ = single(IQ ./ 1000);