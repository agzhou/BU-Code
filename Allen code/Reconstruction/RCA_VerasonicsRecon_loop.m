
%% Description
% Verasonics reconstruction of RF data with the RC15gV probe that loops
% over all the files in a folder

% Output of IQ data for all subframes
% C-R and R-C pairs of TX-RX

% This version uses the same Recon structure (assumes the same Receive
% parameters) for each subframe

% Last updated on 01/29/2024 and working for an input of multiple
% superframes (goes through a folder of raw data)

%% TO DO
% make into a function
% add adjustable pixel spacing
% Add some regexp thing to automatically get the # of raw data files in the folder
%% Activate the Verasonics folder
clearvars

addpath(fullfile(cd))
cd 'C:\Users\BOAS-US\Desktop\Vantage-4.9.7-2505271400' % Lab PC (bme-boas-19)
addpath 'C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Allen code\Processing'
% cd 'C:\Users\agzhou\Vantage-4.9.7-2505271400' % Office PC (bme-boas-27)

activate

%% Load parameters and RcvData, then reshape the RcvData
% [paramsFilename, paramsPath] = uigetfile;
% [RcvDataFilename, RcvDataPath] = uigetfile;
% 
% load([paramsPath, paramsFilename])
% load([RcvDataPath, RcvDataFilename])

%% Load parameters, create save path, choose some options for recon
clearvars

% Mcr_datapath = 'G:\Allen\Data\03-17-2025 AZ02 ULM\RC15gV\run 2 right eye\';
Mcr_datapath = uigetdir('G:\Allen\Data\', 'Select the raw data path');
Mcr_datapath = [Mcr_datapath, '\'];
load([Mcr_datapath, 'params.mat']) % load acquisition parameters

% Mcr_savepath = [Mcr_datapath, 'IQ Data - Verasonics Recon\'];
% Mcr_savepath = 'F:\Allen\Data\03-17-2025 AZ02 ULM\RC15gV\run 2 right eye\';
Mcr_savepath = uigetdir(Mcr_datapath, 'Select the folder to save reconstructed data to');
Mcr_savepath = [Mcr_savepath, '\'];
% Mcr_savepath = ['K:\Allen data backup\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\', 'IQ Data with pages - Verasonics Recon\']; 
% mkdir(Mcr_savepath)

% Prompt for parameter user input
parameterPrompt = {'Start file number', 'End file number'};
parameterDefaults = {'1', ''};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

Mcr_startFile = str2double(parameterUserInput{1}); % File to start reconstructing from
Mcr_endFile = str2double(parameterUserInput{2});   % File to stop reconstructing on

Mcr_filenameStructure = ['RF-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];
Mcr_IQfilenameStructure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

%%%%%%%%%%
saveAllAngles = 0; % choose if you want to save the matrix with pages for each angle or not

%% autorun VSX flag
% Uncomment/initialize this variable to some arbitrary value if you want
% the script to autoquit runAcq

Mcr_AutoScriptTest = 1;
% Mcr_GuiHide = 1;

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

% load the first one to get this size
load([Mcr_datapath, Mcr_filenameStructure, '1']);
rpf = size(RcvData, 1);
clear RcvData
Resource.RcvBuffer.rowsPerFrame = rpf;

Resource.RcvBuffer.numFrames = numFramesPerBuffer;
Resource.RcvBuffer.lastFrame = 1; % reset the counter
Resource.Parameters.simulateMode = 2; % Enable mode 2, which processes data in the buffers

Resource.Parameters.verbose = 2; % Describe errors in varying levels
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

%% New Event structure

% Flow:
% 1. Transmit (TX)
% 2. Receive (Receive)
% 3. Reconstruction (Recon)
% 4. Processing (Process)
% 5. Control (SeqControl)


SeqControl(1).command = 'noop'; % VSX errors if there is no SeqControl structure

n = 0;
Event = struct('info', {}, 'tx', {}, 'rcv', {}, 'recon', {}, 'process', {}, 'SeqControl', {});

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
Mcr_filename = 'RC15gV_Allen_recon.mat';

save(fullfile(currentDir{1:find(contains(currentDir,"Vantage"),1)})+"\MatFiles\"+Mcr_filename, "Event", "SeqControl", "ReconInfo", "Recon", "Resource", "Media", "PData", "Receive", "TGC", "TPC", "Trans", "TW", "TX");


%% Run VSX automatically and make parameter structure for RF file naming


for Mcr_filenum = Mcr_startFile:Mcr_endFile
    tic

    load([Mcr_datapath, Mcr_filenameStructure, num2str(Mcr_filenum)]);
    disp(strcat("Raw data file ", num2str(Mcr_filenum), " loaded."))
    
    % Put RcvData into a cell array for VSX
    r = RcvData;
    clearvars RcvData;
    
    RcvData{1} = r;
    clear r;

    filename = Mcr_filename; % VSX clears variables without the Mcr_ prefix, so redefine "filename" so VSX can autorun
    Mcr_AutoScriptTest = 1;  % VSX also specially clears this, so redefine it

%     disp("running VSX_auto")
    VSX_auto % this is in the Verasonics folder
%     VSX
    VsClose  % close the GUI window. runAcq stops automatically after one loop.

    pause(10)

    disp(strcat("IQ file ", num2str(Mcr_filenum), " reconstructed."))
%     IQ = squeeze(IData{1} + 1i .* QData{1});                                         % Merge the I and Q into one variable
    IData = IData{1};
    QData = QData{1};
%     IQ = squeeze(IData + 1i .* QData);
    IQ = squeeze(complex(IData, QData));

    % savefast([Mcr_savepath, Mcr_IQfilenameStructure, num2str(Mcr_filenum)], 'IData', 'QData')
    save([Mcr_savepath, Mcr_IQfilenameStructure, num2str(Mcr_filenum)], 'IQ', '-v7.3', '-nocompression')
%     save([Mcr_savepath, Mcr_IQfilenameStructure, num2str(Mcr_filenum)], 'IQ', '-v7.3')
    disp(strcat("IQ file ", num2str(Mcr_filenum), " saved."))
    
    toc

    ixc = calcIXC_simple(IQ);
    %     figure; plot(abs(ixc)); xlabel('Frame'); ylabel('|Cross correlation of images|')
    save([Mcr_savepath, 'ixc-', num2str(Mcr_filenum)], 'ixc')

    clearvars IQ IData QData RcvData ImgData ImgDataP
%     clearvars RcvData ImgData ImgDataP
    
    pause(5) % Pause for safety of inter-superframe memory issues
end

save([Mcr_savepath, 'PData'], 'PData') % Save the PData structure

%% Make a plot of all the ixcs
ixc_allfiles = [];
for test_fn = Mcr_startFile:Mcr_endFile
% for test_fn = 146:346
    load([Mcr_savepath, 'ixc-', num2str(test_fn), '.mat'])
    ixc_allfiles = cat(2, ixc_allfiles, ixc);
end

figure; plot(abs(ixc_allfiles)); xlabel("Frame"); ylabel('|Cross correlation of images|')

%% saving speed test
% tic
% test = IData{1};
% savefast([Mcr_savepath, Mcr_IQfilenameStructure, num2str(Mcr_filenum)], 'test')
% save([Mcr_savepath, Mcr_IQfilenameStructure, num2str(Mcr_filenum)], 'IQ', '-v7.3')
% toc
%% random stuff and plotting
% volumeViewer(squeeze(abs(IQ(:, :, :, 2))))
% IQ = IData{1} + 1i .* QData{1};
% figure; imagesc(abs(squeeze(IQ(40, :, :, 1, 1)))')
% figure; imagesc(abs(squeeze(IQ(:, 40, :, 1, 1)))')
% saveRcvData(RcvData{1})
% ImgData_temp = ImgData{1};
% figure; imagesc(squeeze(ImgData_temp(40, :, :, 2))')
% IQ = single(IQ ./ 1000);