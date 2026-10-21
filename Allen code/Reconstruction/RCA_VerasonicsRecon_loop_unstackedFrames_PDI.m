
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
% Add some regexp thing to automatically get the # of raw data files in the folder

%% Add the Reconstruction folder to path
clearvars

codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
AllenReconCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Reconstruction\");
addpath(genpath(AllenReconCodePath))

%% Activate the Verasonics folder

cd 'C:\Users\BOAS-US\Desktop\Vantage-5.0.0-p1'
% cd 'C:\Users\BOAS-US\Desktop\Vantage-4.9.5-2409181500'
% cd 'C:\Users\BOAS-US\Desktop\Vantage-4.9.7-2505271400'
% cd '\\ad\eng\users\a\g\agzhou\Desktop\Vantage-4.9.7-2505271400'
% cd 'C:\Users\agzhou\Vantage-4.9.7-2505271400'
% addpath 'C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Allen code\Reconstruction'

% cd 'G:\My Drive\Verasonics files\Vantage-4.9.2-2308102000'
% cd 'C:\Users\agzhou\Vantage-4.9.7-2505271400'
% cd 'C:\Users\agzhou\Vantage-4.9.5-2409181500'
% addpath '\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Reconstruction';

activate

%% Load parameters, create save path, choose some options for recon

% Mcr_datapath = 'G:\Allen\Data\03-17-2025 AZ02 ULM\RC15gV\run 2 right eye\';
Mcr_datapath = uigetdir('G:\Allen\Data\', 'Select the raw data path');
Mcr_datapath = [Mcr_datapath, '\'];
if ~exist('P', 'var')
    load([Mcr_datapath, 'params.mat']) % load acquisition parameters
end

Mcr_savepath = uigetdir(Mcr_datapath, 'Select the folder to save reconstructed data to');
Mcr_savepath = [Mcr_savepath, '\'];

% Prompt for parameter user input
parameterPrompt = {'Start file number', 'End file number'};
parameterDefaults = {'1', ''};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

Mcr_startFile = str2double(parameterUserInput{1}); % File to start reconstructing from
Mcr_endFile = str2double(parameterUserInput{2});   % File to stop reconstructing on

Mcr_filenameStructure = ['RF-', num2str(round(P.maxAngle)), '-', num2str(P.na), '-', num2str(round(P.frameRate)), '-', num2str(P.numFramesPerBuffer), '-1-'];
Mcr_IQfilenameStructure = ['IQ-', num2str(round(P.maxAngle)), '-', num2str(P.na), '-', num2str(round(P.frameRate)), '-', num2str(P.numFramesPerBuffer), '-1-'];

%%%%%%%%%%
% Output mode (prompted). The options are mutually exclusive.
%   Mcr_saveAllAngles: save the matrix with pages for each angle
%   Mcr_saveRCCR:      save separate coherently-compounded C-R and R-C volumes (for XDoppler) instead of the OPW sum.
%                      Saves IQ_CR and IQ_RC, each [y, x, z, frame], to files with the 'IQ-RCCR-' prefix
% (Mcr_ prefix so VSX doesn't clear these)
Mcr_outputOptions = {'OPW sum (single IQ volume)', 'Separate C-R and R-C volumes (XDoppler)', 'All angles (pages for each angle)'};
Mcr_outputChoice = questdlg('Select the IQ output type:', 'Output options', ...
                            Mcr_outputOptions{1}, Mcr_outputOptions{2}, Mcr_outputOptions{3}, Mcr_outputOptions{1});
if isempty(Mcr_outputChoice)
    error('No output type selected.')
end
Mcr_saveAllAngles = strcmp(Mcr_outputChoice, Mcr_outputOptions{3});
Mcr_saveRCCR      = strcmp(Mcr_outputChoice, Mcr_outputOptions{2});
disp(strcat("Output type: ", Mcr_outputChoice))

if Mcr_saveRCCR
    Mcr_IQfilenameStructure = ['IQ-RCCR-', Mcr_IQfilenameStructure(4:end)]; % strip the 'IQ-' prefix and replace
end

%% Get the updated parameter structure for unstacked data
[P_unstacked] = updateParams_unstackedFrames(P);

%% autorun VSX flag
% Uncomment/initialize this variable to some arbitrary value if you want
% the script to autoquit runAcq

Mcr_AutoScriptTest = 1;
% Mcr_GuiHide = 1;

%% Assign variables and structures

pair = 2; % The R-C and C-R pair of acquisitions per angle
numChannels = P_unstacked.Resource.Parameters.numTransmit;

% assignFromParameterStructure;
assignStructVars(P_unstacked);                                                        % Assign parameters from P structure into the workspace

maxAcqLength = maxAcqLength_adjusted;
Trans_acq = Trans;
TX_acq = TX;
TW_acq = TW;
% Resource = P_unstacked.Resource;
% Receive_acq = Receive;
% Receive = Receive(1:numFramesPerBuffer * na * pair);
% Receive = P_unstacked.Receive;
clear Event Process Recon ReconInfo SeqControl Trans TW TX



%% Resource, define system parameters
% Resource.Parameters.numTransmit = numChannels;                              % number of transmit channels
% Resource.Parameters.numRcvChannels = numChannels;                           % number of receive channels
% % Resource.Parameters.connector = 1;                                        % transducer connector to use since the current plate for the 256 bit system is split into two 128 bit connectors. 1 is left and 2 is right
% Resource.Parameters.speedOfSound = Resource_acq.Parameters.speedOfSound;    % speed of sound in m/s
% 
% Resource.RcvBuffer = Resource_acq.RcvBuffer(1);
% 
% % load the first one to get this size
% load([Mcr_datapath, Mcr_filenameStructure, '1']);
% rpf = size(RcvData, 1);
% clear RcvData
% Resource.RcvBuffer.rowsPerFrame = rpf;
% 
% Resource.RcvBuffer.numFrames = numFramesPerBuffer;
Resource.RcvBuffer.lastFrame = 1; % reset the counter
Resource.Parameters.simulateMode = 2; % Enable mode 2, which processes data in the buffers

% Resource.Parameters.verbose = 2; % Describe errors in varying levels
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

%% PData structure (Pixel Data --> image reconstruction range)
% For 2D scans and slices of 3D scans, it's always a rectangular area at a
% fixed location in the transducer coord system

numElements = Trans.numelements./2; % the structure gives # row elements + # column elements

PData.PDelta = [Trans.spacing, Trans.spacing, 0.5]; % Spacing between pixels in x, y, z, in wavelengths
% PData.PDelta = [0.5, 0.5, 0.5]; % Spacing between pixels in x, y, z, in wavelengths

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

if Mcr_saveAllAngles
    Resource.InterBuffer(1).pagesPerFrame = na * pair;
elseif Mcr_saveRCCR
    Resource.InterBuffer(1).pagesPerFrame = pair; % page 1: C-R (TX column, RX row), page 2: R-C (TX row, RX column)
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
    if Mcr_saveAllAngles
        ReconInfo(rii).pagenum = n;
        ReconInfo(rii).mode = 'replaceIQ'; % replace IQ data
    elseif Mcr_saveRCCR
        % TX 1:na are C-R (acquired first in each frame), TX na+1:2na are R-C.
        % Each group is coherently compounded over its na angles into its own page.
        ReconInfo(rii).pagenum = 1 + (n > na);
    end
%     ReconInfo(1).regionnum = 1; %1 for the whole volume, 5 for the slices

end

if Mcr_saveRCCR
    ReconInfo(1).mode = 'replaceIQ';        % start of C-R page
    ReconInfo(na + 1).mode = 'replaceIQ';   % start of R-C page
    ReconInfo(end).mode = 'accumIQ_replaceIntensity'; % at the last acquisition, update the ImgData
elseif ~Mcr_saveAllAngles
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
Mcr_filename = 'RC15gV_Allen_recon_unstackedFrames.mat';

save(fullfile(currentDir{1:find(contains(currentDir,"Vantage"),1)})+"\MatFiles\"+Mcr_filename, "Event", "SeqControl", "ReconInfo", "Recon", "Resource", "Media", "PData", "Receive", "TGC", "TPC", "Trans", "TW", "TX");


%% Run VSX automatically and reconstruct/save each file
if ~exist('Mcr_P', 'var')
    Mcr_P = P; % VSX_auto will clear the P, so set it to not be cleared because unstackFrames needs it
end

save([Mcr_savepath, 'PData'], 'PData') % Save the PData structure

for Mcr_filenum = Mcr_startFile:Mcr_endFile
% for Mcr_filenum = Mcr_endFile:-1:Mcr_startFile
% for Mcr_filenum = 142:Mcr_endFile
% for Mcr_filenum = [37, 110, 111]
    tic

    load([Mcr_datapath, Mcr_filenameStructure, num2str(Mcr_filenum)], 'RcvData');
    disp(strcat("Raw data file ", num2str(Mcr_filenum), " loaded."))
    
    % Put RcvData into a cell array for VSX
    r = unstackFrames(RcvData, Mcr_P);
    clearvars RcvData;
    
    RcvData{1} = r;
    clear r;

    filename = Mcr_filename; % VSX clears variables without the Mcr_ prefix, so redefine "filename" so VSX can autorun
    Mcr_AutoScriptTest = 1;  % VSX also specially clears this, so redefine it

%     disp("running VSX_auto")
    VSX_auto % this is in the Verasonics folder
    VsClose  % close the GUI window. runAcq stops automatically after one loop.

    pause(10)

    disp(strcat("IQ file ", num2str(Mcr_filenum), " reconstructed."))
%     IQ = IData{1} + 1i .* QData{1};                                         % Merge the I and Q into one variable
    IData = IData{1};
    QData = QData{1};

    if Mcr_saveRCCR
        % IData/QData are [y, x, z, page, frame]. Page 1 = C-R, page 2 = R-C.
        % Index the pages before complex() so the full 5D array isn't duplicated.
        IQ_CR = squeeze(complex(IData(:, :, :, 1, :), QData(:, :, :, 1, :)));
        IQ_RC = squeeze(complex(IData(:, :, :, 2, :), QData(:, :, :, 2, :)));
        clear IData QData
        % save([Mcr_savepath, Mcr_IQfilenameStructure, num2str(Mcr_filenum)], 'IQ_CR', 'IQ_RC', '-v7.3', '-nocompression')
        % disp(strcat("IQ (C-R and R-C) file ", num2str(Mcr_filenum), " saved."))

        % ==== Calculate PDI ==== %
        % SVD filter CR and RC IQ volume stacks separately
        [xpo, ypo, zpo, nfo] = size(IQ_CR); % Original sizes
        % figure; imagesc(squeeze(max(abs(IQ_CR(:, :, :, 1)), [], 1)).')
        
        xrange = 1:xpo;
        yrange = 1:ypo;
        % zrange = 30:130;
        zrange = 10:110;
        % zrange = 1:zpo;
        % zrange = zstart:zend;
        voxelRange = {xrange, yrange, zrange};
        xp = length(xrange); yp = length(yrange); zp = length(zrange); nf = size(IQ_CR, 4);
            
        % Mask the IQ to some region (according to the function input)
        xrange = voxelRange{1};
        yrange = voxelRange{2};
        zrange = voxelRange{3};
        IQ_CR = IQ_CR(xrange, yrange, zrange, :);
        IQ_RC = IQ_RC(xrange, yrange, zrange, :);
        
        % [brainMask, maskInfo] = createBrainVoxelMask(sum(abs(IQ), 4), 'MarginVoxels', 0);
        % IQ = IQ .* brainMask;   % zero out non-brain voxels before the SVD sees them
        
        % SVD decluttering
        % [xp, yp, zp, nf] = size(IQ);
        sv_threshold_lower = 30; % Lower singular value threshold
        sv_threshold_upper = nf;
        IQf_RC = SVD_CF(IQ_RC, sv_threshold_lower, sv_threshold_upper);
        IQf_CR = SVD_CF(IQ_CR, sv_threshold_lower, sv_threshold_upper);

        IQ_OPW = IQ_CR + IQ_RC;
        clearvars IQ_RC IQ_CR
        
        % IQf_OPW = IQf_RC + IQf_CR;
        IQf_OPW = SVD_CF(IQ_OPW, sv_threshold_lower, sv_threshold_upper);
        PDI_OPW = mean(abs(IQf_OPW).^2, 4);
        
        % XD = IQf_RC .* conj(IQf_CR);
        % PDI_XD = mean(abs(XD).^2, 4);
        % PDI_XD = 2.*real(mean(XD, 4));
        % PDI_XD = abs(mean(IQf_RC .* conj(IQf_CR) + IQf_CR .* conj(IQf_RC), 4));
        PDI_XD = abs(mean(IQf_RC .* conj(IQf_CR), 4)); % Not exactly the estimator from the paper, but close
        save([Mcr_savepath, 'PDI-', num2str(Mcr_filenum)], 'PDI_OPW', 'PDI_XD', '-v7.3')
    else
%     IQ = squeeze(IData + 1i .* QData);
        IQ = squeeze(complex(IData, QData));

%     savefast([Mcr_savepath, Mcr_IQfilenameStructure, num2str(Mcr_filenum)], 'IData', 'QData')
%     save([Mcr_savepath, Mcr_IQfilenameStructure, num2str(Mcr_filenum)], 'IData', 'QData', '-v7.3')
        save([Mcr_savepath, Mcr_IQfilenameStructure, num2str(Mcr_filenum)], 'IQ', '-v7.3', '-nocompression')
        disp(strcat("IQ file ", num2str(Mcr_filenum), " saved."))
        warning('Does not save PDI for OPW only')

    end

    toc

    % ixc = calcIXC_simple(IQ);
    % %     figure; plot(abs(ixc)); xlabel('Frame'); ylabel('|Cross correlation of images|')
    % save([Mcr_savepath, 'ixc-', num2str(Mcr_filenum)], 'ixc', '-v7.3', '-nocompression')

    % clearvars IQ IQ_CR IQ_RC IData QData RcvData ImgData ImgDataP
%     clearvars RcvData ImgData ImgDataP
    
    pause(10) % Pause for safety of inter-superframe memory issues

end

save([Mcr_savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'voxelRange', 'xp', 'yp', "zp", "nf");

%% Make a plot of all the ixcs
% ixc_allfiles = [];
% % for filenum = Mcr_startFile:Mcr_endFile
% for test_fn = 146:346
%     load([Mcr_savepath, 'ixc-', num2str(test_fn), '.mat'])
%     ixc_allfiles = cat(2, ixc_allfiles, ixc);
% end
% 
% figure; plot(abs(ixc_allfiles)); xlabel("Frame"); ylabel('|Cross correlation of images|')

%% saving speed test
% tic
% test = IData{1};
% savefast([Mcr_savepath, Mcr_IQfilenameStructure, num2str(Mcr_filenum)], 'test')
% save([Mcr_savepath, Mcr_IQfilenameStructure, num2str(Mcr_filenum)], 'IQ', '-v7.3')
% toc
%% XDoppler axial velocity estimate from a saveRCCR file (Leroy et al., Ultrasonics 2026, Eq. 8)
% Run the clutter filter (e.g., SVD) on IQ_CR and IQ_RC separately first.
% Slow-time order is CR(1), RC(1), CR(2), RC(2), ...
% The estimator is symmetric in CR/RC, but the interleaving order matters.
%
% load([Mcr_savepath, Mcr_IQfilenameStructure, '1'])   % IQ_CR, IQ_RC: [y, x, z, frame]
% nf = size(IQ_CR, 4);
% C1 = sum(IQ_CR(:, :, :, 1:nf)   .* conj(IQ_RC(:, :, :, 1:nf)),   4); % CR(i) * conj(RC(i + 0.5))
% C2 = sum(IQ_RC(:, :, :, 1:nf-1) .* conj(IQ_CR(:, :, :, 2:nf)),   4); % RC(i + 0.5) * conj(CR(i + 1))
% % (optionally smooth C1, C2 over a ~3x3x3 voxel window before taking the angle)
% T = 1 / P.frameRate;         % time between consecutive CR (or RC) volumes. This is 2*na*T_PRF only if frames are back-to-back.
% c = P.Resource.Parameters.speedOfSound; f0 = P.Trans.frequency * 1e6; % check these field names for your P
% vz = -c ./ (4*pi*f0*T) .* (angle(C1) + angle(C2));   % sign follows the paper's convention; check against your Kasai code
% pd = sqrt(abs(C1) ./ nf);    % XDoppler power Doppler: sqrt of the RC/CR correlation product (Eq. 2)

%% random stuff and plotting
% volumeViewer(squeeze(abs(IQ(:, :, :, 2))))
% IQ = IData{1} + 1i .* QData{1};
% figure; imagesc(abs(squeeze(IQ(40, :, :, 1, 1)))')
% figure; imagesc(abs(squeeze(IQ(:, 40, :, 1, 1)))')
% saveRcvData(RcvData{1})
% ImgData_temp = ImgData{1};
% figure; imagesc(squeeze(ImgData_temp(40, :, :, 2))')
% IQ = single(IQ ./ 1000);