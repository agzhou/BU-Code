%
% File name: PostProcessVermon1024_8MHz_FlashAngles_Light_superframe.m
%            Post process the data collected by the
%            SetUpVermon1024_8MHz_FlashAngles_light_superframe script
%
% This script is an example on how to use the Verasonics software in simulate mode 2
% to reconstruct all of the collected RF data offline.  This Setup script performs
% the following steps:
%  1. Loads in the .mat used to collect the RF data.
%  2. Adjusts Resource structure and Receive structure to match the an efficient
%     reconstruction loop to beamform all of the RF Data.
%  3. Creates 2 external processes.  One for loading in the RF Data and one for
%     saving the reconstructed ImgData.
%  4. Replaces the Event sequence with an efficient one that performs all of the
%     needed file loading, reconstruction, and saving of the data.
%
%
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage NXT Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% Copyright © 2013-2025 Verasonics, Inc.

% Copyright (C) 2001-2025, Verasonics, Inc.
% All worldwide rights and remedies under all intellectual
% property laws and industrial property laws are reserved.

clear all


load('MatFiles/Vermon1024_8MHz-FlashAngles-Light_superframe.mat','Receive','Trans','Resource','TX','TW','Recon','ReconInfo','TGC','P')
nSubFrames = P.nSubFrames;
nAcqs = P.numAcqPerFrame;
Receive = Receive(1:nAcqs*nSubFrames);
Resource.Parameters.simulateMode = 2;
Resource.RcvBuffer.rowsPerFrame = Resource.RcvBuffer.rowsPerFrame/nSubFrames;
Resource.RcvBuffer.numFrames = nSubFrames;
Resource.ImageBuffer.numFrames = 1;
Receive = Receive(1:nAcqs);
pgain = 88;


%PData(1).PDelta = [1,1,1]*1;
PData(1).PDelta = [Trans.spacing, Trans.spacing, 0.5];

if ~exist('extent','var'), extent = max(max(Trans.ElementPosWL(:,1)),max(Trans.ElementPosWL(:,2))); end
zApex = -extent/tan(P.viewAngle);

PData(1).Size(1) = ceil(2.0*(P.endDepth-zApex)*tan(P.viewAngle)/PData(1).PDelta(2));  if mod(PData(1).Size(1),2)==0, PData(1).Size(1) =  PData(1).Size(1)+1; end
PData(1).Size(2) = PData(1).Size(1);
PData(1).Size(3) = ceil((P.endDepth)/PData(1).PDelta(3));
PData(1).Origin = [-((PData(1).Size(2)-1)/2)*PData(1).PDelta(1), ((PData(1).Size(1)-1)/2)*PData(1).PDelta(1), 0];

PData(1).Region(1) = struct('Shape',struct('Name','Pyramid','Position',[0,0,zApex],'angle',P.viewAngle,'z1',P.startDepth,'z2',P.endDepth));
PData(1).Region(2) = struct('Shape',struct('Name','Slice','Orientation','yz','oPAIntersect',PData.Origin(1)+PData.PDelta(2)*(PData.Size(2)-1)/2,'andWithPrev',1));
PData(1).Region(3) = struct('Shape',struct('Name','Slice','Orientation','xz','oPAIntersect',PData.Origin(2)-PData.PDelta(1)*(PData.Size(1)-1)/2,'andWithPrev',1));
%PData(1).Region(4) = struct('Shape',struct('Name','Slice','Orientation','xy','oPAIntersect',XYDisp,'andWithPrev',1));

% [PData(1).Region] = computeRegions(PData(1));
%
% PData.Region(5).PixelsLA = unique([PData.Region(2).PixelsLA; PData.Region(3).PixelsLA; PData.Region(4).PixelsLA]);
% PData.Region(5).Shape.Name = 'Custom';
% PData.Region(5).numPixels = length(PData.Region(5).PixelsLA);
[PData(1).Region] = computeRegions(PData(1));

PData.Region(4).PixelsLA = unique([PData.Region(2).PixelsLA; PData.Region(3).PixelsLA]);
PData.Region(4).Shape.Name = 'Custom';
PData.Region(4).numPixels = length(PData.Region(4).PixelsLA);
%% if we want to reconstruct the full 3d volume

for a = 1:length(ReconInfo)
    ReconInfo(a).regionnum = 1;  %complete volume
    %ReconInfo(a).regionnum = 4;  %orthogonal planes
    %ReconInfo(a).regionnum = 2;  %YZ Plane only
    %ReconInfo(a).regionnum = 3;  %XZ Plane only
end
%%%%%%%%


Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).Title = '3D FlashAngles-Spiral-SRAC Image - XZ plane';
Resource.DisplayWindow(1).pdelta = 0.25;
Resource.DisplayWindow(1).Position = [0,580, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).Orientation = 'xz';
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0.0,0.0];
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).AxesUnits = 'wavelengths';
Resource.DisplayWindow(1).numFrames = 20;

Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).Title = '3D FlashAngles-Spiral-SRAC Image - YZ plane';
Resource.DisplayWindow(2).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(2).Position = [430,580, ...
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
Resource.DisplayWindow(2).Orientation = 'yz';
Resource.DisplayWindow(2).ReferencePt = [0,-PData(1).Origin(2),0];
Resource.DisplayWindow(2).Colormap = gray(256);
Resource.DisplayWindow(2).AxesUnits = 'wavelengths';
Resource.DisplayWindow(2).numFrames = 20;

%%
% Specify Process structure arrays.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...
    'framenum',-1,...
    'pdatanum',1,...
    'srcData','intensity3D',...
    'pgain', pgain,...
    'persistMethod','none',...
    'persistLevel',30,...
    'interpMethod','4pt',...
    'compressMethod','log',...
    'compressFactor',55,...
    'mappingMethod','full',...
    'display',1,...
    'displayWindow',1};

Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',1,...
    'framenum',-1,...
    'pdatanum',1,...
    'srcData','intensity3D',...
    'pgain', pgain,...
    'persistMethod','none',...
    'persistLevel',30,...
    'interpMethod','4pt',...
    'compressMethod','log',...
    'compressFactor',55,...
    'mappingMethod','full',...
    'display',1,...
    'displayWindow',2};

Process(3).classname = 'External';
Process(3).method = 'loadRF';
Process(3).Parameters = { ...
    'srcbuffer','none',... % name of buffer to export.
    'dstbuffer','receive',... % name of buffer to export.
    'dstbufnum',    1,    ... % buffer number of buffer to export
    };

Process(4).classname = 'External';
Process(4).method = 'saveImageDataP';
Process(4).Parameters = { ...
    'srcbuffer','image',... % name of buffer to export.
    'srcbufnum',    1,    ... % buffer number of buffer to export
    'srcframenum',  1,   ... % process the last frame in a RcvBuffer
    'dstbuffer','none',... % name of buffer to export.
    };

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'returnToMatlab';

frameRateFactor = 1; % Factor for converting sequenceRate to frameRate.
nsc = length(SeqControl)+1; % Count of SeqControl objects.

n = 1; % Event structure arrays.
saveEvents = [];

Event(n).info = 'Load RF SuperFrame';
Event(n).tx = 0;    % No TX structure.
Event(n).rcv = 0;   % No Rcv structure of frame.
Event(n).recon = 0; % Reconstruction of the slices in real time or the whole volume depending on ReconRegion
Event(n).process = 3;   % no processing
Event(n).seqControl = 0;
n = n+1;

for a = 1:nSubFrames
    Event(n).info = ['3D Reconstruct Frame ' num2str(j)];
    Event(n).tx = 0;    % No TX structure.
    Event(n).rcv = 0;   % No Rcv structure of frame.
    Event(n).recon = 1; % Reconstruction of the slices in real time or the whole volume depending on ReconRegion
    Event(n).process = 0;   % no processing
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['Process XZ - Frame ' num2str(j)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['Process YZ - Frame ' num2str(j)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'save ImgDataP + make movie';
    Event(n).tx = 0;    % No TX structure.
    Event(n).rcv = 0;   % No Rcv structure of frame.
    Event(n).recon = 0; % Reconstruction of the slices in real time or the whole volume depending on ReconRegion
    Event(n).process = 4;   % no processing
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Return to Matlab';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 2;
    n=n+1;
end
Event(n).info = 'Jump';
Event(n).tx = 0;         % no TX structure.
Event(n).rcv = 0;        % no Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    %
Event(n).seqControl = 1; %

%% User specified UI Control Elements

EF(1).Function = vsv.seq.function.ExFunctionDef('loadRF', @loadRF);
EF(2).Function = vsv.seq.function.ExFunctionDef('saveImageDataP', @saveImageDataP);


% Save all the structures to a .mat file.
filename = 'MatFiles/M3dV-FlashAngles-Light_superframe_post.mat';
save(filename);


%% Function for loading in RF Data and for offline reconstruction
function RFDataReshaped = loadRF()
persistent rfDims dirListSorted fileListSorted fileIndex subFrameIndex init fid startTime

startFileIndex = 1;
nSubFrames = evalin('base','P.nSubFrames');
nSampsPerVol = evalin('base','Receive(P.numAcqPerFrame).endSample');

if isempty(init)
    startTime = tic();
    rfDims(1) = evalin('base','Resource.RcvBuffer.rowsPerFrame');
    rfDims(2) = evalin('base','Resource.RcvBuffer.colsPerFrame');
    rfDims(3) = nSubFrames;
    if isunix()
        dataDir = '/media/verasonics/WD';
        fileList1 = dir(sprintf('%s%i/*.rf',dataDir,1));
        fileList2 = dir(sprintf('%s%i/*.rf',dataDir,2));
        fileList3 = dir(sprintf('%s%i/*.rf',dataDir,3));
        fileList4 = dir(sprintf('%s%i/*.rf',dataDir,4));
    elseif ispc()
        fileList1 = dir('D:\*.rf');
        fileList2 = dir('E:\*.rf');
        fileList3 = dir('F:\*.rf');
        fileList4 = dir('G:\*.rf');
    end
    fileList = [fileList1; fileList2; fileList3; fileList4];
    [fileListSorted, sortIndex] = sort({fileList.name});
    dirList = {fileList.folder};
    dirListSorted = dirList(sortIndex);
    fileIndex = startFileIndex;
    subFrameIndex = 1;
    init = 1;
end

fprintf("Processing superframe %d subframe %d ", fileIndex, subFrameIndex);
tic()
filepath = [dirListSorted{fileIndex} '/' fileListSorted{fileIndex}];
fprintf("loading superFrame File: %s    ", filepath);
fid = fopen(filepath, 'rb');
RFData = fread(fid, [rfDims(1)*nSubFrames, rfDims(2)], '*int16');
RFDataReshaped = zeros(rfDims,'int16'); %Initialize the buffer
for a = 1:nSubFrames
    RFDataReshaped(1:nSampsPerVol,:,a) = RFData((a-1)*nSampsPerVol + (1:nSampsPerVol),:);
end
readTime = toc();
fprintf('file read/reshape speed: %g GB/s\n', numel(RFData)*2/1e9/readTime);
fclose(fid);

fileIndex = fileIndex + 1;

% Progress reporting
elapsedTime = seconds(toc(startTime));
elapsedTime.Format = 'hh:mm:ss';
complete = fileIndex/length(dirListSorted);
totalTime = elapsedTime / complete;
totalTime.Format = 'hh:mm:ss';
remainingTime = (totalTime - elapsedTime);
remainingTime.Format = 'hh:mm:ss';
fprintf("percent complete = %g    Elapsed Time: %s  Total Time: %s  Time Remaining: %s\n",round(complete*100,2), elapsedTime, totalTime, remainingTime)

end


%% Function for saving the processed image data into a new directory
function saveImageDataP(ImgDataP)

persistent fileIndex
if isempty(fileIndex)
    fileIndex = 1;
end

dataDir = '/media/verasonics/WD1/imgDataP/';
if ~exist(dataDir, 'dir')
    mkdir(dataDir)
end

fname = sprintf('%simgDataP%06d.mat', dataDir, fileIndex);
fprintf("saving: %s\n", fname);
save(fname,"ImgDataP")

fileIndex = fileIndex + 1;
end
