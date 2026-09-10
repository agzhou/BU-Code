%
% File name: PostProcessL22_14vXFlashAngles_fastsave.m
%            Post process the data collected by the
%            SetUpL22_14vXFlashAngles_fastSave script
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

clear all


%load('MatFiles/L22-14vXFlashAngles_fastsave.mat','Receive','Trans','Resource','TX','TW','Recon','ReconInfo','TGC','P')
load('MatFiles/L22-14vXFlashAngles_fastsave.mat')
nSubFrames = P.nSubFrames;
nAcqs = P.numAcqPerFrame;
Receive = Receive(1:nAcqs*nSubFrames);
Resource.Parameters.simulateMode = 2;
Resource.RcvBuffer.rowsPerFrame = ceil(Resource.RcvBuffer.rowsPerFrame/nSubFrames);
Resource.RcvBuffer.numFrames = nSubFrames;
Resource.ImageBuffer.numFrames = 1;
Receive = Receive(1:nAcqs);
pgain = 88;

%Resource.Parameters.initializeOnly = 1;
clear Event SeqControl
%%
Process(2).classname = 'External';
Process(2).method = 'loadRF';
Process(2).Parameters = { ...
    'srcbuffer','none',... % name of buffer to export.
    'dstbuffer','receive',... % name of buffer to export.
    'dstbufnum',    1,    ... % buffer number of buffer to export
    };

Process(3).classname = 'External';
Process(3).method = 'saveImageDataP';
Process(3).Parameters = { ...
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
Event(n).process = 2;   % no processing
Event(n).seqControl = 0;
n = n+1;

for a = 1:nSubFrames
    Event(n).info = ['Reconstruct Frame ' num2str(j)];
    Event(n).tx = 0;    % No TX structure.
    Event(n).rcv = 0;   % No Rcv structure of frame.
    Event(n).recon = 1; % Reconstruction of the slices in real time or the whole volume depending on ReconRegion
    Event(n).process = 0;   % no processing
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['Process Frame ' num2str(j)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;


    Event(n).info = 'save ImgDataP + make movie';
    Event(n).tx = 0;    % No TX structure.
    Event(n).rcv = 0;   % No Rcv structure of frame.
    Event(n).recon = 0; % Reconstruction of the slices in real time or the whole volume depending on ReconRegion
    Event(n).process = 3;   % no processing
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
nSampsPerFrame = evalin('base','Receive(P.numAcqPerFrame).endSample');

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
    RFDataReshaped(1:nSampsPerFrame,:,a) = RFData((a-1)*nSampsPerFrame + (1:nSampsPerFrame),:);
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
