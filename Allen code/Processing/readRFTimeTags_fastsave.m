% Check Time tags

clearvars
% close all

%% Load parameters

FilePath = uigetdir('D:\', 'Select the RF data path');
FilePath = string([FilePath, '\']);

RFPath = FilePath;

load(FilePath + "params.mat") % Load acquisition parameters

[RFName, ~, ~] = uigetfile('*.rf', 'Select the first RF file', [RFPath]);
RFName = string(RFName);

RFcount = countFiles(RFName, RFPath); % Count # of RF files in the path

%% Main Loop
fileInfo = strsplit(RFName,'-');

temp = strsplit(fileInfo(end), '.');
startFileInd = str2double(temp(1));

% Framewise time tags
% RFtimeTags = zeros(str2double(fileInfo{5}), P.Resource.Parameters.numRcvChannels, RFcount);

% Superframe-wise time tags (stacked frames)
% RFtimeTags = zeros(1, P.Resource.Parameters.numRcvChannels, RFcount);
RFtimeTags = zeros(P.numFramesPerBuffer, P.Resource.Parameters.numRcvChannels, RFcount);

% Get the start and end sample indices for the RF, for each subframe
startSamples = zeros(P.numFramesPerBuffer, 1);
endSamples = zeros(P.numFramesPerBuffer, 1);
% warning('Per-frame code for stacked RF timetags only works for linear array right now; need to do P.na*2 for RCA...')
for iFrame = 1:P.numFramesPerBuffer
    if isequal(P.Trans.name, 'L22-14v')
        startSamples(iFrame) = P.Receive((iFrame-1)*P.na + 1).startSample;
        endSamples(iFrame) = P.Receive((iFrame-1)*P.na + 1).endSample;
    elseif isequal(P.Trans.name, 'RC15gV')
        startSamples(iFrame) = P.Receive((iFrame-1)*P.na*2 + 1).startSample;
        endSamples(iFrame) = P.Receive((iFrame-1)*P.na*2 + 1).endSample;
    end
end
%%
for iFile = startFileInd:RFcount + startFileInd - 1
% for iFile = 1:4

    iFileInfo = fileInfo;
    % iFileInfo{end} = [sprintf( '%05d', iFile - 1), '.rf']; % fastsave indexes from 0 and zero pads to 5 digits
    iFileInfo{end} = [sprintf( '%05d', iFile), '.rf']; % fastsave indexes from 0 and zero pads to 5 digits
    iFileName = strjoin(iFileInfo, '-');
    
    % Load the .rf binary file
    disp(['Loading data: ', iFileName]);
    fid = fopen(FilePath+iFileName, 'rb'); 
    RFData = fread(fid, [P.Resource.RcvBuffer(1).rowsPerFrame, P.Resource.RcvBuffer(1).colsPerFrame], '*int16');
    fclose(fid);
    disp('Data loaded!');
    
    for iFrame = 1:P.numFramesPerBuffer
        RFtimeTags(iFrame, :, iFile - startFileInd + 1) = readTimeTags(RFData(startSamples(iFrame):endSamples(iFrame), :));
    end
end

RFTimeTags_allSFStacked = reshape(permute(RFtimeTags, [2, 1, 3]), [P.Resource.Parameters.numRcvChannels, P.numFramesPerBuffer*RFcount]); % Save a version of the RF timetags where all the frames are stacked

%% Plot Time tags
% RFTimeTags_raw = timeTags(1:end);
% figure; plot(RFTimeTags_raw)
% title('Raw TimeTags')
% RFTimeTags = RFTimeTags_raw - RFTimeTags_raw(1);
% figure
% plot(RFTimeTags)
% title('TimeTags from zero')
% RFTimeTags_diff = diff(RFTimeTags_raw);
% figure; plot(RFTimeTags_diff); title('difference')
% % chk4 = find(chk3 > 0.1)



% genSliderV2(RFtimeTags) % Plot the time tags per channel, for each superframe
% figure; imagesc(squeeze(RFtimeTags(:, :, 1))); colorbar; xlabel('Channel index'); ylabel('Frame')
% figure; imagesc(diff(squeeze(RFtimeTags), 1, 2).')
figure; imagesc(RFTimeTags_allSFStacked); xlabel('Frame number'); ylabel('Channel number'); cb = colorbar; cb.Label.String = 'Time tag [s]';
figure; plot(RFTimeTags_allSFStacked(1, :)); xlabel('Frame number'); ylabel('Time tag for channel 1 [s]')

%% Save the (RF) frame timing data
% save(FilePath + "RFTimeTags.mat", 'RFTimeTags_raw', 'RFTimeTags', "RFTimeTags_diff", 'RFcount')
save(FilePath + "RFTimeTagsPerChannel.mat", 'RFtimeTags', 'RFTimeTags_allSFStacked', 'RFcount')

%% Test
% getTimeStamp(double(RFData(1:2,1,12)))/4e4

%% Helper Functions
function [fileCount] = countFiles(fileName,filePath)

    fileInfo = strsplit(fileName, '-');

    % This keeps everything except the last numeric component
    prefix = strjoin(fileInfo(1:end-1), '-'); 
    
    % Construct the search pattern
    searchPattern = fullfile(filePath, prefix + "-*.rf"); % Wildcard for different numbers
    
    % Get a list of matching files
    fileList = dir(searchPattern);
    
    % Count the number of matching files
    fileCount = numel(fileList);


end

% Output the RF time tags, for each channel
function [timeTags] = readTimeTags(RFData)
    
    timeTags = zeros(size(RFData, 3), size(RFData, 2)); % # frames x # channels timetag matrix
    for frmCount = 1:size(RFData, 3)
        for channelInd = 1:size(RFData, 2)
            timeStamp = getTimeStamp(double(RFData(1:2, channelInd, frmCount)));
            % the 32 bit time tag counter increments every 25 usec, so we have to scale
            % by 25 * 1e-6 to convert to a value in seconds
    
            timeTags(frmCount, channelInd) = timeStamp/4e4;
        end
    end
end

function [tStamp] = getTimeStamp(W)

    % get time tag from first two samples
    % time tag is 32 bit unsigned interger value, with 16 LS bits in sample 1
    % and 16 MS bits in sample 2.  Note RDatain is in signed INT16 format so must
    % convert to double in unsigned format before scaling and adding
    for i=1:2
        if W(i) < 0
            % translate 2's complement negative values to their unsigned integer
            % equivalents
            W(i) = W(i) + 65536;
        end
    end
    tStamp = W(1) + 65536 * W(2);

end

