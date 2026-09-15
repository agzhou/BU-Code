% Check Time tags

% clearvars
% close all

%% Load parameters

FilePath = uigetdir('F:\', 'Select the RF data path');
FilePath = string([FilePath, '\']);

RFPath = FilePath;

load(FilePath + "params.mat") % Load acquisition parameters

[RFName, ~, ~] = uigetfile('*.mat', 'Select the first RF file', [RFPath]);
RFName = string(RFName);

RFcount = countFiles(RFName, RFPath); % Count # of RF files in the path

%% Main Loop
fileInfo = strsplit(RFName,'-');

% Framewise time tags
% RFtimeTags = zeros(str2double(fileInfo{5}), P.Resource.Parameters.numRcvChannels, RFcount);

% Superframe-wise time tags (stacked frames)
RFtimeTags = zeros(1, P.Resource.Parameters.numRcvChannels, RFcount);

for iFile = 1:RFcount
% for iFile = 1:4

    iFileInfo = fileInfo;
    iFileInfo{end} = [num2str(iFile), '.mat'];
    iFileName = strjoin(iFileInfo, '-');
    
    % Load IQ Data
    disp(['Loading data: ', iFileName]);
    RFData = load(fullfile(RFPath, iFileName),'RcvData').('RcvData');
    disp('Data loaded!');

    RFtimeTags(:, :, iFile) = readTimeTags(RFData);

end

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

genSliderV2(RFtimeTags) % Plot the time tags per channel, for each superframe
figure; imagesc(squeeze(RFtimeTags(:, :, 1))); colorbar; xlabel('Channel index'); ylabel('Frame')
figure; imagesc(diff(squeeze(RFtimeTags), 1, 2).')

%% Save the (RF) frame timing data
% save(FilePath + "RFTimeTags.mat", 'RFTimeTags_raw', 'RFTimeTags', "RFTimeTags_diff", 'RFcount')
save(FilePath + "RFTimeTagsPerChannel.mat", 'RFtimeTags', 'RFcount')

%% Test
% getTimeStamp(double(RFData(1:2,1,12)))/4e4

%% Helper Functions
function [fileCount] = countFiles(fileName,filePath)

    fileInfo = strsplit(fileName, '-');

    % This keeps everything except the last numeric component
    prefix = strjoin(fileInfo(1:end-1), '-'); 
    
    % Construct the search pattern
    searchPattern = fullfile(filePath, prefix + "-*.mat"); % Wildcard for different numbers
    
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

