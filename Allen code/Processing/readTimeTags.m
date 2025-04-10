% Check Time tags - the ones that are saved per superframe/buffer
% adapting code from Nikunj

% clearvars
% close all

%% Load parameters

datapath = "G:\Allen\Data\04-01-2025 functional acq testing\";

% RFName = "RF-5-11-2000-500-1-1.mat";
RFName = "RF-5-11-1000-180-1-1.mat";
% RFName = "RF-5-11-500-500-1-1.mat";
% RFName = "RF-5-11-1000-500-1-1.mat";
% RFName = "RF-5-11-400-500-1-1.mat";

%% Define Load path and save path for RF and IQ respectively
RFcount = countFiles(RFName,datapath);

%% Main Loop
fileInfo = strsplit(RFName,'-');
timeTags = zeros(RFcount, 1);

startTimeTag = load(datapath + 'startTimeTag.mat').('timetag'); % Time tag at the very start of the experiment

for iFile = 1:RFcount
% for iFile = 1:1
    
    iFileInfo = fileInfo;
    iFileInfo{end} = [num2str(iFile), '.mat'];
    iFileName = strjoin(iFileInfo, '-');
    
    % Load IQ Data
    disp(['Loading data: ', iFileName]);
%     RFData = load(fullfile(datapath, iFileName),'RcvData').('RcvData');
    iTimetag = load(fullfile(datapath, iFileName),'timetag').('timetag');
    disp('Data loaded!');

    timeTags(iFile) = calculateElapsedTime(iTimetag, startTimeTag);

end
% timeTags = [0; timeTags]; % add the zero timepoint

%% Plot Time tags
figure; plot(timeTags, '-o'); title('TimeTags')

timeTagDiff = diff(timeTags);
figure; plot(timeTagDiff); title('difference')

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

function [elapsedTime] = calculateElapsedTime(timetag, startTimeTag)
    elapsedTime = timetag - startTimeTag;
%     elapsedTime.Format = 'hh:mm:ss.SSS';
    elapsedTime = seconds(elapsedTime); % Convert from the duration object into seconds
    
end

