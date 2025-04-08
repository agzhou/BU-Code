% Check Time tags - the ones that are saved per superframe/buffer
% adapting code from Nikunj

% clearvars
% close all

%% Load parameters

FilePath = "G:\Allen\Data\04-01-2025 functional acq testing\";
% FilePath = "G:\Nikunj\test1\";

RFPath = FilePath;
RFName = "RF-5-11-500-100-1-1.mat";
% RFName = "RF-5-11-100-100-1-1.mat";

%% Define Load path and save path for RF and IQ respectively
RFcount = countFiles(RFName,RFPath);


%% Main Loop
fileInfo = strsplit(RFName,'-');
timeTags = zeros(RFcount, 1);
% for iFile = 1:RFcount
for iFile = 1:1
    

    iFileInfo = fileInfo;
    iFileInfo{end} = [num2str(iFile), '.mat'];
    iFileName = strjoin(iFileInfo, '-');
    
    % Load IQ Data
    disp(['Loading data: ', iFileName]);
    RFData = load(fullfile(RFPath, iFileName),'RcvData').('RcvData');
    iTimetag = load(fullfile(RFPath, iFileName),'timetag').('timetag');
    disp('Data loaded!');

%     timeTags(iFile) = convertTimeTags(RFData);

end

%% Plot Time tags
chk = timeTags(101:end);
figure; plot(chk)
title('TimeTags')
chk2 = chk - chk(1);
figure
plot(chk2)
title('TimeTags from zero')

chk3 = diff(chk2);
figure; plot(chk3); title('difference')
chk4 = find(chk3 > 0.1)

%% Test
getTimeStamp(double(RFData(1:2,1,12)))/4e4

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

function [timeTags] = convertTimeTags(RFData)
    
    timeTags = zeros(size(RFData,3),1);
    for frmCount = 1:size(RFData,3)
    
        timeStamp = getTimeStamp(double(RFData(1:2,1,frmCount)));
        % the 32 bit time tag counter increments every 25 usec, so we have to scale
        % by 25 * 1e-6 to convert to a value in seconds

        timeTags(frmCount) = timeStamp/4e4;
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

