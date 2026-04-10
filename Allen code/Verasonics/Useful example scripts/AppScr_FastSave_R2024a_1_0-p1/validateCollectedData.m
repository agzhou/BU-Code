%
% File name: validateCollectedData.m
%            A simple script to review the collected RF data from fastsave
%            and verify that the interval between timetags is consistent
%            and matches the expected interval.  The data is plotted and
%            any points that are offset from the red line indicate a
%            delayed or dropped frame.
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage NXT Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% Copyright © 2013-2025 Verasonics, Inc.

close all
VsClose

load('MatFiles/L22-14vXFlashAngles.mat','Resource','P')
P.numAcqPerFrame = 7;
sampsPerTR = Resource.RcvBuffer.rowsPerFrame/P.numAcqPerFrame/P.nSubFrames;

load superFrameFrameRate
expectedPeriod = 1/superFrameRate;

% --- Load RF File Names ---%
if isunix()
    fsDriveNames = {'/media/verasonics/WD1','/media/verasonics/WD2','/media/verasonics/WD3','/media/verasonics/WD4'};
elseif ispc()
    fsDriveNames = {'E:\','F:\','G:\','H:\'};
end

fileList = [];
for a = 1:length(fsDriveNames)
    dirlist = dir(sprintf('%s/*.rf', fsDriveNames{a}));
    fileList = [fileList; dirlist];
end
if (length(fileList) < 1)
    error('no files in directory')
end

% --- Load RF Data Time Tags ---%
acqindex = 1;
index = 1;
% timeTag is a variable containing the timetags at the start of each file
% timeTagAcq is a variable containing the timtag for every acquisition
for i = 1:length(fileList)
    fprintf('file:%i\n',i);
    filename = fullfile(fileList(i).folder,fileList(i).name);
    fileinfo = dir(filename);
    if (fileinfo.bytes > 0)
        fid = fopen(filename);

        % Check time tags between each super frame (timeTag)
        % Check time tags between each acquisition (timeTagAcq)
        for j = 1:(P.numAcqPerFrame*P.nSubFrames)
            try
                timeTagAcq(acqindex) = double(fread(fid,1,'int32'));

                if j == 1
                    timeTag(index) = timeTagAcq(acqindex);
                    index = index + 1;
                end
            catch
                warning('could not read timetag on this file')
                timeTagAcq(acqindex) = 0;
            end
            acqindex = acqindex + 1;
        end
        fclose(fid);
    else
        warning('empty file')
    end
end

timeTagSorted = sort(timeTag);

timeTag_sec = (timeTagSorted-timeTagSorted(1))*25e-6;% Counter increments every 25 usec, so we have to scale by 25 * 1e-6 to get seconds
timeTagPeriod = timeTag_sec(2:end) - timeTag_sec(1:end-1);
timeTagPeriod = round(timeTagPeriod,9);  %round to the nearest nanosecond to eliminate numerical noise in value

%find 1st missed frame
firstMissedFrame = find(timeTagPeriod~=expectedPeriod,1,'first');
if isempty(firstMissedFrame)
    firstMissedFrame = length(timeTagPeriod);
end


writeDuration_sec = seconds(datetime(datestr(fileList(end).datenum))-datetime(datestr(fileList(1).datenum))); %includes the extra time at the end of writing

acquisitionDuration_sec = max(timeTag_sec)-min(timeTag_sec);
acquisitionDurationNoMissed_sec =timeTag_sec(firstMissedFrame)-min(timeTag_sec);

totalSize_bytes = sum([fileList.bytes]);
totalSize_Tbytes = totalSize_bytes/1e12;
totalSizeNoMissed_bytes = sum([fileList(1:firstMissedFrame).bytes]);
totalSizeNoMissed_Tbytes = totalSizeNoMissed_bytes/1e12;

avgDataThroughput_GBpS = totalSize_bytes/acquisitionDuration_sec/1e9;
avgDataThroughputNoMissed_GBpS = totalSize_bytes/acquisitionDurationNoMissed_sec/1e9;
totalSpace = 0;
for a = 1:length(fsDriveNames)
    fileObj = java.io.File(fsDriveNames{a});
    totalSpace = totalSpace + double(fileObj.getTotalSpace())/1e9;
end
percentStorageUsed = 100*1000*totalSize_Tbytes/totalSpace;
noMissedPercentStorageUsed = 100*1000*totalSizeNoMissed_Tbytes/totalSpace;

fprintf("-------------------------------------------------------------------\n")
fprintf(" Total # files captured:  %g files\n", length(fileList));
fprintf(" Total duration of acquisition:  %d seconds\n", round(acquisitionDuration_sec));
fprintf(" Total duration of acquisition to 1st missed frame:  %d seconds\n", round(acquisitionDurationNoMissed_sec));
fprintf(" Total size of collected data:  %.2g TB\n", totalSize_Tbytes);
fprintf(" Total size of collected data to 1st missed frame:  %.2g TB\n", totalSizeNoMissed_Tbytes);
fprintf(" Percent of Storage Used:  %.2g \n", percentStorageUsed);
fprintf(" Percent of Storage Used No Missed:  %.2g \n", noMissedPercentStorageUsed);
fprintf(" Average Data Throughput:  %.2g GB/s \n", avgDataThroughput_GBpS);
fprintf(" Average Data Throughput No Missed:  %.2g GB/s \n", avgDataThroughputNoMissed_GBpS);
fprintf("-------------------------------------------------------------------\n")


%% Plot data
subplot(121)
plot(timeTag_sec(1:end-1),timeTagPeriod,'.')
xlabel('Acquisition Time (sec)')
ylabel('Period Between Super-Frames (sec)')
line([0 timeTag_sec(end)],[expectedPeriod expectedPeriod],'color',[1 0 0])
legend('super-frame period','expected period')
title('Period Between Super-Frame Acquisitions')

subplot(122)
plot((1:length(timeTagPeriod))*fileinfo.bytes/1e9,timeTagPeriod,'.')
xlabel('total data written GB')
ylabel('Period Between Super-Frames (sec)')
line([0 timeTag_sec(end)],[expectedPeriod expectedPeriod],'color',[1 0 0])
legend('super-frame period','expected period')
title('Period Between Super-Frame Acquisitions')
