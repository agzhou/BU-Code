% Description:
%   Choose and read a mp4 video. Note that the timetag feature here only
%   works with how Basler saves the start of the capture in its filename.

% Output:
%   videoData: struct with the video's frames, timestamps, and start time


function [videoData] = readMP4
    %% Load the video data
    [videoFilePathFN, videoFilePath] = uigetfile(['*.mp4'], 'Select the video data');
    videoFilePath = [videoFilePath, videoFilePathFN];
    % load(videoFilePath)
    
    v = VideoReader(videoFilePath);
    
    %% Store the frames and time stamps
    % someFrames = read(v, [1, 10]);
    % someFramesGS = rgb2gray(someFrames); % Convert RGB to Grayscale
    % tic
    videoFrames = zeros(v.Height, v.Width, v.NumFrames); % Store all the video's frames in a matrix
    videoTimeStamps = zeros(v.NumFrames, 1);             % Time stamps of each video frame, relative to the start of the video capture [s]
    for fi = 1:v.NumFrames % Go through each video frame
    % for fi = 1:10
        vfi = rgb2gray(read(v, fi)); % Video frame at frame fi
        videoTimeStamps(fi) = v.CurrentTime;
        videoFrames(:, :, fi) = vfi;
    end
    % clearvars fi vfi
    % toc
    
    % Store the frames in the struct
    videoData.frames = videoFrames;
    videoData.timestamps = videoTimeStamps;

    %%
    videoTimeStampDiffs = diff(videoTimeStamps); % Interval between each successive video frame [s]
    
    %% Get the start time of the video capture from the filename
    vfp_ss = strsplit(videoFilePath, "_"); % Video file path, string split
    videoStartDateGlobalChar = vfp_ss{end - 1}; % String representing the date of the video start
    videoStartTimeGlobalChar = vfp_ss{end}; % String representing the PC's local time of the video start
    
    % Store video start info in a struct
    videoData.startYear = str2num(videoStartDateGlobalChar(1:4));
    videoData.startMonth = str2num(videoStartDateGlobalChar(5:6));
    videoData.startDay = str2num(videoStartDateGlobalChar(7:8));
    
    videoData.startHour = str2num(videoStartTimeGlobalChar(1:2));
    videoData.startMinute = str2num(videoStartTimeGlobalChar(3:4));
    videoData.startSecond = str2num(videoStartTimeGlobalChar(5:6));
    videoData.startMillisecond = str2num(videoStartTimeGlobalChar(7:9));
    
    % Create a timetag object for the video's start
    videoData.startTimetag = datetime(videoData.startYear, videoData.startMonth, videoData.startDay, videoData.startHour, videoData.startMinute, videoData.startSecond, videoData.startMillisecond, 'Format', 'dd-MMM-uuuu HH:mm:ss.SSS');
end
