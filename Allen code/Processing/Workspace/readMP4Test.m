%% Load the video data
[videoFilePathFN, videoFilePath] = uigetfile(['*.mp4'], 'Select the video data');
videoFilePath = [videoFilePath, videoFilePathFN];
% load(videoFilePath)

v = VideoReader(videoFilePath)

%%
% someFrames = read(v, [1, 10]);
% someFramesGS = rgb2gray(someFrames); % Convert RGB to Grayscale
tic
videoFrames = zeros(v.Height, v.Width, v.NumFrames); % Store all the video's frames in a matrix
videoTimeStamps = zeros(v.NumFrames, 1);             % Time stamps of each video frame, relative to the start of the video capture [s]
for fi = 1:v.NumFrames % Go through each video frame
% for fi = 1:10
    vfi = rgb2gray(read(v, fi)); % Video frame at frame fi
    videoTimeStamps(fi) = v.CurrentTime;
    videoFrames(:, :, fi) = vfi;
end
% clearvars fi vfi
toc

%%
videoTimeStampDiffs = diff(videoTimeStamps); % Interval between each successive video frame [s]