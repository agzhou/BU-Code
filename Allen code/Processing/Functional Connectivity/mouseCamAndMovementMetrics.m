videoData = readMP4();

%% Get timestamp of final block
dt = (bs * (1-bo))/P.frameRate; % Time step calculated from the block size, overlap, and frame rate
% t = 0:dt:(numBlocks - 1)*dt; % Time stamps of each block

t_lastblock = numBlocks * dt; % [s]

%%
sf_to_video_delay = seconds(TD.sfStartTimetag - videoData.startTimetag);

videoCropped = videoData.frames;
us_cam_mask = videoData.timestamps > sf_to_video_delay & videoData.timestamps < t_lastblock;
videoCropped = videoCropped(:, :, us_cam_mask);
% test = resample(videoData.frames, t, dimension = 3);

%%

vw = VideoWriter()