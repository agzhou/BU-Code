%% Description: Plot some metrics for motion against a behavioral video

% Should already have the FC data processed using analyzeFC3D_continuous.m

%% Get timestamp of final block
dt = (bs * (1-bo))/P.frameRate; % Time step calculated from the block size, overlap, and frame rate
% t = 0:dt:(numBlocks - 1)*dt; % Time stamps of each block

t_lastblock = numBlocks * dt; % [s]

%% Load video data
videoData = readMP4();

%% Load the SV across blocks data

load([data_dirpath, 'fUSdata-', num2str(1), '.mat'], 'SVs')
SVsallBlocks = zeros([length(SVs), numBlocks]); % Matrix with the CBVi for every superframe
SVsallBlocks(:, 1) = SVs;

for bn = 1:numBlocks
%     load([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI')
    load([savepath, 'fUSdata-', num2str(bn), '.mat'], 'SVs')

    SVsallBlocks(:, bn) = SVs;
end

%% MIP video of the PDI
generateTiffStack_acrossframes(PDIallBlocks_reg .^ 0.5, [264, 160, 228] .* 50e-3, 'hot', 1:size(PDIallBlocks_reg, 1), t)

%% Plot SV magnitude over blocks vs. the movement metrics

vw = VideoWriter([data_dirpath, 'SV_magnitude.mp4'], 'MPEG-4');
vw.Quality = 100;
% Set the frame rate equal to the original video's
vw.FrameRate = 10; 
tl = tiledlayout(1, 2);
tl.TileSpacing = 'loose';
tl.Padding = 'loose';

open(vw)
fh = figure;
% fh.Position = [0, 0, 1920, 1080];
fh.Position = [0, 0, 2560, 1440];

ylims = [min(SVsallBlocks, [], 'all'), max(SVsallBlocks, [], 'all')];
for fi = 1:size(SVsallBlocks, 2)
% for fi = 1
    nexttile(1)
    semilogy(SVsallBlocks(:, fi), 'LineWidth', 4)
    ylim(ylims)
    xlabel('SV number')
    ylabel('SV magnitude')

    nexttile(2)
    plot(TD.daqTimeTags, accel_zm ./ max(accel_zm), '--') % Plot normalized, zero-meaned accelerometer magnitude
    hold on
    plot(t, GVTD ./ max(GVTD), 'LineWidth', 2); xlabel("Time [s]"); ylabel("GVTD")
    plot(t, (PDI_reg_global_mean - mean(PDI_reg_global_mean)) ./ max(PDI_reg_global_mean))
    xline(t(fi), 'g-', 'LineWidth', 3)
    ylabel("Accelerometer amplitude")
    legend("Accelerometer magnitude", "GVTD", "PDI global mean", "Current time")
    hold off
    
    frame = getframe(fh);
    im = frame2im(frame);
    writeVideo(vw, im)
end
close(vw)
clearvars fi ylims frame im

%% Edit video so it lines up with the ultrasound acq start time
% Get timestamp of final block
dt = (bs * (1-bo))/P.frameRate; % Time step calculated from the block size, overlap, and frame rate
% t = 0:dt:(numBlocks - 1)*dt; % Time stamps of each block

t_lastblock = numBlocks * dt; % [s]

sf_to_video_delay = seconds(TD.sfStartTimetag - videoData.startTimetag); % Time delay between video start and ultrasound acquisition start

% Crop video so it lies within the ultrasound times
videoCropped = videoData.frames;
videoOrigTimestamps = videoData.timestamps;
us_cam_mask = videoData.timestamps > sf_to_video_delay & videoData.timestamps < t_lastblock;
videoCroppedTimes = videoData.timestamps(us_cam_mask) - sf_to_video_delay;
videoCropped = videoCropped(:, :, us_cam_mask);% test = resample(videoData.frames, t, dimension = 3);
clearvars videoData % Save memory

%% Plot GVTD and accelerometer with the behavioral camera
vw = VideoWriter([data_dirpath, 'cam_movement_metrics.mp4'], 'MPEG-4');
vw.Quality = 100;
% Set the frame rate equal to the original video's
vw.FrameRate = round(1/mean(diff(videoOrigTimestamps))); 

open(vw)
fh = figure;
% fh.Position = [0, 0, 1920, 1080];
fh.Position = [0, 0, 2560, 1440];
% subplot(1, 2, 1)
tl = tiledlayout(1, 2);
tl.TileSpacing = 'loose';
tl.Padding = 'loose';

for fi = 1:size(videoCropped, 3)
% for fi = 1:4
    % subplot(1, 2, 1)
    nexttile(1)
    imagesc(videoCropped(:, :, fi)); colormap gray
    axis equal
    axis tight

    % subplot(1, 2, 2)
    nexttile(2)
    % if fi == 1
        plot(TD.daqTimeTags, accel_zm ./ max(accel_zm), '--') % Plot normalized, zero-meaned accelerometer magnitude
        hold on
        plot(t, GVTD ./ max(GVTD), 'LineWidth', 2); xlabel("Time [s]"); ylabel("GVTD")
        plot(t, (PDI_reg_global_mean - mean(PDI_reg_global_mean)) ./ max(PDI_reg_global_mean))
        xline(videoCroppedTimes(fi), 'g-', 'LineWidth', 3)
        ylabel("Accelerometer amplitude")
        legend("Accelerometer magnitude", "GVTD", "PDI global mean", "Current time")
        hold off
    % else
    %     xline(videoCroppedTimes(fi), 'g-', 'LineWidth', 3)
    % end
    
    frame = getframe(fh);
    im = frame2im(frame);
    writeVideo(vw, im)
end
close(vw)
