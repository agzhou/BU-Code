% Stitch two .avi videos side by side and save
% (from https://stackoverflow.com/a/30786910)

vid1 = VideoReader('video1.avi');
vid2 = VideoReader('video2.avi');

videoPlayer = vision.VideoPlayer;

% new video
outputVideo = VideoWriter('newvideo.avi');
outputVideo.FrameRate = vid1.FrameRate;
open(outputVideo);

while hasFrame(vid1) && hasFrame(vid2)
    img1 = readFrame(vid1);
    img2 = readFrame(vid2);

    imgt = horzcat(img1, img2);

    % play video
    step(videoPlayer, imgt);

    % record new video
    writeVideo(outputVideo, imgt);
end

release(videoPlayer);
close(outputVideo);