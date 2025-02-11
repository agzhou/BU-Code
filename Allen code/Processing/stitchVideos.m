% Stitch two .avi videos side by side and save
% (adapted from https://stackoverflow.com/a/30786910)

vid1 = VideoReader('F:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processed Data\video_10_80.avi');
vid2 = VideoReader('F:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processed Data with NLM\video_10_80.avi');

% videoPlayer = vision.VideoPlayer;

% new video
outputVideo = VideoWriter('F:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processed Data with NLM\video_10_80_stitched.avi');
outputVideo.FrameRate = vid1.FrameRate;
% outputVideo.Quality = vid1.Quality;
open(outputVideo);

while hasFrame(vid1) && hasFrame(vid2)
    img1 = readFrame(vid1);
    img2 = readFrame(vid2);

    if any(size(img1) ~= size(img2))

        img2rs = zeros(size(img1));
        for c = 1:size(img1, 3) % have to resize rgb dimension individually
            img2rs(:, :, c) = imresize(img2(:, :, c), [size(img1, 1), size(img1, 2)]);
        end
        imgt = horzcat(img1, img2rs);
    else
        imgt = horzcat(img1, img2);
    end

    % play video
%     step(videoPlayer, imgt);

    % record new video
    writeVideo(outputVideo, imgt);
end

% release(videoPlayer);
close(outputVideo);