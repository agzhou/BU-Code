%
% File name: makeImgDataPMovie.m
%            simple script for loading in a folder of imgData files and generating
%            a movie from the images
%
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage NXT Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% Copyright © 2013-2025 Verasonics, Inc.

clear
close all

FRAMERATE = 83.3; %1/(240e-6 * 50)

dirname = "/media/verasonics/WD1/imgDataP/*.mat";

dirlist = dir(dirname);

for a = 1:length(dirlist)
    load([dirlist(a).folder '/' dirlist(a).name])

    imagesc((squeeze(ImgDataP.^.2)))
    caxis([1 22])
    colormap(gray)
    title(num2str(a))
    F(a)= getframe(gcf);
    drawnow
end
disp('--- Writing Video File ---')
writerObj = VideoWriter('replayVideo.avi');
writerObj.FrameRate = FRAMERATE;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
