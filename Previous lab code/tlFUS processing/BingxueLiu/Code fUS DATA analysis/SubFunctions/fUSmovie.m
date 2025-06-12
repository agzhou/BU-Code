function fUSmovie(Xmovi,ImgInfo)
% generate a movie for continous fUS acquisitions within t
% Moviie.avi auto save and override to the data directory
% 
% ImgInfo.savePathName=[FilePath,FileName];
% ImgInfo.t = 500/1000;
%
% Last modified: 1/5/2021, Bingxue Liu

[nx,ny,nf] = size(Xmovi);

cmin = min(min(min(Xmovi)));
cmax = max(max(max(Xmovi)));

figure;
Moviie = VideoWriter([ImgInfo.savePathName,'-AVI.avi']);
Moviie.Quality = 100;
Moviie.FrameRate = 10;
open(Moviie);

for i = 1 : nf
    % make a movie
    imagesc(Xmovi(:,:,i)); axis image; colormap(hot); colorbar; hold on;
    caxis([cmin, cmax]); 
    hColorbar = colorbar;
    set(hColorbar, 'Ticks', sort([hColorbar.Limits, hColorbar.Ticks]));
        text(20,20,['Time = ' num2str((i)*ImgInfo.t/nf),'s'],...
            'Position',[150,20],...
            'Units','pixels',...
            'FontSize',12,'Color',[0,1,0])
        hold off;
%    pause (1)
    vframe = getframe(gcf); 
    writeVideo(Moviie, vframe);
   figure(gcf)
end
close(Moviie);