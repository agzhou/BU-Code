function MakeAvimPAT(Xmn, mWL, filepath)

figure;
Moviie = VideoWriter([filepath,'Moviie.avi']);
Moviie.Quality = 100;
Moviie.FrameRate = 2;
open(Moviie);
cmin = min(Xmn(:));
cmax = max(Xmn(:));
 
for i = 1:size(mWL,2)
    imagesc(Xmn(:,:,i)); axis image; colorbar; hold on;
    caxis([cmin cmax]); 
        text(20,20,['Wavelength ' num2str(mWL(i)),'nm'],...
            'Position',[150,20],...
            'Units','pixels',...
            'FontSize',12,'Color',[1,0,0])
        hold off;
   pause (1)
    vframe = getframe(gcf); 
    writeVideo(Moviie, vframe);
   figure(gcf)
end
close(Moviie);
