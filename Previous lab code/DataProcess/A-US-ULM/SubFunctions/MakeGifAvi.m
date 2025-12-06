%% Make .gif and .avi file
function MakeGifAvi(IMGdata, ImgInfo,BKdata)
% ImgInfo: xCoor, yCoor,cMap,cAxis,tIntFrame,delayTime,savePathName,
% textX,textZ;
%% example
% [nz,nx]=size(dIQ(:,:,1));
% for it=1:100
%     IMGdata(:,:,it)=imresize(conv2(BB(:,:,it),PRSSinfo.sysPSF),[nz,nx]);
%     BKdata(:,:,it)=abs(dIQ(:,:,it));
% end
% FilePath='G:\PROJ-R-Stroke\1113Stroke7\Baseline\vUS\';
% FileName='SVD-Speckle-Phase-ZoomIn';
% ImgInfo.xCoor=[1:size(Data,2)]*PixSize;
% ImgInfo.yCoor=[1:size(Data,1)]*PixSize;
% ImgInfo.cMap='parula';
% ImgInfo.textX=ImgInfo.xCoor(floor(end/3));
% ImgInfo.textZ=ImgInfo.yCoor(floor(end/16));
% ImgInfo.cAxis=[-4 4]*1e0;
% %ImgInfo.cAxis=[mean(Data(:))-std(Data(:)) max(Data(:))*2/4]*1e0;
% ImgInfo.tIntFrame=0.0002;
% ImgInfo.delayTime=0.1;
% ImgInfo.savePathName=[FilePath,FileName];
% MakeGifAvi((IMGdata(:,:,:)), ImgInfo, BKdata)
% or:
% MakeGifAvi((IMGdata(:,:,:)), ImgInfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nx,ny,nf]=size(IMGdata);
[VzCmap, VzCmapUp, VzCmapDn]=Colormaps_fUS;
for iSeg=1:ceil(nf/10)
    hF=figure;
    set(hF,'Position',[300 400 700 400])
    for ifile=(iSeg-1)*10+1:min(iSeg*10,nf)
        if nargin<3
            hAxes1 = axes;
            colormap(ImgInfo.cMap);
            img=squeeze(IMGdata(:,:,ifile));
            h1=imagesc(ImgInfo.xCoor,ImgInfo.yCoor,img);
            caxis(hAxes1,ImgInfo.cAxis);
            axis equal
            axis tight
            axis off
            colorbar
            textshow=text(ImgInfo.textX,ImgInfo.textZ,['iFrame=',num2str(ifile),', t=',num2str(ImgInfo.tIntFrame*ifile),' s']);
%             textshow=text(ImgInfo.textX,ImgInfo.textZ,['t=0-',num2str(ImgInfo.tIntFrame*ifile),' s']);
            textshow.Color='red';
            textshow.FontSize=16;
            textshow.FontWeight='bold';
            % draw and save the gif file
            drawnow;
            axis off
            pause(0.002);
            set(gcf,'color','white')
        else
            if size(BKdata,3)>1
                ImgBK=BKdata(:,:,ifile);
            else
                ImgBK=BKdata;
            end
            ImgOLP=squeeze(IMGdata(:,:,ifile));
            % background
            hAxes1 = axes;
            imagesc(ImgInfo.xCoor,ImgInfo.yCoor,ImgBK);
            caxis(hAxes1,[0 1.5])
            colormap(hAxes1,gray)            
%             caxis(hAxes1,[-30 30])
%             colormap(hAxes1,VzCmap)
            colorbar('Ticks',[-10 10], 'TickLabels',{[],[]} );
            axis equal tight
            axis(hAxes1,'off')
            
            hold on;
            hAxes2 = axes;
            h2=imagesc(ImgInfo.xCoor,ImgInfo.yCoor,ImgOLP);
%             set(h2,'AlphaData',(ImgOLP>0.5))
            set(h2,'AlphaData',(ImgOLP>0.25))
            colormap(hAxes2,ImgInfo.cMap);
            caxis(hAxes2,ImgInfo.cAxis);
            colorbar
            hold off
            axis equal tight
            axis(hAxes2,'off')
            
            linkaxes([hAxes1,hAxes2])
%             textshow=text(ImgInfo.textX,ImgInfo.textZ,['iFrame=',num2str(ifile),', t=',num2str(ImgInfo.tIntFrame*ifile),' s']);
            textshow=text(ImgInfo.textX,ImgInfo.textZ,['t=','0-', num2str(ImgInfo.tIntFrame*ifile),' s']);
            textshow.Color='red';
            textshow.FontSize=16;
            textshow.FontWeight='bold';
            axis off
            drawnow; pause(0.05);
        end
        frames(ifile)=getframe(gcf);
    end
    close (hF)
    pause(0.2)
end
outfileGIF = [ImgInfo.savePathName,'-GIF.gif'];
newVid = VideoWriter([ImgInfo.savePathName,'-AVI.avi']);
% newVid = VideoWriter([ImgInfo.savePathName,'Uncompressed AVI']);
newVid.FrameRate = 1/ImgInfo.delayTime;
newVid.Quality = 100;
open(newVid)
for ifile=1:nf
    % save GIF 
    im = frame2im(frames(ifile));
    [imind,cm] = rgb2ind(im,256);  
    if ifile==1
        imwrite(imind,cm,outfileGIF,'gif','DelayTime',ImgInfo.delayTime,'loopcount',inf);
    else
        imwrite(imind,cm,outfileGIF,'gif','DelayTime',ImgInfo.delayTime,'writemode','append');
    end
    % save AVI video
    writeVideo(newVid,frames(ifile).cdata)%within the for loop saving one frame at a tim
end
close(newVid)