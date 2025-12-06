%% Brain tissue shift imaging
clear all;
defaultpath='H:\PJ - USI\';
%addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions');
[FileName,FilePath]=uigetfile(defaultpath);  % read data of a small part of the brain cortex (IQR matrix)  
disp('Loading data...');
load ([FilePath, FileName]);
disp('Data loaded!');
% IQ0=IQData{1};
% dataRAW=squeeze(IQ0(:,:,1,1,:));
% cIQ=hilbert(squeeze(sum(IQ,3))); % Coherence Compounding IQ data, then Hilbert transform to get complext Coherence compounded IQ data (envelop+phase)
cIQ=IQ;
[nz0,nx0,nt0]=size(cIQ);
clear IQ IQData;
ntP=nt0;
plotTrue=0;
%% SVD process
IQR=cIQ(:,:,1:min(nt0,ntP)); % IQ data used for data processing
[nz,nx,nt]=size(IQR);
rank=[25:300];
S=reshape(IQR,[nz*nx,nt]);
S_COVt=(S'*S);
[V,D]=eig(S_COVt); % V is the right singular Vector of S/eigenvector; D is the eigenvalue/square of Singular value
for it=1:nt 
    Ddiag(it)=abs(sqrt(D(it,it)));
end
Ddiag=20*log10(Ddiag/max(Ddiag)); % singular value in db
[Ddesc, Idesc]=sort(Ddiag,'descend');
% figure,plot(Ddesc);
for it=1:nt
    Vdesc(:,it)=V(:,Idesc(it));
end
UDelta=S*Vdesc;
%%%% Noise equalization 
Vnoise=zeros(size(Vdesc));
Vnoise(:,end)=Vdesc(:,end);
sNoise=reshape(UDelta*Vnoise',[nz,nx,nt]);
sNoiseMed=medfilt2(abs(squeeze(mean(sNoise,3))),[50 50],'symmetric');
sNoiseMedNorm=sNoiseMed/min(sNoiseMed(:));
%% SVD- based power doppler
Vrank=zeros(size(Vdesc));
Vrank(:,rank)=Vdesc(:,rank);
sBlood0=reshape(UDelta*Vrank',[nz,nx,nt]);
sBlood0=sBlood0./repmat(sNoiseMedNorm,[1,1,nt]); % blood signal obtained by SVD 
% sBlood0=diff(IQR,1,3);
% sBlood=sBlood0./repmat(sNoiseMedNorm,[1,1,nt-1]); % blood signal obtained by SVD 
% PDISVD=mean((abs(sBlood)).^2,3); 
sBloodavg=MoveAvg(abs(sBlood0),10,3);
sBlood=abs(sBlood0)./sBloodavg+median(sBloodavg,3);
PDISVD=mean((sBlood).^2,3); 
figure,imagesc(PDISVD);colormap hot
%% cross coefficient 
refIQ=squeeze(IQR(:,:,1));
for it=1:ntP
    iIQ=squeeze(IQR(:,:,it));
    xCoefRaw(it)=sum(iIQ(:).*conj(refIQ(:)))/(sqrt(sum(abs(iIQ(:)).^2))*sqrt(sum(abs(refIQ(:)).^2)));
end
tCoor=linspace(1/P.CCFR,ntP/P.CCFR,ntP)*1e3;
figure;
subplot(2,1,1);
plot(tCoor,abs(xCoefRaw));
xlabel('t [ms]');
ylabel('xCoef magnitude')
subplot(2,1,2);
plot(tCoor,angle(xCoefRaw));
xlabel('t [ms]');
ylabel('xCoef Phase')
%% calculate all ROIs
ROIx=[40:100:200];
ROIz=[30:100:150];
nxROI=length(ROIx);
nzROI=length(ROIz);
nRfnImg=10; % IQ image refine
nRfnCoef=50; % cross correlation map refine
nOrgPix=0.5;
for ixROI=1:nxROI-1
    tic
    for izROI=1:nzROI-1
        IQROI=cIQ(ROIz(izROI):min(ROIz(izROI)+3,nz0),ROIx(ixROI):min(ROIx(ixROI)+3,nx0),:);
        % IQ=abs(IQ); % image shift compensation for image intensity/magnitude
        [nz,nx,nt]=size(IQROI);
        
        IQRfn=imresize(squeeze(IQROI(:,:,1)),nRfnImg);
        RefIQ=IQRfn(nRfnImg+1:end-nRfnImg,nRfnImg+1:end-nRfnImg);
        for it=1:ntP
            iIQ0=squeeze(IQROI(:,:,it));
            iIQrfn=imresize(iIQ0,nRfnImg);
            for iShiftDx=-nOrgPix*nRfnImg:nOrgPix*nRfnImg
                iStartX=max(iShiftDx,1);iEndX=min(nx*nRfnImg+iShiftDx,nx*nRfnImg);
                iNx=iEndX-iStartX+1;
                for iShiftDz=-nOrgPix*nRfnImg:nOrgPix*nRfnImg
                    iIQ=iIQrfn(nRfnImg+1+iShiftDz:nRfnImg+iShiftDz+(nz-2)*nRfnImg,nRfnImg+1+iShiftDx:nRfnImg+iShiftDx+(nx-2)*nRfnImg);
                    Xcoef(iShiftDz+nOrgPix*nRfnImg+1,iShiftDx+nOrgPix*nRfnImg+1,it)=sum(iIQ(:).*conj(RefIQ(:)))/(sqrt(sum(abs(iIQ(:)).^2))*sqrt(sum(abs(RefIQ(:)).^2)));
                end
            end
        end
        % XcoefOrg=squeeze(Xcoef(nOrgPix*nRfnImg+1,nOrgPix*nRfnImg+1,:));
        XcoefOrg=squeeze(max(max(Xcoef,[],1),[],2));
        % plot xCoef map
        iXcoef=abs(squeeze(Xcoef(:,:,1)));
        zCoorCoef=[-nOrgPix*nRfnImg:nOrgPix*nRfnImg]*P.dzImg*1e3/nRfnImg;
        xCoorCoef=[-nOrgPix*nRfnImg:nOrgPix*nRfnImg]*P.dxImg*1e3/nRfnImg;
%         figure,imagesc(zCoorCoef,xCoorCoef,iXcoef)
%         colormap(jet)
%         colorbar
%         xlabel('X shift [um]')
%         ylabel('Z shift [um]')
%         axis equal tight
        %% refine cross correlation map and find the image shift
        zCoefRfn=linspace(-nOrgPix*nRfnImg,nOrgPix*nRfnImg,(2*nOrgPix*nRfnImg+1)*nRfnCoef);
        xCoefRfn=linspace(-nOrgPix*nRfnImg,nOrgPix*nRfnImg,(2*nOrgPix*nRfnImg+1)*nRfnCoef);
        for it=1:ntP
            iXcoef=abs(squeeze(Xcoef(:,:,it)));
            iXcoefFilt=imgaussfilt(iXcoef,2);
            iXcoefRfn=imresize(iXcoefFilt,nRfnCoef);
            [zMax(it),xMax(it)]=find(iXcoefRfn==max(iXcoefRfn(:)));
        end
        zShiftPix=zMax-(nOrgPix*nRfnImg+1/2)*nRfnCoef;
        xShiftPix=xMax-(nOrgPix*nRfnImg+1/2)*nRfnCoef;
        zShiftPeak(izROI,ixROI,:)=zCoefRfn(zMax)*P.dzImg*1e3/nRfnImg;
        xShiftPeak(izROI,ixROI,:)=xCoefRfn(xMax)*P.dxImg*1e3/nRfnImg;
    end
    toc
    ixROI
end
% %% save data
% save('H:\PJ - USI\PROJ - fUS Velocimetry\0611InvivoTest\Cardiac&Respiratory\10-5-CP3-xShiftPeak1000.mat','xShiftPeak');
% save('H:\PJ - USI\PROJ - fUS Velocimetry\0611InvivoTest\Cardiac&Respiratory\10-5-CP3-zShiftPeak1000.mat','zShiftPeak');
% save('H:\PJ - USI\PROJ - fUS Velocimetry\0611InvivoTest\Cardiac&Respiratory\10-5-CP3-xCoef1000.mat','xCoefRaw');
%% plot animation
if plotTrue==1
    %% plot result and save the movie, zShift and magnitude of xCoef
    xCoor=P.xCoor;zCoor=P.zCoor;
    [VzNegCmap, VzCmap, TshiftCmap]=Colormaps_fUS;
    figure,
    for it=1:ntP
        axes('Position', [0.1 0.4 0.8 0.5]);
        imagesc(xCoor,zCoor,squeeze(zShiftPeak(:,:,it)));
        colorbar;colormap(TshiftCmap);caxis([-4 4]);
        %     title(num2str(it));
        textshow=text(3.5,0.8,['t=',num2str(tCoor(it)),'ms']);
        textshow.Color='white';
        axis equal tight;
        xlabel('x [mm]');
        ylabel('z [mm]')
        
        axes('Position', [0.1 0.1 0.8 0.2]);
        plot(tCoor,abs(xCoefRaw))
        hold on;
        plot(tCoor(it),abs(xCoefRaw(it)),'ro')
        xlabel('t [ms]')
        ylim([0.99 1.01])
        
        drawnow;
        pause(0.2);
        frames(it)=getframe(gcf);
    end
    
    outfileGIF = [FilePath,FileName(1:end-4),'-GIFzMag.gif'];
    newVid = VideoWriter([FilePath,FileName(1:end-4),'-AVIzMag.avi']);
    newVid.FrameRate = 4;
    newVid.Quality = ntP;
    open(newVid)
    for ifile=1:ntP
        % save GIF
        im = frame2im(frames(ifile));
        [imind,cm] = rgb2ind(im,256);
        if ifile==1
            imwrite(imind,cm,outfileGIF,'gif','DelayTime',0.1,'loopcount',inf);
        else
            imwrite(imind,cm,outfileGIF,'gif','DelayTime',0.1,'writemode','append');
        end
        % save AVI video
        writeVideo(newVid,frames(ifile).cdata)%within the for loop saving one frame at a tim
    end
    close(newVid)
    %% plot result and save the movie, zShift and phase of xCoef
    figure,
    for it=1:ntP
        axes('Position', [0.1 0.4 0.8 0.5]);
        imagesc(xCoor,zCoor,squeeze(zShiftPeak(:,:,it)));
        colorbar;colormap(TshiftCmap);caxis([-4 4]);
        %     title(num2str(it));
        textshow=text(3.5,0.8,['t=',num2str(tCoor(it)),'ms']);
        textshow.Color='white';
        axis equal tight;
        xlabel('x [mm]');
        ylabel('z [mm]')
        
        axes('Position', [0.1 0.1 0.8 0.2]);
        plot(tCoor,angle(xCoefRaw))
        hold on;
        plot(tCoor(it),angle(xCoefRaw(it)),'ro')
        xlabel('t [ms]')
        ylim([-0.12 0.04])
        
        drawnow;
        pause(0.2);
        frames(it)=getframe(gcf);
    end
    
    outfileGIF = [FilePath,FileName(1:end-4),'-GIFzPhase.gif'];
    newVid = VideoWriter([FilePath,FileName(1:end-4),'-AVIzPhase.avi']);
    newVid.FrameRate = 4;
    newVid.Quality = ntP;
    open(newVid)
    for ifile=1:ntP
        % save GIF
        im = frame2im(frames(ifile));
        [imind,cm] = rgb2ind(im,256);
        if ifile==1
            imwrite(imind,cm,outfileGIF,'gif','DelayTime',0.1,'loopcount',inf);
        else
            imwrite(imind,cm,outfileGIF,'gif','DelayTime',0.1,'writemode','append');
        end
        % save AVI video
        writeVideo(newVid,frames(ifile).cdata)%within the for loop saving one frame at a tim
    end
    close(newVid)
    
    %% plot result and save the movie, xShift and magnitude of xCoef
    figure,
    for it=1:ntP
        axes('Position', [0.1 0.4 0.8 0.5]);
        imagesc(xCoor,zCoor,squeeze(xShiftPeak(:,:,it)));
        colorbar;colormap(TshiftCmap);caxis([-4 4]);
        %     title(num2str(it));
        textshow=text(3.5,0.8,['t=',num2str(tCoor(it)),'ms']);
        textshow.Color='white';
        axis equal tight;
        xlabel('x [mm]');
        ylabel('z [mm]')
        
        axes('Position', [0.1 0.1 0.8 0.2]);
        plot(tCoor,abs(xCoefRaw))
        hold on;
        plot(tCoor(it),abs(xCoefRaw(it)),'ro')
        xlabel('t [ms]')
        ylim([0.99 1.01])
        
        drawnow;
        pause(0.2);
        frames(it)=getframe(gcf);
    end
    
    outfileGIF = [FilePath,FileName(1:end-4),'-GIFxMag.gif'];
    newVid = VideoWriter([FilePath,FileName(1:end-4),'-AVIxMag.avi']);
    newVid.FrameRate = 4;
    newVid.Quality = ntP;
    open(newVid)
    for ifile=1:ntP
        % save GIF
        im = frame2im(frames(ifile));
        [imind,cm] = rgb2ind(im,256);
        if ifile==1
            imwrite(imind,cm,outfileGIF,'gif','DelayTime',0.1,'loopcount',inf);
        else
            imwrite(imind,cm,outfileGIF,'gif','DelayTime',0.1,'writemode','append');
        end
        % save AVI video
        writeVideo(newVid,frames(ifile).cdata)%within the for loop saving one frame at a tim
    end
    close(newVid)
    %% plot result and save the movie, xShift and phase of xCoef
    figure,
    for it=1:ntP
        axes('Position', [0.1 0.4 0.8 0.5]);
        imagesc(xCoor,zCoor,squeeze(xShiftPeak(:,:,it)));
        colorbar;colormap(TshiftCmap);caxis([-4 4]);
        %     title(num2str(it));
        textshow=text(3.5,0.8,['t=',num2str(tCoor(it)),'ms']);
        textshow.Color='white';
        axis equal tight;
        xlabel('x [mm]');
        ylabel('z [mm]')
        
        axes('Position', [0.1 0.1 0.8 0.2]);
        plot(tCoor,angle(xCoefRaw))
        hold on;
        plot(tCoor(it),angle(xCoefRaw(it)),'ro')
        xlabel('t [ms]')
        ylim([-0.12 0.04])
        
        drawnow;
        pause(0.2);
        frames(it)=getframe(gcf);
    end
    
    outfileGIF = [FilePath,FileName(1:end-4),'-GIFxPhase.gif'];
    newVid = VideoWriter([FilePath,FileName(1:end-4),'-AVIxPhase.avi']);
    newVid.FrameRate = 4;
    newVid.Quality = ntP;
    open(newVid)
    for ifile=1:ntP
        % save GIF
        im = frame2im(frames(ifile));
        [imind,cm] = rgb2ind(im,256);
        if ifile==1
            imwrite(imind,cm,outfileGIF,'gif','DelayTime',0.1,'loopcount',inf);
        else
            imwrite(imind,cm,outfileGIF,'gif','DelayTime',0.1,'writemode','append');
        end
        % save AVI video
        writeVideo(newVid,frames(ifile).cdata)%within the for loop saving one frame at a tim
    end
    close(newVid)
end