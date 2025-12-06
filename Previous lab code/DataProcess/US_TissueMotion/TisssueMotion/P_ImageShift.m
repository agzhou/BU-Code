%% load BB result
FilePath='Z:\US-DATA\PROJ-R-STROKE\0808Stroke1\Ischemia\Ischemic\MB-CP7\'; 
FileName='BB-6-3-500-1000-1-CP-20';
load([FilePath,FileName]);
%% plot BB result
figure,imagesc(sum(BB,3));colormap(hot);caxis([0 3]);
%% load IQ data
load('Z:\US-DATA\PROJ-R-STROKE\0808Stroke1\Ischemia\Ischemic\MB-CP7\IQ-6-3-500-1000-1-CP-20.mat')
%% SVD filtering to get bulk tissue signal
cIQ=IQ(10:40,50:125,:);
[bulkIQ, Noise]=SVDfilter(cIQ,[1 10]); % sIQ: signal IQ
dIQ=bulkIQ-mean(bulkIQ,3);
[nz,nx,nt]=size(bulkIQ);
%% plot time seriese bulkIQ
% figure;
% for it=1:100
%     imagesc(abs((dIQ(:,:,it))));colormap(hot);caxis([0.1 10]*1e6)
%     title(num2str(it))
%     drawnow;pause(0.2);
% end
ImgInfo.xCoor=1:nx;
ImgInfo.yCoor=1:nz;
ImgInfo.cMap='hot';
% ImgInfo.cAxis=[0.01 0.1];
ImgInfo.cAxis=[0.2 5]*1e6;
ImgInfo.tIntFrame=0.002;
ImgInfo.delayTime=0.1;
ImgInfo.savePathName=[FilePath,'dIQ',FileName(3:end-4)];
MakeGifAvi(abs(dIQ(:,:,1:100)), ImgInfo)
%% image shift calculation - global phase 
for it=1:500
    B=abs(imresize(bulkIQ(:,:,it),20));
    C=abs(imresize(bulkIQ(:,:,501),20));
    FB=fftshift(fft2(B));
    FC=fftshift(fft2(C));
    D=FB.*conj(FC)./abs(FB.*conj(FC));
    F=(ifftshift(ifft2(D)));
    [M,I]=max(abs(F(:)));
    [zGPS(it),xGPS(it)]=find(abs(F)==M);
end
figure,subplot(2,1,1); plot(xGPS); title('xShift')
subplot(2,1,2); plot(zGPS); title('zShift')
%% image shift calculation - cross correlation
nxROI=1;nzROI=1;
[nPzROI, nPxROI, nt]=size(bulkIQ);
for ixROI=1:nxROI
    tic
    for izROI=1:nzROI
        iIQ=dIQ((izROI-1)*nPzROI+1:izROI*nPzROI,(ixROI-1)*nPxROI+1:ixROI*nPxROI,:);
        [xShift(:,izROI,ixROI), zShift(:,izROI,ixROI), Xcoef(:,izROI,ixROI)]=ImgSftCmp(iIQ, P,2,5,30);
    end
    toc
end
figure,subplot(2,1,1); plot(xShift); title('xShift')
subplot(2,1,2); plot(zShift); title('zShift')
%% remove very large shift
smallShift=find(abs(xShift)<10);
figure,imagesc(sum(BB(:,:,smallShift),3));colormap(hot);caxis([0 3]);

%% selcte ROI to remove
BW=roipoly(sum(BB,3));
%%
[nz,nx,nt]=size(BB);
for it=1:nt
    iROI=squeeze(BB(:,:,it)).*(uint8(BW));
    shifted(it)=squeeze(sum(iROI(:)));
end
figure,plot(shifted)
%%
smallShift=find(abs(shifted)<1);
figure,imagesc(sum(BB(:,:,smallShift),3));colormap(hot);caxis([0 3]);

%%
Path='Z:\US-DATA\PROJ-R-STROKE\0808Stroke\Ischemia\Ischemic\MB-CP7\';
NameBase='BB-6-3-500-1000-1-CP-';
for ifile=1:50
    iName=[NameBase,num2str(ifile)];
    load([Path,iName]);
    [nz,nx,nt]=size(BB);
    for it=1:nt
        iROI=squeeze(BB(:,:,it)).*(uint8(BW));
        shifted(it)=squeeze(sum(iROI(:)));
    end
    ShiftFrame=find(abs(shifted)>0);
    BB(:,:,ShiftFrame)=0;
    BBRmv(:,:,ifile)=sum(BB,3);
end


    
    