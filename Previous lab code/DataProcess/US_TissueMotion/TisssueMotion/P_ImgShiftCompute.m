%% Bulk shift calculation using cross correlation 
%% ultrasound localized microscopy data processing
clear all;
defaultpath='Z:\US-DATA\0803PMPmouse\CP6-MB';

addpath('D:\OneDrive\Work\PROJ - FUS\CODE\Functions') % Path on JTOPTICS
[FileName,FilePath]=uigetfile(defaultpath);  % read data of a small part of the brain cortex (IQR matrix) 
fileInfo=strsplit(FileName,'-');
%% Select ROI
disp('Loading data...');
load ([FilePath, FileName]);
disp('Data loaded!');
SignalRank=[35 300];
cIQ=(IQ(:,:,1:400));
% cIQ=single(ID(:,:,1:400));
[nz,nx,nt]=size(cIQ);
[sIQ, Noise]=SVDfilter(cIQ,SignalRank); % sIQ: signal IQ
%  figure,imagesc(Noise)
% eqnIQ=sIQ./repmat(sNoiseMedNorm,[1,1,nfRpt]);
nsIQ=sIQ./repmat(Noise,[1,1,nt]); % blood signal obtained by SVD 
nsIQavg=MoveAvg(abs(nsIQ),10,3);
nsIQ=abs(nsIQ)./nsIQavg+median(nsIQavg,3);
PDI=mean((nsIQ).^2,3); 
PDIdb=10*log10(PDI./max(PDI(:))); % SVD-based PD image in dB

figure
imagesc(PDIdb);
caxis(MyCaxis(PDIdb));
colormap(hot);
colorbar;
title ('SVD-based PDI')
xlabel('X [PIX]'); ylabel('Z [PIX]');

[Xslt Zslt]=ginput(2);
Xslt=round(Xslt);
Zslt=round(Zslt);
ROI(:,1)=Zslt;
ROI(:,2)=Xslt;
%% data processing parameters
prompt={'File Start ','number of files', 'nSubROI-X','nSubROI-Z','SVD rank low', 'SVD rank high','ROI_x1','ROI_x2','ROI_z1','ROI_z2'};
name='CSD data processing';
defaultvalue={num2str(fileInfo{8}(1:end-4)),'50','3','3', '1', '5',num2str(ROI(1,2)),num2str(ROI(2,2)),num2str(ROI(1,1)),num2str(ROI(2,1))};
numinput=inputdlg(prompt,name, 1, defaultvalue);
startFile=str2num(numinput{1});
nFile=str2num(numinput{2});
nxROI=str2num(numinput{3});
nzROI=str2num(numinput{4});
SignalRank(1)=str2num(numinput{5});
SignalRank(2)=str2num(numinput{6});
ROI(1,2)=str2num(numinput{7});
ROI(2,2)=str2num(numinput{8});
ROI(1,1)=str2num(numinput{9});
ROI(2,1)=str2num(numinput{10});
nPxROI=floor(abs(diff(Xslt))/nxROI); % number of Xpixel for each ROI
nPzROI=floor(abs(diff(Zslt))/nzROI); % number of Zpixel for each ROI
%% Inner block image shift 
for iFile=startFile:startFile+nFile-1
    clear IQ;
    iFileInfo=fileInfo;
    iFileInfo{8}=[num2str(iFile) iFileInfo{8}(end-3:end)];
    iFileName=strjoin(iFileInfo,'-');
    disp('Loading data...');
    load ([FilePath, iFileName]);
    disp('Data loaded!');
    [nz,nx,nt]=size(IQ);
    cIQ=(IQ(min(ROI(:,1)):max(ROI(:,1)),min(ROI(:,2)):max(ROI(:,2)),:));
    SignalRank=[SignalRank(1) SignalRank(2)];
    [bulkIQ, Noise]=SVDfilter(cIQ,SignalRank); % sIQ: signal IQ
    RefPlane(:,:,iFile)=squeeze(bulkIQ(:,:,round(nt/2))); % in block reference plane
    for ixROI=1:nxROI
        tic
        for izROI=1:nzROI
            iIQ=bulkIQ((izROI-1)*nPzROI+1:izROI*nPzROI,(ixROI-1)*nPxROI+1:ixROI*nPxROI,:);
            [xShift(:,izROI,ixROI), zShift(:,izROI,ixROI), Xcoef(:,izROI,ixROI)]=ImgSftCmp(iIQ, P,1,5,30);
        end
        toc
    end
    iSaveName=['IS',iFileName(3:end)];
    save([FilePath,iSaveName],'xShift','zShift','Xcoef','ROI','P');
    disp([iSaveName]);
end
%% Inter block image shift calculation
for ixROI=1:nxROI
        tic
        for izROI=1:nzROI
            iIQ=RefPlane((izROI-1)*nPzROI+1:izROI*nPzROI,(ixROI-1)*nPxROI+1:ixROI*nPxROI,:);
            [xBLKShift(:,izROI,ixROI), zBLKShift(:,izROI,ixROI), XBLKcoef(:,izROI,ixROI)]=ImgSftCmp(iIQ, P,1,8,30);
        end
        toc
end
save([FilePath,'BLKIS.mat'],'xBLKShift','zBLKShift','RefPlane','ROI','P');
            
    
    