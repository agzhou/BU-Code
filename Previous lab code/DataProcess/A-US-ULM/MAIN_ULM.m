%% ultrasound localized microscopy data processing
clear all;
defaultpath='Z:\US-DATA\PROJ-D-vUS\invivo-vUSvsULM-20181113-Acute-ISO\DATA-CP6-vUSBB\ULM-Baseline-CP6';
addpath('.\SubFunctions')
[FileName,FilePath]=uigetfile(defaultpath);  % read data of a small part of the brain cortex (IQR matrix) 
fileInfo=strsplit(FileName,'-');
myFile=matfile([FilePath,FileName]);
P=myFile.P;
PRSSinfo=P;
%% data processing parameters
prompt={'File Start ','number of files', 'SVD rank low', 'SVD rank high',...
    'ULM refine pixel size [um]','PSF xCoef Threshold', 'P-P distance criteria [um]',...
    'Number of Tracking frame','minTrackable Frame','Final Image Pixel size [um]'};
name='IQ2ULM data processing';
defaultvalue={num2str(fileInfo{8}(1:end-4)),'50', '5', '200',...
    '10','0.6','100','10','4','5'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
startFile=str2num(numinput{1});
nFile=str2num(numinput{2});
PRSSinfo.SVDrank=[str2num(numinput{3}),min(str2num(numinput{4}),P.numCCframes)];
PRSSinfo.lPix=str2num(numinput{5});
PRSSinfo.PSFtd=str2num(numinput{6});
PRSSinfo.dCrit=str2num(numinput{7});
PRSSinfo.nTrack=str2num(numinput{8});
PRSSinfo.thdTrk=str2num(numinput{9});
PRSSinfo.IntPixSize=str2num(numinput{10});
PRSSinfo.NEQ=1;
PRSSinfo.HPfC=25;
PRSSinfo.rFrame=P.CCFR;
%% US system point spread function
PRSSinfo.rfn(1)=round(P.dzImg*1e3/PRSSinfo.lPix); % z refine scale
PRSSinfo.rfn(2)=round(P.dxImg*1e3/PRSSinfo.lPix); % x refine scale
PRSSinfo.Dim=[numel(P.zCoor),numel(P.xCoor)].*PRSSinfo.rfn;
FWHM_X=240; % lateral resolutin, FWHM-Amplitude, um
FWHM_Z=100;  % axial resolutin, FWHM-Amplitude, um
Sigma_X=FWHM_X/(2*sqrt(2*log(2)));
Sigma_Z=FWHM_Z/(2*sqrt(2*log(2)));
xPSF0=-30:30; zPSF0=xPSF0; % pixels
[xPSF,zPSF]=meshgrid(xPSF0,zPSF0);
PRSSinfo.sysPSF=exp(-((xPSF/(Sigma_X/PRSSinfo.lPix)).^2+(zPSF/(Sigma_Z/PRSSinfo.lPix)).^2)/2);

pathInfo=strsplit(FilePath,'\');
fileInfo=strsplit(FileName,'-');
SavePath=[strjoin(pathInfo(1:end-2),'\'),'\RESULTdmIQ-',pathInfo{end-1},'\'];
if ~exist(SavePath)
    mkdir(SavePath);
end
%% ULM calculation
for iFile=startFile:startFile+nFile-1
    clear BB BB0 BBV BBVz BBPD;
    iFileInfo=fileInfo;
    iFileInfo{8}=[num2str(iFile) iFileInfo{8}(end-3:end)];
    iFileName=strjoin(iFileInfo,'-');
    %% load data
    disp(['Loading data: ',num2str(iFileName)]);
    load ([FilePath, iFileName]);
    disp('Data loaded!');
    [nz0,nx0,nt]=size(IQ);
    nz=nz0*PRSSinfo.rfn(1);
    nx=nx0*PRSSinfo.rfn(2);
    %% T1. brain ROI selection
%     IQm=mean(abs(IQ(:,:,1:100)),3);
%     figure,imagesc(IQm);
%     caxis([0 5e7])
%     [loc_x, loc_y]=ginput(15);  % no image rotation,
%     Bd=fnplt(cscvn([loc_y';loc_x']));
%     BW=roipoly(IQm,Bd(2,:),Bd(1,:));
%     hold on, imagesc(BW.*IQm);
    %% T2. Motion correction, not really help
%     bIQ=SVDfilter(IQ,[1 2]); % tissue
%     refIQ(:,:,iFile)=log(abs(bIQ(:,:,floor(nt/2)+1)));
%     for it=1:nt
%         [D(:,:,:,it),~] = imregdemons(log(abs(bIQ(:,:,it))),refIQ(:,:,iFile), 50);
%         regIQ(:,:,it)=imwarp((IQ(:,:,it)),D(:,:,:,it));
%     end
%     bIQ=SVDfilter(IQ,[1 2]); % tissue
%     refIQ(:,:,iFile)=log(abs(bIQ(:,:,floor(nt/2)+1)));
%     for it=1:nt
%         [D(it)] = imregcorr(log(abs(bIQ(:,:,it))), refIQ(:,:,iFile),'rigid');
%         regIQ(:,:,it)=imwarp((IQ(:,:,it)),D(it));
%     end
    %% 1. IQ to BB
    % frame-to-frame substraction
    %% 1. dIQ
%     dIQ=diff(IQ,1,3);
%     dIQ(:,:,nt)=IQ(:,:,end)-IQ(:,:,end-2);
%     [dCoorBB]=dIQ2BB(abs(dIQ), PRSSinfo);
    %% 2.  dmIQ
    mIQ=movmean(IQ,35,3);
    dIQ=IQ-mIQ;
    [dCoorBB]=dIQ2BB(abs(dIQ), PRSSinfo);
    %% 3. sIQ
    [sIQ, sIQHP, sIQHHP, eqNoise]=IQ2sIQ(IQ,PRSSinfo); % 0: no noise equalization
%     [dCoorBB]=dIQ2BB(abs(sIQ), PRSSinfo);
    %% 4. dmsIQ
%     [sIQ, sIQHP, sIQHHP, eqNoise]=IQ2sIQ(IQ,PRSSinfo); % 0: no noise equalization
%     msIQ=movmean(abs(sIQ),201,3);
%     dIQ=abs(sIQ)-msIQ;
%     [dCoorBB]=dIQ2BB(abs(dIQ), PRSSinfo);
    %% 
    CoorBB=dCoorBB;
    % 1.1 Coordinate to 3D image
    for it=1:nt
        BBcnt=zeros(nz, nx);
        itCoorBB=CoorBB{it};
        Ind=sub2ind(size(BBcnt),max(round(itCoorBB(:,1)),1),max(round(itCoorBB(:,2)),1));
        BBcnt(Ind)=1;
        BB0(:,:,it)=BBcnt;
    end
    figure,imagesc(sum(BB0(:,:,:),3));colormap(hot);caxis([0 3])
     % 1.2 show individual frame  
%     figure
%     for it=1:1:300
%         CB=ones(5,5);
%         BBshow=convn(BB0(:,:,it),CB,'same');
%         imagesc(BBshow);
%         colormap(hot)
%         axis equal tight
%         caxis([0 2]);
%         drawnow;
%         pause(0.5);
%     end
    %% 2. BB pair and track
    disp(['Pair and track BB...']);
    [BBPD]=BBPT(CoorBB,PRSSinfo);
    %% 3. vULM calculation
    disp(['vULM calculation']);
    [BB,BBV, BBVz]=BBPD2vBB1(BBPD,PRSSinfo);
    %% 4. Plot and save results
    PRSSinfo.Dim=[nz,nx,nt];
    Coor.z=[1:PRSSinfo.Dim(1)]/PRSSinfo.lPix;
    Coor.x=[1:PRSSinfo.Dim(2)]/PRSSinfo.lPix;
    PRSSinfo.Coor=Coor;
    SaveName=['ULM-',strjoin(fileInfo(2:end-1),'-'),'-',num2str(iFile)];
    
    BBup=-1*(BB(:,:,1)).^(1/2);
    BBdn=(BB(:,:,2)).^(1/2);
    
%     figure;
%     Fuse2Images(BBV(:,:,1),BBV(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
%     title('V, vULM [mm/s]')
%     
%     FigULM=figure;
%     Fuse2Images(BBup,BBdn,[-5 5],[-5 5],Coor.x,Coor.z,1);
%     title('ULM')

    FigULM=figure('visible','off');
%     FigULM=figure;
    Fuse2Images(BBup,BBdn,[-5 5],[-5 5],Coor.x,Coor.z,1);
    title('ULM')
    FigvULM=figure('visible','off');
%     FigvULM=figure;
%     set(FigvULM,'Position',[400 400 1600 450])
%     subplot(1,2,1)
    Fuse2Images(BBV(:,:,1),BBV(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
    title('V, vULM [mm/s]')
%     
%     subplot(1,2,2)
%     Fuse2Images(BBVz(:,:,1),BBVz(:,:,2),[-25 25],[-25 25],Coor.x,Coor.z,2.5);
%     title('V_z, vULM [mm/s]')
%     
    saveas(FigULM,[SavePath, SaveName,'.tif'],'tif');
    saveas(FigULM,[SavePath, SaveName,'.fig'],'fig');
    saveas(FigvULM,[SavePath, 'v',SaveName,'.tif'],'tif');
    saveas(FigvULM,[SavePath, 'v',SaveName,'.fig'],'fig');
    save([SavePath,SaveName,'.mat'],'-V7.3','BB','BBV','BBVz','BBPD','CoorBB','PRSSinfo','Coor')
%     disp([SaveName, ' processed and saved']);
end 
%% Check BB raw data, BB identification, and after paired and tracked BB
%% 1. Plot and Save BB raw data, dIQ/sIQ gif&avi data
clear BKdata IMGdata
[nz0,nx0,nt]=size(IQ);
% % 1. sIQ obtained with SVD or with HP
% PRSSinfo.SVDrank=[15,500];
% PRSSinfo.NEQ=1;
% PRSSinfo.HPfC=25;
% PRSSinfo.rFrame=P.CCFR;
% [sIQ, sIQHP, sIQHHP, eqNoise]=IQ2sIQ(IQ,PRSSinfo); % 0: no noise equalization
% % 2. dIQ obtained with IQ(t)=IQ(t+1)-IQ(t);
% % dIQ=diff(IQ,1,3);
% % dIQ(:,:,nt)=IQ(:,:,end)-IQ(:,:,end-2);
% % 3. dIQ obtained with IQ(t)-mean(IQ(t-10:t+10))
% % IQMov=movmean(IQ,10,3);
% % dIQ=IQ-IQMov;

cIQ=dIQ;
vNorm=max(abs(cIQ(:)));
nRfn=1;
for it=1:200
    IMGdata(:,:,it)=imresize(abs(cIQ(:,:,it))/vNorm*20,[nz0*nRfn,nx0*nRfn]);
%     IMGdata(:,:,it)=imresize((cIQ(:,:,it)),[nz0*nRfn,nx0*nRfn]);
end

PixSize=0.1; % mm
gifSavePath=FilePath;
SaveName='dmIQ';
ImgInfo.xCoor=[1:size(IMGdata,2)]*PixSize;
ImgInfo.yCoor=[1:size(IMGdata,1)]*PixSize;
ImgInfo.cMap='gray';
% ImgInfo.cAxis=[0.01 0.1];
ImgInfo.textX=ImgInfo.xCoor(floor(end*2/10));
ImgInfo.textZ=ImgInfo.yCoor(floor(end*15/16));
ImgInfo.cAxis=[0.05 2.5];
ImgInfo.tIntFrame=0.002;
ImgInfo.delayTime=0.15;
ImgInfo.savePathName=[gifSavePath,SaveName];
MakeGifAvi(IMGdata(:,:,:), ImgInfo)
% % plot a 2D image at a single time point
% % Fig=figure; set(Fig, 'Position',[300 400 900 600]); imagesc(IMGdata(:,:,20)); axis equal tight; colormap(gray); caxis([0 1.5]); colorbar
%% 2, Plot and Save dIQ vs sIQ
% [VzCmap, VzCmapUp, VzCmapDn]=Colormaps_fUS;
% clear BKdata IMGdata
% vNorm=max(abs(dIQ(:)));
% % dIQ vs sIQ
% for it=1:300
%     BKdata(:,:,it)=abs(dfIQ(:,:,it))/vNorm*20;
%     IMGdata(:,:,it)=abs(dIQ(:,:,it))/vNorm*15;
% end
% 
% PixSize=0.05; % mm
% gifSavePath=FilePath;
% SaveName='dIQ vs dfIQ';
% ImgInfo.xCoor=[1:size(IMGdata,2)]*PixSize;
% ImgInfo.yCoor=[1:size(IMGdata,1)]*PixSize;
% ImgInfo.cMap='hot';
% % ImgInfo.cAxis=[0.01 0.1];
% ImgInfo.textX=ImgInfo.xCoor(floor(end*6/8));
% ImgInfo.textZ=ImgInfo.yCoor(floor(end*15/16));
% ImgInfo.cAxis=[0 1];
% ImgInfo.tIntFrame=0.002;
% ImgInfo.delayTime=0.2;
% ImgInfo.savePathName=[gifSavePath,SaveName];
% MakeGifAvi(IMGdata(:,:,:), ImgInfo,BKdata(:,:,:))
%% 3, Plot and Save dIQ/sIQ-based background and BB identified gif&avi
[VzCmap, VzCmapUp, VzCmapDn]=Colormaps_fUS;
clear BKdata IMGdata
vNorm=max(abs(sIQ(:)));
for it=1:300
    BKdata(:,:,it)=imresize(abs(sIQ(:,:,it))/vNorm*15,[nz,nx]);
    IMGdata(:,:,it)=conv2(BB0(:,:,it),PRSSinfo.sysPSF,'same');
end

PixSize=0.05; % mm
gifSavePath=FilePath;
SaveName='sIQBB vs sIQ';
ImgInfo.xCoor=[1:size(IMGdata,2)]*PixSize;
ImgInfo.yCoor=[1:size(IMGdata,1)]*PixSize;
ImgInfo.cMap='hot';
% ImgInfo.cAxis=[0.01 0.1];
ImgInfo.textX=ImgInfo.xCoor(floor(end*1/8));
ImgInfo.textZ=ImgInfo.yCoor(floor(end*15/16));
ImgInfo.cAxis=[0 1];
ImgInfo.tIntFrame=3;
ImgInfo.delayTime=0.2;
ImgInfo.savePathName=[gifSavePath,SaveName];
MakeGifAvi(IMGdata(:,:,:), ImgInfo,BKdata(:,:,:))
% % 
%% 3, Plot and Save dIQ/sIQ-based background and BB paired and tracked gif&avi
% clear BKdata IMGdata
% vNorm=max(abs(sIQ(:)));
% pdBB=zeros(size(CoorBB(:,:,1:400)));
% for it=1:400
% %     pdBB=zeros(nz,nx);
% %     iBBPD=BBPD{it};
% %     iBBPD(:,:,iBBPD(1,1,:)==0)=[];
% %     iBBPD=round(iBBPD);
% %     pdBB(sub2ind(size(pdBB),iBBPD(1,1,:),iBBPD(2,1,:)))=1;
% %     IMGdata(:,:,it)=conv2(pdBB,PRSSinfo.sysPSF,'same');
% %     BKdata(:,:,it)=imresize(abs(sIQ(:,:,it))/vNorm*15,[nz,nx]);
%     IMGdata(:,:,it)=conv2(CoorBB(:,:,it),PRSSinfo.sysPSF,'same');
%     BKdata(:,:,it)=conv2(sBB(:,:,it),PRSSinfo.sysPSF,'same');
% end
% PixSize=0.25; % mm
% gifSavePath=FilePath;
% SaveName='sBBvsdBB';
% ImgInfo.xCoor=[1:size(IMGdata,2)]*PixSize;
% ImgInfo.yCoor=[1:size(IMGdata,1)]*PixSize;
% ImgInfo.cMap='hot';
% % ImgInfo.cAxis=[0.01 0.1];
% ImgInfo.textX=ImgInfo.xCoor(floor(end/3));
% ImgInfo.textZ=ImgInfo.yCoor(floor(end/18));
% ImgInfo.cAxis=[0.2 1.2];
% ImgInfo.tIntFrame=0.002;
% ImgInfo.delayTime=0.1;
% ImgInfo.savePathName=[gifSavePath,SaveName];
% MakeGifAvi(IMGdata(:,:,:), ImgInfo,BKdata(:,:,:))
% % 
% % %% 4. BB accumulation
for it=1:200
    IMGdata(:,:,it)=sum(BB0(:,:,1:it*5),3);
end
PixSize=0.5; % mm
gifSavePath=FilePath;
SaveName='t20-sIQ-BBacum';
ImgInfo.xCoor=[1:size(IMGdata,2)]*PixSize;
ImgInfo.yCoor=[1:size(IMGdata,1)]*PixSize;
ImgInfo.cMap='hot';
% ImgInfo.cAxis=[0.01 0.1];
ImgInfo.textX=ImgInfo.xCoor(floor(end/3));
ImgInfo.textZ=ImgInfo.yCoor(floor(end/16));
ImgInfo.cAxis=[0 3];
ImgInfo.tIntFrame=0.01;
ImgInfo.delayTime=0.1;
ImgInfo.savePathName=[gifSavePath,SaveName];
MakeGifAvi(IMGdata(:,:,:), ImgInfo)
% % % plot a 2D image at a single time point
% % Fig=figure; set(Fig, 'Position',[300 400 900 600]); imagesc(IMGdata(:,:,20)); axis equal tight; colormap(gray); caxis([0 1.5]); colorbar