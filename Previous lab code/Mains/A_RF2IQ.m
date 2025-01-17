%% Coherence compounding Beamform RF data
%% RF2IQ data process
clear all;
DeftPath='/projectnb/npbfus/NonStash/Jianbo/DATA/0405Resolution/XZ Resolution/X0/';
[filename,datapath]=uigetfile(DeftPath);
addpath('D:\CODE\Functions');

% % dataFile=matfile([datapath, filename]);
% % P=dataFile.P;
load([datapath, filename]);
RF0=reshape(RFRAW,[P.actZsamples,P.numAngles,P.numCCframes,P.nCh]);
%% check RF data 
% figure,imagesc(abs(squeeze(RF0(1:end-50,6,1,:))))
% Pres=max(abs(squeeze(RF0(100:600,6,1,:))));
% x=linspace(0 ,12.8,128);
% HW=fwhm(x,Pres);
% figure,plot(Pres)
% % title(['FWHM=',num2str(HW),' mm'])
% xlabel('Transducer element')
%%
P.nRFref=2;
if P.nRFref==1
    RF=RF0;
else
    for iCC=1:P.numCCframes
        for iAgl=1:P.numAngles
            iRF0=squeeze(RF0(:,iAgl,iCC,:));
            RF(:,iAgl,iCC,:)=imresize(iRF0,[P.actZsamples*P.nRFref,P.nCh]);
        end
    end
end
clear RF0
fileInfo=strsplit(filename,'-'); % File info
startDepthDAQ=P.startDepth;
% wavelength=P.vSound/(P.TWfrequency*1e3); % mm
wavelength=P.vSound/(18*1e3); % mm
zDelay=startDepthDAQ*wavelength/2;
%% calculate the reference matrix
xCoor=linspace(0,P.pitch*P.nCh, P.nCh*2);
zCoor=linspace(zDelay,(zDelay+15), 15/wavelength*2);
nx=length(xCoor);
nz=length(zCoor);
IndCtriMatrix=zeros(nz,nx,P.nCh,P.numAngles);
ApodChn=zeros(nz,nx,P.nCh,P.numAngles);
P.startDepthRec=2*P.startDepth-15; % calibrated startDepth value
[IndCtriMatrix,ApodChn]=RefIndMatrix(P,xCoor,zCoor, 0.28);
%% calculate one compounded image for reference
for iCC=1 %:P.numCCframes
    for iAgl=1:P.numAngles
        iRF=squeeze(RF(:,iAgl,iCC,:));
        iBF=zeros(nz,nx);
        for ix=1:nx
            for iz=1:nz
                [~, ctrChn]=find(ApodChn(iz,ix,:,iAgl)>0);
                for iCh=ctrChn(1):ctrChn(end)
                    iBF(iz,ix)=iBF(iz,ix)+iRF(IndCtriMatrix(iz,ix,iCh,iAgl),iCh);
                end
            end
        end
        IQ(:,:,iAgl,iCC)=iBF;%hilbert(iBF);
    end
    disp(['iCCFrame ', num2str(iCC),' is processed.'])
end
fig=figure;
set(fig,'Position',[200 500 400 300]);
imagesc(xCoor,zCoor,abs(hilbert(squeeze(sum(IQ(:,:,:,1),3)))));
colormap(gray)
xlabel('x [mm]');
ylabel('z [mm]');
axis equal tight
%% save beamformed data, only for save a few frames reconstructed using the script above
% P.xCoor=xCoor;
% P.zCoor=zCoor;
% fileInfo=strsplit(filename,'-');
% NameIQ0=[strjoin(fileInfo(1:7),'-'),'-','IQ0'];   % save name for coherence compounded IQ data
% disp('Data saving...');
% save([datapath,NameIQ0,'.mat'],'IQ','P');
% disp('Data saved!');

%% Beamform image parameter
prompt={['zImgStart, zDelay=', num2str(zDelay),' (mm)'],...
    ['zImg, zDelay=', num2str(zDelay),' (mm)'],...
    ['dzImg, wavelength=',num2str(wavelength),' (mm)'],...
    ['xImg, probeLength=',num2str(P.nCh*P.pitch),' (mm)'],...
    ['dxImg, pitch=',num2str(P.pitch),' (mm)'],'nRF refine','NA',...
    'Processing start file','number of files to process'};
name='Beamforming';
defaultvalue={num2str(zDelay),'10', num2str(wavelength/2), num2str(P.nCh*P.pitch),num2str(P.pitch/2),num2str(P.nRFref),'0.28',fileInfo{7},'1'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
zDelay=str2num(numinput{1});
zImg=str2num(numinput{2});
dzImg=str2num(numinput{3});
xImg=str2num(numinput{4});
dxImg=str2num(numinput{5});
P.nRFref=str2num(numinput{6});
NA=str2num(numinput{7});
startFile=str2num(numinput{8});
nFile=str2num(numinput{9});
clear IndCtriMatrix ApodChn xCoor zCoor
xCoor=[0:dxImg:xImg];
zCoor=[zDelay:dzImg:(zDelay+zImg)];
nx=length(xCoor);
nz=length(zCoor);
P.nSmplPerWvlnth=4;
nRFref=P.nRFref;
close (fig)
%% calculate the reference matrix

IndCtriMatrix=zeros(nz,nx,P.nCh,P.numAngles);
ApodChn=zeros(nz,nx,P.nCh,P.numAngles);
% P.startDepth=2*P.startDepth-15; % calibrated startDepth value
[IndCtriMatrix,ApodChn]=RefIndMatrix(P,xCoor,zCoor, NA);
xElemCoor=(0.5:P.nCh)*P.pitch;   
for iFile=startFile:startFile-1+nFile
    ifilename=[strjoin(fileInfo(1:6),'-'),'-',num2str(iFile),'-',strjoin(fileInfo(8:end),'-')];
    disp('Loading Data');
    load([datapath, ifilename]);
    RF0=reshape(RFRAW,[P.actZsamples,P.numAngles,P.numCCframes,P.nCh]);
    if nRFref==1
        RF=RF0;
    else
        for iCC=1:P.numCCframes
            for iAgl=1:P.numAngles
                iRF0=squeeze(RF0(:,iAgl,iCC,:));
                RF(:,iAgl,iCC,:)=imresize(iRF0,[P.actZsamples*nRFref,P.nCh]);
            end
        end
    end
    clear RF0
    %% remove high reflection signal
%     for iAgl=1:P.numAngles
%         for iCh=1:P.nCh
%             RF(abs(squeeze(mean(RF(:,iAgl,:,iCh),3)))>mean(abs(squeeze(mean(RF(:,iAgl,:,iCh),3))))*2+std(abs(squeeze(mean(RF(:,iAgl,:,iCh),3))))*2,iAgl,:,iCh)=0;
%         end
%     end
%     disp('Data loaded');
    
    tic
    % for iSup=1:P.numSupFrames
    IQ=single(zeros(nz,nx,P.numAngles,P.numCCframes));
    for iAgl=1:P.numAngles
        iRF=squeeze(RF(:,iAgl,:,:));
        for ix=1:nx
            for iz=1:nz
                [~, ctrChn]=find(ApodChn(iz,ix,:,iAgl)>0);
                for iCh=ctrChn(1):ctrChn(end)
                    WT=exp(-(xCoor(ix)-xElemCoor(iCh))^2/(2*(zCoor(iz)*NA)^2));
                    IQ(iz,ix,iAgl,:)=squeeze(IQ(iz,ix,iAgl,:))'+single(squeeze(iRF(IndCtriMatrix(iz,ix,iCh,iAgl),:,iCh)))*WT;
                end
            end
        end
        disp(['iAngle ', num2str(iAgl),' is processed.'])
    end
    % end
    toc
    %% Coherence compounded IQ data
    P.xCoor=xCoor;
    P.dxImg=dxImg;
    P.zCoor=zCoor;
    P.dzImg=dzImg;
    P.zDelayImg=zDelay;
    P.NA=NA;
    %% data saving
    %% IQ and beamformed data name
    NameIQ=[strjoin(fileInfo(1:6),'-'),'-',num2str(iFile),'-','IQ'];   % save name for coherence compounded IQ data
    P.nRFref=nRFref;
    disp('Data saving...');
    save([datapath,NameIQ,'.mat'],'IQ','P');
    disp('Data saved!');
    NameIQ
end
figure,imagesc(xCoor,zCoor,abs(hilbert(squeeze(sum(IQ(:,:,:,1),3)))));
colormap(gray)
xlabel('x [mm]');
ylabel('z [mm]');
axis equal tight