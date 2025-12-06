%% load and plot Power Doppler results
%% for multiple files PDI process
clear all;
defaultpath='D:\OneDrive\Work\PROJ - FUS\PROJ - US velocimetry\Data - In vivo Validation\1113CP6\';
addpath('.\SubFunctions') % Path on JTOPTICS
[FileName,FilePath]=uigetfile(defaultpath);  % read data of a small part of the brain cortex (IQR matrix)  
fileInfo=strsplit(FileName(1:end-4),'-');
startCP0=regexp(fileInfo{7},'\d*','Match');
cpName=regexp(fileInfo{7},'\D*','Match');
if isempty(startCP0)
    startCP1='0';
else
    startCP1=startCP0{1};
end
%% data processing parameters
prompt={'Start Repeat', 'Number of Repeats','Start file (CP)','Number of files (CPs)','nRfn'};
name='File info';
defaultvalue={fileInfo{8},'1',num2str(startCP1), '1','1'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
startRpt=str2num(numinput{1});
nRpt=str2num(numinput{2});          % number of repeat for each coronal plane
startCP=str2num(numinput{3});
nCP=str2num(numinput{4});          % number of coronal planes
nRfn=str2num(numinput{5});          % image refine scale
%% LOAD AND PLOT
[VzCmap, VzCmapDn,VzCmapUp]=Colormaps_fUS;
for iCP=startCP:startCP+nCP-1
    for iRpt=startRpt:startRpt+nRpt-1
        iFileInfo=fileInfo;
        if nCP>1
            iFileInfo{7}=[cpName{1},num2str(iCP)];
        end              
        iFileInfo{8}=num2str(iRpt);
        iFileName=[strjoin(iFileInfo,'-'),'.mat'];
        if exist([FilePath,iFileName],'file')==2
            load([FilePath,iFileName]);
            myFile=matfile([FilePath,iFileName]);
            R=myFile.R;
            Mf=myFile.Mf;
            Vx=myFile.Vx;
            Vz=myFile.Vz;
%             Vz0=myFile.Vz0;
            PDI=myFile.PDI;
%             PDIHHP=myFile.PDIHHP;
            %                 PDISVD=myFile.PDISVD;
%             eqNoise=myFile.eqNoise;
            disp([iFileName,' is loaded!'])
        else
            disp(['Skipped ', iFileName])
        end
        [nz,nx,nPDI]=size(Vz);
        V=sign(Vz(:,:,1:2)).*sqrt(Vx.^2+Vz(:,:,1:2).^2);
        %% Refine image and plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         PDI=(PDI./eqNoise.^1.8).^0.5;
        [nzPDI,nxPDI,~]=size(PDI);
        PDIn=zeros(size(PDI));
        PDIrfn=zeros(nz*nRfn,nx*nRfn,2);
        Vrfn=zeros(nz*nRfn,nx*nRfn,2);
        Vzrfn=zeros(nz*nRfn,nx*nRfn,2);
        for iPDI=1:2
            PDItemp=squeeze(PDI(:,:,iPDI));
            PDIn(:,:,iPDI)=(PDItemp-min(PDItemp(:)))/(max(PDItemp(:))-min(PDItemp(:)));  % Normalized PDI, [0 1]
            PDIrfn(:,:,iPDI)=imresize(PDIn(:,:,iPDI),nRfn*nz/nzPDI);  % resize PDI for result plot
            
            Vrfn(:,:,iPDI)=imresize(V(:,:,iPDI),nRfn);
            Vzrfn(:,:,iPDI)=imresize(Vz(:,:,iPDI),nRfn);
        end
        %% x,z coordinates after image refined
        dz=P.dzImg*numel(P.zCoor)/nz;             % IQ dz 
        dx=P.dxImg*numel(P.xCoor)/nx;             % IQ dx
        Coor.x=linspace(0,nx*dx,nx*nRfn);
        Coor.z=linspace(P.zCoor(1),P.zCoor(1)+nz*dz,nz*nRfn);
        %% Result plot
        FigV=figure;
        set(FigV,'Position',[300 400 1000 800])
        Fuse2Images(V(:,:,1),V(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
        title(['V-',iFileName(1:end-4)]);
        saveas(FigV,[FilePath,'V',iFileName(4:end-4),'.png'],'png');
        saveas(FigV,[FilePath,'V',iFileName(4:end-4),'.fig'],'fig');
        FigVz=figure;
        set(FigVz,'Position',[300 400 1000 800])
        Fuse2Images(Vz(:,:,1),Vz(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
        title(['Vz-',iFileName(1:end-4)]);
        saveas(FigVz,[FilePath,'Vz',iFileName(4:end-4),'.png'],'png');
        saveas(FigVz,[FilePath,'Vz',iFileName(4:end-4),'.fig'],'fig');
        %% plot PDI-based HSV velocity map
%         figure;
%         PLOTwtV(V,(PDI).^0.4,Coor,[-30 30])
%         title('vUS-V [mm/s]');
        %% save figure without show
%         Coor.x=xCoor; Coor.z=zCoor;
%         f1=figure('visible','off');
%         PLOTwtV(Vrfn,PDIrfn,Coor,[-30 30])
%         title(['V-',iFileName(1:end-4)]);
%         saveas(f1,[FilePath,'VPDI',iFileName(4:end-4),'.png'],'png');
%         saveas(f1,[FilePath,'VPDI',iFileName(4:end-4),'.fig'],'fig');
% 
%         f2=figure('visible','off');
%         Fuse2Images(Vz(:,:,1),Vz(:,:,2),[-30 30],[-30 30],xCoor,zCoor,2.5);
%         title(['Vz-',iFileName(1:end-4)]);
%         saveas(f2,[FilePath,'Vz',iFileName(4:end-4),'.png'],'png');
%         saveas(f2,[FilePath,'Vz',iFileName(4:end-4),'.fig'],'fig');
    end
end


