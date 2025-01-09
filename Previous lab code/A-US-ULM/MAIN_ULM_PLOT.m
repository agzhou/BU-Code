%% load and plot paired microbubble and velocity map
clear all;
defaultpath='Z:\US-DATA\PROJ-D-iCDfUS\D-invivoPMP16-20190924-ISO-vUS\CP8\RESULTsIQ-ULM-CP8-V[5  40]-dt[4ms]';
addpath('.\SubFunctions')
[FileName,FilePath]=uigetfile(defaultpath);  % read data of a small part of the brain cortex (IQR matrix) 
myFile=matfile([FilePath,FileName]);
% PRSSinfo=myFile.PRSinfo;  % SCC processing
PRSSinfo=myFile.PRSSinfo; % from local processing
fileInfo=strsplit(FileName,'-');
prompt={'File Start ', 'number of files'};
name='ULM data processing';
defaultvalue={num2str(fileInfo{8}(1:end-4)),'50'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
startFile=str2num(numinput{1});
nFile=str2num(numinput{2});
%% for multiple CP processing
fileInfo=strsplit(FileName,'-');
%%
[VzCmap, VzCmapUp, VzCmapDn]=Colormaps_fUS;
BBAll=zeros(PRSSinfo.Dim(1),PRSSinfo.Dim(2),2,nFile); 
BBVAll=zeros(PRSSinfo.Dim(1),PRSSinfo.Dim(2),2,nFile);
BBVzAll=zeros(PRSSinfo.Dim(1),PRSSinfo.Dim(2),2,nFile);
for iFile=startFile:startFile+nFile-1
    iFileInfo=fileInfo;
    iFileInfo{8}=[num2str(iFile) iFileInfo{8}(end-3:end)];
    iFileName=strjoin(iFileInfo,'-');
    if exist([FilePath, iFileName],'file')
        load([FilePath, iFileName]);
        disp([iFileName, ' Loaded!']);
        BBAll(:,:,:,iFile-startFile+1)=BB;
        BBVAll(:,:,:,iFile-startFile+1)=BBV;
        BBVzAll(:,:,:,iFile-startFile+1)=BBVz;
%         nBB(iFile)=numel(find(BB0>0));
    else
        disp([iFileName, ' skipped']);
    end
end
% save([FilePath, 'BBAll.mat'],'BBAll','BBVAll','BBVzAll','PRSSinfo','Coor');
%% Plot PVelocity map
[nz, nx,nt]=size(BBVAll(:,:,1,:));
startT=5;
endT=nt;
% All BB
BBup=-1*(sum(BBAll(:,:,1,startT:endT),4)).^(1/3);
BBdn=(sum(BBAll(:,:,2,startT:endT),4)).^(1/3);
BBall=squeeze(sum(sum(BBAll(:,:,:,startT:endT),4),3));
clear BB
BB(:,:,1)=BBup;
BB(:,:,2)=BBdn;
BB(:,:,3)=BBall;
% All V & Vz
BBV=zeros(nz,nx,3);
BBVz=zeros(nz,nx,3);
% BBV(:,:,1)=min(BBVupAll(:,:,startT:endT),[],3);
% BBV(:,:,2)=max(BBVdnAll(:,:,startT:endT),[],3);
% BBV(:,:,1)=median(BBVupAll(:,:,startT:endT),3);
% BBV(:,:,2)=median(BBVdnAll(:,:,startT:endT),3);
for iz=1:nz
    for ix=1:nx
        %% V
        upVtemp=squeeze(BBVAll(iz,ix,1,startT:endT));
        dnVtemp=squeeze(BBVAll(iz,ix,2,startT:endT));
        upV=median(upVtemp(abs(upVtemp)>0.5));
        dnV=median(dnVtemp(abs(dnVtemp)>0.5));
%         upV=mean(upVtemp(abs(upVtemp)>3));
%         dnV=mean(dnVtemp(abs(dnVtemp)>3));
        if ~isnan(upV)
            BBV(iz,ix,1)=upV; % upFlow
        end
        if ~isnan(dnV)
            BBV(iz,ix,2)=dnV; % downFlow
        end
        %% Vz
        upVztemp=squeeze(BBVzAll(iz,ix,1,startT:endT));
        dnVztemp=squeeze(BBVzAll(iz,ix,2,startT:endT));
        upVz=median(upVztemp(abs(upVztemp)>0.5));
        dnVz=median(dnVztemp(abs(dnVztemp)>0.5));
%         upVz=mean(upVztemp(abs(upVztemp)>3));
%         dnVz=mean(dnVztemp(abs(dnVztemp)>3));
        if ~isnan(upVz)
            BBVz(iz,ix,1)=upVz; % upFlow
        end
        if ~isnan(dnVz)
            BBVz(iz,ix,2)=dnVz; % downFlow
        end
    end
end
BBV(:,:,3)=max(abs(BBV),[],3);
BBVz(:,:,3)=max(abs(BBVz),[],3);

% BBmsk=zeros(size(BB));
% BBmsk(abs(BB)>2)=1;
BBmsk=1;
BBV=BBV.*BBmsk;
BBVz=BBVz.*BBmsk;
%% plot V
vCrange=[-30 30];
fig1=figure;
h1=axes;
imagesc(BBV(:,:,1)); % up flow
% BBVup=(LWMEAN(BBV(:,:,1),5));
% imagesc(BBVup); % up flow
colormap(h1,VzCmap);
caxis(vCrange);
colorbar
axis equal tight
colorbar
hold on;
h2=axes;
imagesc(BBV(:,:,2)); % down flow
% BBVdn=(LWMEAN(BBV(:,:,2),5));
% imagesc(BBVdn); % up flow
alpha(h2,1*double(abs(BBV(:,:,2))>1))
colormap(h2,VzCmap);
caxis(vCrange);
axis equal tight
colorbar
axis off
linkaxes([h1,h2]);
title ('V')
saveas(fig1,[FilePath,'BBV11-20.png'],'png');
saveas(fig1,[FilePath,'BBV11-20.fig'],'fig');
%% plot Vz
fig2=figure;
h1=axes;
imagesc(BBVz(:,:,1)); % up flow
colormap(h1,VzCmap);
caxis(vCrange);
colorbar
axis equal tight
colorbar
hold on;
h2=axes;
imagesc(BBVz(:,:,2)); % down flow
alpha(h2,1*double(abs(BBVz(:,:,2))>0.2))
colormap(h2,VzCmap);
caxis(vCrange);
axis equal tight
colorbar
axis off
linkaxes([h1,h2]);
title ('Vz')
saveas(fig2,[FilePath,'BBVz.png'],'png');
saveas(fig2,[FilePath,'BBVz.fig'],'fig');
%% plot directional flow map
Crange=[-5 5];
fig3=figure;
h1=axes;
imagesc(BB(:,:,1)); % up flow
% BBup=abs(LWMEAN(BB(:,:,1),3)).^(1/2)*(-1);
% imagesc(BBup); % up flow
colormap(h1,VzCmap);
caxis(Crange);
% colorbar
axis equal tight
hold on;
h2=axes;
imagesc(BB(:,:,2)); % down flow
% BBdn=abs(LWMEAN(BB(:,:,2),3)).^(1/2);
% imagesc(BBdn); % up flow
alpha(h2,1*double(abs(BB(:,:,2))>0.2))
colormap(h2,VzCmap);
caxis(Crange);
axis equal tight
% colorbar
axis off
linkaxes([h1,h2]);
title ('BB')
saveas(fig3,[FilePath,'BB.png'],'png');
saveas(fig3,[FilePath,'BB.fig'],'fig');
%% save data
save([FilePath, 'BB.mat'],'BB','BBV','BBVz','PRSSinfo');%,'Coor');