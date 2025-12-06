%% load and plot Power Doppler results
clear all; clc
defaultpath='Z:\US-DATA\PROJ-D-vUS\invivo-vUSvsULM-20181113-Acute-ISO\DATA-CP6-vUSBB\vUS-Baseline-CP6\';
addpath('.\SubFunctions') % sub function Path
[FileName,FilePath]=uigetfile(defaultpath);  % read data of a small part of the brain cortex (IQR matrix)  
fileInfo=strsplit(FileName(1:end-4),'-');
startCP0=fileInfo{7}(3:end);
if isempty(startCP0)
    startCP0=0;
end
%% data processing parameters
prompt={'Start file (CP)','Number of files (CPs)','Start Repeat', 'Number of Repeats'};
name='File info';
defaultvalue={num2str(startCP0), '1',fileInfo{8},'1','5'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
startCP=str2num(numinput{1});
nCP=str2num(numinput{2});          % number of coronal planes
startRpt=str2num(numinput{3});
nRpt=str2num(numinput{4});          % number of repeat for each coronal plane
%% plot and save figures
[VzCmap,VzCmapDn, VzCmapUp, PhtmCmap, fUSCmapUp, fUSCmapDn]=Colormaps_fUS;
for iCP=startCP:startCP+nCP-1
    for iRpt=startRpt:startRpt+nRpt-1
        iFileInfo=fileInfo;
        if nCP>1
            iFileInfo{7}(3:end)=num2str(iCP);
        end            
        iFileInfo{8}=num2str(iRpt);
        iFileName=[strjoin(iFileInfo,'-'),'.mat'];
        load([FilePath,iFileName]);
        Coor.x=PRSSinfo.xCoor;
        Coor.z=PRSSinfo.zCoor;
        PDIdb=log10(PDI./eqNoise.^1);
%         PDIdb=log10(PDIHHP./eqNoise.^1);
        for iPN=1:3
            PDIdb(:,:,iPN)=(PDIdb(:,:,iPN)-min(min(PDIdb(:,:,iPN))))/((max(max(PDIdb(:,:,iPN)))-min(min(PDIdb(:,:,iPN))))*0.9);
        end
%         Fig=figure('visible','off');
        Fig=figure;
        set(Fig, 'Position',[300 400 1700 450]);
        subplot(1,3,1)
        h1=imagesc(Coor.x,Coor.z,PDIdb(:,:,1).^1.5*(1.5)); % up flow
        colormap(fUSCmapUp);
        caxis([0 1]);
        colorbar
        axis equal tight;
        xlabel('x [mm]')
        ylabel('z [mm]')
        title(['PDI f(positive)']);
        
        subplot(1,3,2)
        h=gca;
        h2=imagesc(Coor.x,Coor.z,PDIdb(:,:,2).^1.5*1.4); % down flow
        colormap(gca,fUSCmapDn);
        caxis([0 1]);
        colorbar
        hold off
        axis equal tight;
        xlabel('x [mm]')
        ylabel('z [mm]')
        title(['PDI f(negative)']);
        subplot(1,3,3)
        h=gca;
        imagesc(Coor.x,Coor.z,PDIdb(:,:,3));
        colormap(h,gray);
        colorbar
        axis equal tight
        xlabel('x [mm]')
        ylabel('z [mm]')
        title(['PDI all frequency']);
        
        saveas(Fig,[FilePath, 'PDI',iFileName(4:end-4),'.tif'],'tif');
        saveas(Fig,[FilePath, 'PDI',iFileName(4:end-4),'.fig'],'fig');
        disp(['Results are saved! - ', datestr(datetime('now'))]);
    end
end
%% positive frequency PDI overlapped on negative frequency PDI
% figure;
% ha1=axes;
% imagesc(Coor.x,Coor.z,PDIdb(:,:,1).^1.5*(1.5)); % up flow
% colormap(ha1,fUSCmapUp);
% caxis([0 1]);
% colorbar
% axis equal tight;
% xlabel('x [mm]')
% ylabel('z [mm]')
% 
% hold on
% ha2=axes;
% h2=imagesc(Coor.x,Coor.z,PDIdb(:,:,2).^1.5*1.4); % down flow
% AlphaMsk=(abs(PDIdb(:,:,2))/0.5).^5;
% AlphaMsk(AlphaMsk>1)=1;
% AlphaMsk(AlphaMsk<0.3)=0;
% set(h2,'AlphaData',AlphaMsk*0.6)
% colormap(VzCmap);
% colormap(ha2,fUSCmapDn);
% caxis([0 1]);
% colorbar
% hold off
% axis equal tight off;
% title(['PDI f_p overlapped on f_n']);

