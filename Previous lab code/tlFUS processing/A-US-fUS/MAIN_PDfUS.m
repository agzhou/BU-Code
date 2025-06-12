%% IQ to fUS data processing, GPU optional data processing
% input:
    % IQ data (nz,nx,nt)
    % PRSSinfo: data acquistion information, including
        % PRSSinfo.SVDrank: SVD rank [low high]
        % PRSSinfo.HPfC:  High pass filtering cutoff frequency, Hz
        % PRSSinfo.NEQ: do noise equalization? 0: no noise equalization; 1: apply noise equalization
        % PRSSinfo.rfnScale: spatial refind scale
% output:
    % PDI: Power Doppler based fUS, [nz,nx,3], 3: [up, down, all]
% Functions:
    % [sIQ, sIQHP, sIQHHP, eqNoise]=IQ2sIQ(IQ,PRSSinfo)
        % [sIQ, Noise]=SVDfilter(IQ,SignalRank)
    % [PDI]=sIQ2PDI(sIQ)
clear all;
addpath('.\SubFunctions') % sub function Path
load('D:\CODE\Mains\DAQParameters.mat');
defaultpath=DAQInfo.savepath;
[FileName,FilePath]=uigetfile(defaultpath);  % read data of a small part of the brain cortex (IQR matrix) 
fileInfo=strsplit(FileName(1:end-4),'-');
startCP0=regexp(fileInfo{7},'\d*','Match');
cpName=regexp(fileInfo{7},'\D*','Match');
if isempty(startCP0)
    startCP1='0';
else
    startCP1=startCP0{1};
end
%% Load IQ data DAQ information
myFile=matfile([FilePath,FileName]);
P=myFile.P;
%% data processing parameters
prompt={'startRPT', 'nRPT','File (CP) Start','number of files (CPs)', ...
    'Frame Rate','SVD rank low', 'SVD rank high','HPcut','XZ refine'};
name='CSD data processing';
defaultvalue={fileInfo{8},'1',startCP1,'1',...
     num2str(P.CCFR),'25', '300','25', '1'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
startRPT=str2num(numinput{1});
nRPT=str2num(numinput{2});
startCP=str2num(numinput{3});
nCP=str2num(numinput{4});
PRSSinfo.rFrame=str2num(numinput{5}); % sIQ frame rate, Hz
PRSSinfo.SVDrank=[str2num(numinput{6}),min(str2num(numinput{7}),P.numCCframes)];
PRSSinfo.HPfC=str2num(numinput{8});
PRSSinfo.rfnScale=str2num(numinput{9});

PRSSinfo.dzImg=P.dzImg/PRSSinfo.rfnScale;
PRSSinfo.dxImg=P.dxImg/PRSSinfo.rfnScale;
PRSSinfo.xCoor=interp(P.xCoor,PRSSinfo.rfnScale);
PRSSinfo.zCoor=interp(P.zCoor,PRSSinfo.rfnScale);
PRSSinfo.NEQ=0; % no noise equalization
[VzCmap,VzCmapDn, VzCmapUp, PhtmCmap, fUSCmapUp, fUSCmapDn]=Colormaps_fUS;
%% data processing
for iCP=startCP:startCP+nCP-1
    for iRPT=startRPT:startRPT+nRPT-1
        tic;
        iFileInfo=fileInfo;
        if nCP>1
            iFileInfo{7}=[cpName{1},num2str(iCP)];
        end    
        iFileInfo{8}=[num2str(iRPT)];
        iFileName=[strjoin(iFileInfo,'-'),'.mat'];
        disp(['Loading data: ',num2str(iFileName),', ', datestr(datetime('now'))]);
        load ([FilePath, iFileName]);
        %% Clutter rejection
        disp(['Clutter Rejection - ', datestr(datetime('now'))]);
        [sIQ, sIQHP, sIQHHP, eqNoise]=IQ2sIQ(IQ,PRSSinfo); % 0: no noise equalization
        [nz,nx,nt]=size(sIQ);
        clear IQ
        disp(['Power Doppler Processing - ', datestr(datetime('now'))]);
        [PDI0]=sIQ2PDI(sIQ);  % PDI processing
%         [PDIHP0]=sIQ2PDI(sIQHP);  % PDI processing
%         [PDIHHP0]=sIQ2PDI(sIQHHP);  % PDI processing
        if PRSSinfo.rfnScale>1
            for iD=1:3
                PDI(:,:,iD)=imresize(PDI0(:,:,iD),[nz,nx]*PRSSinfo.rfnScale,'bilinear');
%                 PDIHP(:,:,iD)=imresize(PDIHP0(:,:,iD),[nz,nx]*PRSSinfo.rfnScale,'bilinear');
%                 PDIHHP(:,:,iD)=imresize(PDIHHP0(:,:,iD),[nz,nx]*PRSSinfo.rfnScale,'bilinear');
            end
            eqNoise=imresize(eqNoise,[nz,nx]*PRSSinfo.rfnScale,'bilinear');
        else
            PDI=PDI0;
%             PDIHP=PDIHP0;
%             PDIHHP=PDIHHP0;
        end
        %% plot and save fig
        Coor.x=PRSSinfo.xCoor;
        Coor.z=PRSSinfo.zCoor;
        
        PDIdb=log10(PDI./eqNoise.^1.5);
%         PDIdb=log10(PDI);
        for iPN=1:3
            PDIdb(:,:,iPN)=(PDIdb(:,:,iPN)-min(min(PDIdb(:,:,iPN))))/((max(max(PDIdb(:,:,iPN)))-min(min(PDIdb(:,:,iPN))))*0.9);
        end
        Fig=figure(1);
        set(Fig,'Position',[300 300 900 600])
        h=gca;
        imagesc(Coor.x,Coor.z,PDIdb(:,:,3));
        caxis([0.2 1.2])
        colormap(h,hot);
        colorbar
        axis equal tight
        xlabel('x [mm]')
        ylabel('z [mm]')
        title(['PDI all frequency']);
        
        save([FilePath,'PDI',iFileName(3:end)],'-v7.3','PDI','eqNoise','PRSSinfo','P');
        saveas(Fig,[FilePath, 'PDI',iFileName(4:end-4),'.tif'],'tif');
%         saveas(Fig,[FilePath, 'PDI',iFileName(4:end-4),'.fig'],'fig');
        disp(['Results are saved! - ', datestr(datetime('now'))]);
    end
end