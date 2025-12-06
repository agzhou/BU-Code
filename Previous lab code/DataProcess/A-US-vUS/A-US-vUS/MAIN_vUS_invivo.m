%% IQ to vUS data processing, GPU-based data processing
% input:
    % IQ data (nz,nx,nt)
    % PRSSinfo: data acquistion information, including
        % PRSSinfo.FWHM: (X, Y, Z) spatial resolution, Full Width at Half Maximum of point spread function, m
        % PRSSinfo.rFrame: DAQ frame rate, Hz
        % PRSSinfo.f0: Transducer center frequency, Hz
        % PRSSinfo.C: Sound speed in the sample, m/s
        % PRSSinfo.g1nT: g1 calculation sample number
        % PRSSinfo.g1nTau: maximum number of time lag
        % PRSSinfo.SVDrank: SVD rank [low high]
        % PRSSinfo.HPfC:  High pass filtering cutoff frequency, Hz
        % PRSSinfo.NEQ: do noise equalization? 0: no noise equalization; 1: apply noise equalization
        % PRSSinfo.rfnScale: spatial refind scale
        % PRSSinfo.useMsk: 1: use ULM data as spatial mask; 0: no spatial mask
        % PRSSinfo.ulmMsk: ULM-based spatial constrain mask
            % [nz,nx,3], 1: up flow (positive frequency); 2 down flow (negative
            % frequency); 3: all flow 
            % ulmMsk=1 otherwise
% output:
    % PDI: Power Doppler based fUS, [nz,nx,3], 3: [up, down, all]
    % Ms: static component fraction, [nz,nx,2]
    % Mf: dynamic component fraction, [nz,nx,2], 2: [real,imag]
    % Vx: x-direction velocity component, [nz,nx], mm/s
    % Vz: axial-direction velocity component, [nz,nx], mm/s
    % V=sqrt(Vx.^2+Vz.^2), [nz,nx], mm/s
    % pVz: Vz distribution (sigma-Vz), [nz,nx]
    % R: fitting accuracy, [nz,nx]
    % GGf: gg fitting results, [nz,nx, nTau]
    % Vcz: Color Doppler axial velocity, [nz,nx], mm/s
% Functions:
    % [sIQ, sIQHP, sIQHHP, eqNoise]=IQ2sIQ(IQ,PRSSinfo)
        % [sIQ, Noise]=SVDfilter(IQ,SignalRank)
    % [PDI]=sIQ2PDI(sIQ)
    % [Vcz]=ColorDoppler(sIQ,PRSSinfo)
    % [Mf, Vx, Vz, V, pVz ,R, Ms, CR, GGf]=sIQ2vUS_NPDV_GPU(sIQ, PRSSinfo)
        % GG = sIQ2GG(sIQ, PRSSinfo)
        % RotCtr = FindCOR(GG)
        % [Vz, Tvz]=GG2Vz(GG, PRSSinfo, nItp)
        % [Vz,Vx,pVz,Ms,Mf,R,GGf]=GG2vUS_GPU(GG, Vz0, Ms0, MfR0, PRSSinfo)
%clear all; clc
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
%% Use ULM spatial mask or not
useULMmsk = questdlg('Use ULM mask for vUS spatial constrain?', ...
    'Select', ...
    'YES', 'NO', 'Cancel', 'Cancel');
if strcmp(useULMmsk, 'YES')
    PRSSinfo.useMsk=1;
    PRSSinfo.rfnScale=2;
    if ~exist([FilePath,'BB.mat'],'file')
        msgbox('Warning: BB.mat file is not existed in the data folder!')
        return;
    else
        load([FilePath,'\BB.mat']);
    end
    PRSSinfo.ulmMsk=[];
    for iBB=1:3
        PRSSinfo.ulmMsk(:,:,iBB)=(imresize(abs(BB(:,:,iBB)),PRSSinfo.rfnScale/5)>1.5); % use ULM as spatial constrain mask
    end
else
    PRSSinfo.useMsk=0;
    PRSSinfo.ulmMsk=1;
    PRSSinfo.rfnScale=1;
end
%% Use GPU calculation or not
useGPU = questdlg('Use GPU for data processing?', 'Select', ...
    'YES', 'NO', 'Cancel', 'Cancel');
if strcmp(useGPU, 'YES')
    PRSSinfo.useGPU=1;
else
    PRSSinfo.useGPU=0;
end
%% Load IQ data DAQ information
myFile=matfile([FilePath,FileName]);
P=myFile.P;
%% data processing parameters
prompt={'File (CP) Start','number of files (CPs)','startRPT', 'nRPT', ...
    'g1StartT',['g1Nt (nCC=',num2str(P.numCCframes)], 'g1Ntau','Frame Rate',...
    'SVD rank low', 'SVD rank high','HPcut',...
    'XZ refine','FWHM-x (um)','FWHM-z (um)', 'vSound (m/s)', 'fTransCenter [MHz]'};
name='CSD data processing';
defaultvalue={startCP1,'1',fileInfo{8},'1',...
     '1','1000','100',num2str(P.CCFR),...
     '25', '1000','25',...
     num2str(PRSSinfo.rfnScale), '125','100', num2str(P.vSound), num2str(P.frequency)};
numinput=inputdlg(prompt,name, 1, defaultvalue);
startCP=str2num(numinput{1});
nCP=str2num(numinput{2});
startRPT=str2num(numinput{3});
nRPT=str2num(numinput{4});
PRSSinfo.g1StartT=str2num(numinput{5});
PRSSinfo.g1nT=str2num(numinput{6});
PRSSinfo.g1nTau=str2num(numinput{7});
PRSSinfo.rFrame=str2num(numinput{8}); % sIQ frame rate, Hz
PRSSinfo.SVDrank=[str2num(numinput{9}),min(str2num(numinput{10}),P.numCCframes)];
PRSSinfo.HPfC=str2num(numinput{11});
PRSSinfo.rfnScale=str2num(numinput{12});
PRSSinfo.FWHM=[str2num(numinput{13}) str2num(numinput{14})]*1e-6;  % (X, Z) spatial resolution, Full Width at Half Maximum of point spread function, m
% PRSSinfo.FWHM=[125 100]*1e-6; % x,z spacial resolution, amplitude, for angled flow fitting, m            
PRSSinfo.C=str2num(numinput{15});                    % sound speed, m/s
PRSSinfo.f0=str2num(numinput{16})*1e6;               % Transducer center frequency, Hz

PRSSinfo.dzImg=P.dzImg/PRSSinfo.rfnScale;
PRSSinfo.dxImg=P.dxImg/PRSSinfo.rfnScale;
PRSSinfo.xCoor=interp(P.xCoor,PRSSinfo.rfnScale);
PRSSinfo.zCoor=interp(P.zCoor,PRSSinfo.rfnScale);
PRSSinfo.NEQ=0; % no noise equalization
PRSSinfo.MpVz=1; % maximum velocity distribution, sigma-Vz
%% data processing
[VzCmap]=Colormaps_fUS;
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
        [sIQ, sIQHP, sIQHHP, eqNoise]=IQ2sIQ(IQ(:,:,1:PRSSinfo.g1nT),PRSSinfo); % 0: no noise equalization
        [nz,nx,nt]=size(sIQ);
        clear IQ
        disp(['Power Doppler Processing - ', datestr(datetime('now'))]);
        [PDI0]=sIQ2PDI(sIQ);  % PDI processing
        [PDIHHP0]=sIQ2PDI(sIQHHP);  % PDI processing
        for iD=1:3
            PDI(:,:,iD)=imresize(PDI0(:,:,iD),[nz,nx]*PRSSinfo.rfnScale,'bilinear');
            PDIHHP(:,:,iD)=imresize(PDIHHP0(:,:,iD),[nz,nx]*PRSSinfo.rfnScale,'bilinear');
        end
        eqNoise=imresize(eqNoise,[nz,nx]*PRSSinfo.rfnScale,'bilinear');
%         disp(['Color Doppler Processing - ', datestr(datetime('now'))]);
%         Vcz0=(ColorDoppler(sIQ,PRSSinfo)); % color Doppler, all frequency
        disp(['g1-based vUS Processing - ', datestr(datetime('now'))]);
        if strcmp(useGPU, 'YES')
            [Mf, Vz, V, pVz, Vcz, R, CR, Vx, Ms, pnR,GGf]=sIQ2vUS_NPDV_GPU(sIQ, PRSSinfo);
        else
            [Mf, Vz, V, pVz, Vcz, R, CR, Vx, Ms, pnR,GGf,sumGG0,sumGGV0,sumGGFV0, GGfid]=sIQ2vUS_NPDV(sIQ, PRSSinfo);
        end
        %% plot and save fig
        Coor.x=PRSSinfo.xCoor;
        Coor.z=PRSSinfo.zCoor;
        
        Fig=figure('visible','off');
        set(Fig, 'Position',[300 400 1800 450]);
        subplot(1,3,1)
        Fuse2Images(V(:,:,1),V(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
        title(['vUS, V [mm/s]']);
        subplot(1,3,2)
        Fuse2Images(Vz(:,:,1),Vz(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
        title(['vUS, Vz [mm/s]']);
        subplot(1,3,3)
        Fuse2Images(Vcz(:,:,1),Vcz(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
        title(['iCD, Vcz [mm/s]']);
        
        if strcmp(useULMmsk, 'YES')
            save([FilePath,'vUSBB',iFileName(3:end)],'-v7.3','Ms','Mf','Vx','Vz','V','R','CR','pVz','Vcz','PDI','PDIHHP','eqNoise','BB', 'BBV','BBVz','PRSSinfo','P');
            saveas(Fig,[FilePath, 'vUSBB',iFileName(3:end-4),'.tif'],'tif');
            saveas(Fig,[FilePath, 'vUSBB',iFileName(3:end-4),'.fig'],'fig');
        else
            save([FilePath,'vUS',iFileName(3:end)],'-v7.3','Ms','Mf','Vx','Vz','V','R','CR','pVz','Vcz','PDI','PDIHHP','eqNoise','PRSSinfo','P','sumGG0','sumGGV0','sumGGFV0','GGf');
            saveas(Fig,[FilePath, 'vUS',iFileName(3:end-4),'.tif'],'tif');
            saveas(Fig,[FilePath, 'vUS',iFileName(3:end-4),'.fig'],'fig');
        end
        disp(['Results are saved! - ', datestr(datetime('now'))]);
    end
end
%% figure plot
[VzCmap]=Colormaps_fUS;
if strcmp(useULMmsk, 'YES')
    Coor.x=PRSSinfo.xCoor;
    Coor.z=PRSSinfo.zCoor;
    
    Fig=figure;
    set(Fig, 'Position',[300 400 1800 450]);
    subplot(1,3,1)
    Fuse2Images(V(:,:,1),V(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
    title(['vUS, V [mm/s]']);
    subplot(1,3,2)
    Fuse2Images(Vz(:,:,1),Vz(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
    title(['vUS, Vz [mm/s]']);
    subplot(1,3,3)
    Fuse2Images(Vcz(:,:,1),Vcz(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
    title(['iCD, Vcz [mm/s]']);
else
    Coor.x=PRSSinfo.xCoor;
    Coor.z=PRSSinfo.zCoor;
    Fig=figure;
    set(Fig, 'Position',[300 400 1800 450]);
    subplot(1,3,1)
    Fuse2Images(V(:,:,1),V(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
    title(['vUS, V [mm/s]']);
    subplot(1,3,2)
    Fuse2Images(Vz(:,:,1),Vz(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
    title(['vUS, Vz [mm/s]']);
    subplot(1,3,3)
    Fuse2Images(Vcz(:,:,1),Vcz(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
    title(['iCD, Vcz [mm/s]']);
    %% plot PDI-based HSV velocity map
    figure;
    PLOTwtV(V,(PDIHHP).^0.4,Coor,[-30 30])
    title('vUS-V [mm/s]');
    %% plot GGf
    slt = [39,187];%[23,163];%[25,172];
%     figure; 
%     subplot(311);plot(squeeze(GGf(slt(1),slt(2),:,1)));
%     subplot(312);plot(squeeze(GGf(slt(1),slt(2),:,2)));
%     subplot(313);plot(squeeze(GGf(slt(1),slt(2),:,3)));
    
    figure; 
    subplot(311); plot(abs(squeeze(GGf(slt(1),slt(2),:,1))));
    subplot(312);plot(abs(squeeze(GGf(slt(1),slt(2),:,2))));
    subplot(313);plot(abs(squeeze(GGf(slt(1),slt(2),:,3))));
    
    figure; 
    subplot(311);hold on; plot(abs(squeeze(GGfid(slt(1),slt(2),:,1))));
    subplot(312);hold on;plot(abs(squeeze(GGfid(slt(1),slt(2),:,2))));
    subplot(313);hold on;plot(abs(squeeze(GGfid(slt(1),slt(2),:,3))));

end
