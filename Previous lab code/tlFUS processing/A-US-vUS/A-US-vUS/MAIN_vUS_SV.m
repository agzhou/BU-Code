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
% output:
    % PDI: Power Doppler based fUS, [nz,nx,3], 3: [up, down, all]
    % Ms: static component fraction, [nz,nx]
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
    % [Mf, Vx, Vz, V, pVz ,R, Ms, CR, GGf]=sIQ2vUS_SV_GPU(sIQ, PRSSinfo)
        % GG = sIQ2GG(sIQ, PRSSinfo)
        % RotCtr = FindCOR(GG)
        % [Vz, Tvz]=GG2Vz(GG, PRSSinfo, nItp)
        % [Vz,Vx,pVz,Ms,Mf,R, GGf]=GG2vUS(GG, Vz0, Ms0, MfR0, PRSSinfo)
%clear all; clc
defaultpath='Z:\US-DATA\PROJ-D-vUS\Phantom-Validation-20181016-FlowRBC\DATA-FlowRBC-20181016-AngledX\';
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
%% Use GPU calculation or not
useGPU = questdlg('Use GPU for data processing?', 'Select', ...
    'YES', 'NO', 'Cancel', 'Cancel');
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
     '3', '1000','25',...
     '1', '125','100', num2str(P.vSound), '16.625'};
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
PRSSinfo.MpVz=0; % no velocity distribution
%% vUS calculation
for iCP=startCP:startCP+nCP-1
    for iRPT=startRPT:startRPT+nRPT-1
        tic;
        clear V;
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
        [PDI]=sIQ2PDI(sIQ);  % PDI processing
        disp(['Color Doppler Processing - ', datestr(datetime('now'))]);
        Vcz0=(ColorDoppler(sIQ,PRSSinfo)); % color Doppler, all frequency
        disp(['g1-based vUS Processing - ', datestr(datetime('now'))]);
        if strcmp(useGPU, 'YES')
            [Mf, Vx, Vz, V, pVz, R, Ms, CR, GGf]=sIQ2vUS_SV_GPU(sIQ, PRSSinfo);
        else
            [Mf, Vx, Vz, V, pVz, R, Ms, CR, GGf,sumGG0,sumGGV0,sumGGFV0,sumGGV0_]=sIQ2vUS_SV(sIQ, PRSSinfo);
        end
        Vcz=imresize(Vcz0, [nz,nx]*PRSSinfo.rfnScale,'bilinear').*CR;
        PDI=imresize(PDI, [nz,nx]*PRSSinfo.rfnScale,'bilinear').*CR;
        save([FilePath,'vUS',iFileName(3:end)],'-v7.3','Ms','Mf','Vx','Vz','V','R','pVz','Vcz','PDI','PRSSinfo','P','sumGG0','sumGGV0','sumGGFV0','sumGGV0_');
        disp(['Results are saved! - ', datestr(datetime('now'))]);
    end
end
%% mean value
% %Transverse ROI
% Vroi=abs(V(54:57,:));
% %Angled ROI
[CX, CY, C]=improfile(V,[8 247], [13 115]);
for iP=1:240
    Vroi(:,:,iP)=abs(V(floor(CY(iP))+[-3:3],floor(CX(iP))));
end  
Vmean=mean(Vroi(:))
Vstd=std(Vroi(:))
%% figure plot
[VzCmap,VzCmapDn, VzCmapUp, pdiCmapUp, PhtmCmap]=Colormaps_fUS;
Coor.x=PRSSinfo.xCoor;
Coor.z=PRSSinfo.zCoor;
Fig=figure;
set(Fig,'Position',[400 400 1700 350])
subplot(1,3,1)
h1=imagesc(Coor.x,Coor.z,abs(V)); 
colormap(PhtmCmap);
caxis([0 30]);
colorbar
axis equal tight;
xlabel('x [mm]')
ylabel('z [mm]')
title('vUS-V [mm/s]')

subplot(1,3,2)
h2=imagesc(Coor.x,Coor.z,abs(Vz)); 
colormap(PhtmCmap);
caxis([0 30]);
colorbar
axis equal tight;
xlabel('x [mm]')
ylabel('z [mm]')
title('vUS-Vz [mm/s]')

subplot(1,3,3)
h3=imagesc(Coor.x,Coor.z,abs(Vcz)); 
colormap(PhtmCmap);
caxis([0 30]);
colorbar
axis equal tight;
xlabel('x [mm]')
ylabel('z [mm]')
title('Color Doppler-Vz [mm/s]')
