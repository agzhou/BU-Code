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
    % [Mf, Vx, Vz, V, pVz ,R, Ms, CR, GGf]=sIQ2vUS_NPDV_GPU(sIQ, PRSSinfo)
        % GG = sIQ2GG(sIQ, PRSSinfo)
        % RotCtr = FindCOR(GG)
        % [Vz, Tvz]=GG2Vz(GG, PRSSinfo, nItp)
        % [Vz,Vx,pVz,Ms,Mf,R,GGf]=GG2vUS_GPU(GG, Vz0, Ms0, MfR0, PRSSinfo)
clear all; clc
defaultpath='projectnb/npbfus/NonStash/DATA/tl_fUS_data/';
addpath('./SubFunctions') % Path on JTOPTICS
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
InVivoOrPhantom = questdlg('InVivo or Phantom', ...
    'Select', ...
    'InVivo', 'Phantom', 'Cancel', 'Cancel');
if strcmp(InVivoOrPhantom, 'InVivo')
    PRSSinfo.inVivo=1;
else
    PRSSinfo.inVivo=0;
end
%% Use GPU calculation or not
% useGPU = questdlg('Use GPU for data processing?', 'Select', ...
%     'YES', 'NO', 'Cancel', 'Cancel');
% if strcmp(useGPU, 'YES')
%     PRSSinfo.useGPU=1;
% else
%     PRSSinfo.useGPU=0;
% end
PRSSinfo.useGPU=1;
%% Load IQ data DAQ information
myFile=matfile([FilePath,FileName]);
P=myFile.P;
PRSSinfo.rfnScale=1;
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
covB=ones(3,3); covB(3,3)=9; covB=covB/sum(covB(:));
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
        PRSSinfo.Dim=[nz,nx,nt];
        clear IQ
        disp(['Power Doppler Processing - ', datestr(datetime('now'))]);
        [PDI0]=sIQ2PDI(sIQ);  % PDI processing
        [PDIHHP0]=sIQ2PDI(sIQHHP);  % PDI processing
        PDI=zeros(nz*PRSSinfo.rfnScale,nx*PRSSinfo.rfnScale,3);        
        for iD=1:3
            PDI(:,:,iD)=imresize(PDI0(:,:,iD),[nz,nx]*PRSSinfo.rfnScale,'bilinear');
            PDIHHP(:,:,iD)=imresize(PDIHHP0(:,:,iD),[nz,nx]*PRSSinfo.rfnScale,'bilinear');
        end
        eqNoise=imresize(eqNoise,[nz,nx]*PRSSinfo.rfnScale,'bilinear');
        disp(['Color Doppler Processing - ', datestr(datetime('now'))]);        
        %% g1fUS data processing
        disp(['g1-based fUS Processing - ', datestr(datetime('now'))]);                
        if PRSSinfo.inVivo==1
            Tau1=2;
            Tau2=7;
            % negative and positive frequency GG
            disp('npIQ to GG to CBF indices...');
            fIQ=sysNoiseRemove(sIQ,PRSSinfo.rFrame); % ** filtered (denoised) IQ **
            [nz,nx,nt]=size(sIQ);
            npGG=zeros(nz,nx,PRSSinfo.g1nTau); % ** g1 with negative and positive frequencies separated **
            unNorm_npGG=npGG;                  % ** un-normalized g1 **
            Vcz=zeros(nz*PRSSinfo.rfnScale,nx*PRSSinfo.rfnScale,3); % ** initialize the z velocity map **
            Vcz(:,:,3)=(ColorDoppler(sIQ,PRSSinfo)); % color Doppler, all frequency

            % ** go through the negative and positive frequencies **
            for iNP=1:2
                iFIQ=zeros(size(fIQ));
                iFIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP)=fIQ(:,:,(floor(nt/2))*(iNP-1)+1:floor(nt/2)*iNP); % Slice the Fourier-Transformed IQ data
                iIQ=(ifft(iFIQ,nt,3)); % ** inverse FFT across frames to get the filtered IQ with only the negative or positive frequencies **
                Vcz(:,:,iNP)=(ColorDoppler_NP(iIQ,PRSSinfo)); % improved color Doppler, all frequency            
                %% normalized g1 calculation
                npGG=sIQ2GG(iIQ, PRSSinfo); % g1 of p or n frequency signal    
                % g1(1) adjustment
                GG2=reshape(npGG,[PRSSinfo.Dim(1)*PRSSinfo.Dim(2),PRSSinfo.g1nTau]);
                ggCR1=(abs(abs(GG2(:,1))-abs(GG2(:,2)))>2*abs((abs(GG2(:,2))-abs(GG2(:,3)))));
                ggCR2=(abs(GG2(:,1))>0.55).*(abs(GG2(:,2))<0.25).*(abs(GG2(:,2))<abs(GG2(:,3)));
                ggCR0=((ggCR1+ggCR2)>0); % modified by Bingxue Liu;
                %GG2(:,1)=(1-ggCR0).*GG2(:,1)+ggCR0.*(GG2(:,2)+(abs(real(GG2(:,2)-GG2(:,3)))+1i*abs(imag(GG2(:,2)-GG2(:,3))))*1.5);
                GG2temp(:,1)=(1-ggCR0).*GG2(:,1)+ggCR0.*(GG2(:,2)+(abs(real(GG2(:,1)-GG2(:,2)))+1i*(imag(GG2(:,2)-GG2(:,3))))*1);% abs; real GG1-GG2
                GG2temp(find(abs(GG2temp)>1)) = GG2(find(abs(GG2temp)>1),1); %
                GG2(:,1) = GG2temp; % modified by Bingxue Liu
                npGG=reshape(GG2,[PRSSinfo.Dim(1),PRSSinfo.Dim(2),PRSSinfo.g1nTau]);
                % GG to CBFspeed Index and CBV index
                for itau=1:PRSSinfo.g1nTau
                    npGG(:,:,itau)=convn(npGG(:,:,itau),covB,'same');
                end
                GGsmth = smoothdata(npGG, 3, 'sgolay', 9);
                CBF_speedInd(:,:,iNP) = abs(squeeze(sqrt(log(abs(GGsmth(:,:,Tau1))./abs((GGsmth(:,:,Tau2))))))/(sqrt(Tau2^2-Tau1^2)*0.1));
                CBVind(:,:,iNP) = squeeze(mean(abs(GGsmth(:,:,1)),3)./(1-mean(abs(GGsmth(:,:,1)),3)));
                
                %% no normalized g1 calculation
                unNorm_npGG=log(abs(sIQ2GG_unNorm_GPU(iIQ, PRSSinfo))); % g1 of p or n frequency signal                
                % g1(1) adjustment
                GG2=reshape(unNorm_npGG,[PRSSinfo.Dim(1)*PRSSinfo.Dim(2),PRSSinfo.g1nTau]);
                ggCR1=(abs(abs(GG2(:,1))-abs(GG2(:,2)))>2*abs((abs(GG2(:,2))-abs(GG2(:,3)))));
                ggCR2=(abs(GG2(:,1))>0.55).*(abs(GG2(:,2))<0.25).*(abs(GG2(:,2))<abs(GG2(:,3)));
                ggCR0=((ggCR1+ggCR2)>0); % modified by Bingxue Liu;
                %GG2(:,1)=(1-ggCR0).*GG2(:,1)+ggCR0.*(GG2(:,2)+(abs(real(GG2(:,2)-GG2(:,3)))+1i*abs(imag(GG2(:,2)-GG2(:,3))))*1.5);
                GG2temp(:,1)=(1-ggCR0).*GG2(:,1)+ggCR0.*(GG2(:,2)+(abs(real(GG2(:,1)-GG2(:,2)))+1i*(imag(GG2(:,2)-GG2(:,3))))*1);% abs; real GG1-GG2
                GG2temp(find(abs(GG2temp)>1)) = GG2(find(abs(GG2temp)>1),1); %
                GG2(:,1) = GG2temp; % modified by Bingxue Liu
                unNorm_npGG=reshape(GG2,[PRSSinfo.Dim(1),PRSSinfo.Dim(2),PRSSinfo.g1nTau]);
                % g1(1) normalize
%                 unNorm_npGG=unNorm_npGG-prctile(abs(unNorm_npGG(:)),0.2);
                for itau=1:PRSSinfo.g1nTau
                    unNorm_npGG(:,:,itau)=convn(unNorm_npGG(:,:,itau),covB,'same');
                end
                unNorm_npGG = smoothdata(unNorm_npGG, 3, 'sgolay', 9);
                unNorm_npGG=unNorm_npGG-min(abs(unNorm_npGG(:)));
                unNorm_npGG=unNorm_npGG/(max(abs(unNorm_npGG(:))));
%                 unNorm_npGG=unNorm_npGG/(prctile(abs(unNorm_npGG(:)),99.9));
%                 g11=abs(unNorm_npGG(:,:,1));
%                 [indZ,indX]=find(g11>=1);
%                 for iPix=1:length(indZ)
%                     unNorm_npGG(indZ(iPix),indX(iPix),:)=unNorm_npGG(indZ(iPix),indX(iPix),:)/max(g11(:));
%                 end                                
%                 unNorm_CBF_speedInd(:,:,iNP) = abs(squeeze(sqrt(log(abs(unNorm_npGG(:,:,Tau1))./abs(unNorm_npGG(:,:,Tau2))))))/(sqrt(Tau2^2-Tau1^2)*0.1);
                unNorm_CBF_speedInd(:,:,iNP) = abs(squeeze(sqrt((unNorm_npGG(:,:,Tau1))-(unNorm_npGG(:,:,Tau2)))))/(sqrt(Tau2^2-Tau1^2)*0.1);
                unNorm_CBVind(:,:,iNP) = squeeze(mean(abs(unNorm_npGG(:,:,1)),3)./(1-mean(abs(unNorm_npGG(:,:,1)),3)));
            end
            disp('sIQ to GG to CBF indices...');
            GG = sIQ2GG_GPU(sIQ, PRSSinfo);
            for itau=1:PRSSinfo.g1nTau
                GG(:,:,itau)=convn(GG(:,:,itau),covB,'same');
            end
            GG = smoothdata(GG, 3, 'sgolay', 9);
            CBF_speedInd(:,:,3) = abs(squeeze(sqrt(log(abs(GG(:,:,Tau1))./abs(GG(:,:,Tau2)))))/(sqrt(Tau2^2-Tau1^2)*0.1));
            CBVind(:,:,3) = squeeze(mean(abs(GG(:,:,1)),3)./(1-mean(abs(GG(:,:,1)),3)));
            
            unNorm_npGG=log(abs(sIQ2GG_unNorm_GPU(sIQ, PRSSinfo))); % g1 of p or n frequency signal
            for itau=1:PRSSinfo.g1nTau
                unNorm_npGG(:,:,itau)=convn(unNorm_npGG(:,:,itau),covB,'same');
            end
%             unNorm_npGG=unNorm_npGG-prctile(abs(unNorm_npGG(:)),1);
            unNorm_npGG=unNorm_npGG-min(abs(unNorm_npGG(:)));   
            unNorm_npGG=unNorm_npGG/(max(abs(unNorm_npGG(:))));
%             unNorm_npGG=unNorm_npGG/(prctile(abs(unNorm_npGG(:)),99.9));
%             g11=abs(unNorm_npGG(:,:,1));
%             [indZ,indX]=find(g11>=1);
%             for iPix=1:length(indZ)
%                 unNorm_npGG(indZ(iPix),indX(iPix),:)=unNorm_npGG(indZ(iPix),indX(iPix),:)/max(g11(:));
%             end            
            unNorm_npGG = smoothdata(unNorm_npGG, 3, 'sgolay', 9);
            unNorm_CBF_speedInd(:,:,3) = abs(squeeze(sqrt((unNorm_npGG(:,:,Tau1))-(unNorm_npGG(:,:,Tau2)))))/(sqrt(Tau2^2-Tau1^2)*0.1);
            unNorm_CBVind(:,:,3) = squeeze(mean(abs(unNorm_npGG(:,:,1)),3)./(1-mean(abs(unNorm_npGG(:,:,1)),3)));
            
            disp('CBF and CBV indices calculated ...');
            save([FilePath, 'g1fUS',iFileName(3:end-4),'.mat'],'-v7.3','CBF_speedInd','CBVind','unNorm_CBF_speedInd','unNorm_CBVind',...
                'PDI', 'Vcz','eqNoise','P');
        else % phantom data processing
            Vcz=zeros(nz*PRSSinfo.rfnScale,nx*PRSSinfo.rfnScale);
            Vcz(:,:)=(ColorDoppler(sIQ,PRSSinfo)); % color Doppler, all frequency
            % g1fUS
            Tau1=2;
            Tau2=12;
            [nz,nx,nt]=size(sIQ);
            
            disp('sIQ to GG to CBF indices...');
%             sIQn=sIQ+max(abs(sIQ(:)))/2*(rand(size(sIQ))-0.5); % add noise to signal for CBV
            GG = sIQ2GG_GPU(sIQ, PRSSinfo);
            GG=convn(GG,covB,'same');
            CBF_speedInd = abs(squeeze(sqrt(log(abs(GG(:,:,Tau1))./abs(GG(:,:,Tau2)))))/(sqrt(Tau2^2-Tau1^2)*0.1));
            CBVind(:,:) = squeeze(mean(abs(GG(:,:,1:2)),3)./(1-mean(abs(GG(:,:,1:2)),3)));            
           
            disp('CBF and CBV indices calculated ...');
            save([FilePath, 'g1fUS',iFileName(3:end-4),'.mat'],'-v7.3','CBF_speedInd','CBVind', 'PDI', 'Vcz','P');            
            %% plot the normalized g1-based CBF indices
            %         speedAxisRange=[prctile(CBF_speedInd(:),1) prctile(CBF_speedInd(:),99.5)];
            %         CBVAxisRange=[prctile(CBVind(:),12) prctile(CBVind(:),99)];
            speedAxisRange=[0.05 1];
            CBVAxisRange=[0.1 5];
            Fig=figure('visible','on');
            set(Fig,'Position',[100 100 1000 700]);
            ax(1)=subplot(2,2,1);
            imagesc(PDI(:,:,3));
            colormap(ax(1),hot)
            colorbar;
%             caxis(speedAxisRange)
            title('PDI')
            
            ax(2)=subplot(2,2,2);
            imagesc(Vcz(:,:));
            colormap(ax(2),VzCmap)
            caxis([-15 15])
            colorbar;
            title('Vcz')
            
            ax(3)=subplot(2,2,3);
            imagesc(CBF_speedInd(:,:));
            caxis(speedAxisRange)
            colormap(ax(3),jet)
            colorbar;
            title('CBF-speed Index,allFreq')
            
            ax(4)=subplot(2,2,4);
            imagesc(CBVind(:,:));
            colormap(ax(4),jet)
            caxis(CBVAxisRange)
            colorbar;
            title('CBV Index, allFreq')
            
            saveas(Fig,[FilePath, 'Norm_g1fUS',iFileName(3:end-4),'.tif'],'tif');
            saveas(Fig,[FilePath, 'Norm_g1fUS',iFileName(3:end-4),'.fig'],'fig');
            close(Fig)
        end
        %% figure plot and saving
        if PRSSinfo.inVivo==1
            %% plot the normalized g1-based CBF indices
            CBVeqNoise=abs(CBVind).^0.5.*eqNoise.^0.8;
            %         speedAxisRange=[prctile(CBF_speedInd(:),1) prctile(CBF_speedInd(:),99.5)];
            %         CBVAxisRange=[prctile(CBVind(:),12) prctile(CBVind(:),99)];
            speedAxisRange=[0.1 3];
            CBVAxisRange=[0.1 4];
            Fig=figure('visible','on');
            set(Fig,'Position',[100 100 1600 700]);
            ax(1)=subplot(2,3,1);
            imagesc(CBF_speedInd(:,:,1));
            colormap(ax(1),gray)
            caxis(speedAxisRange)
            title('CBF-speed Index, posFreq')
            
            ax(2)=subplot(2,3,2);
            imagesc(CBF_speedInd(:,:,2));
            colormap(ax(2),gray)
            caxis(speedAxisRange)
            title('CBF-speed Index, negFreq')
            
            ax(3)=subplot(2,3,3);
            imagesc(CBF_speedInd(:,:,3));
            caxis(speedAxisRange)
            colormap(ax(3),gray)
            title('CBF-speed Index,allFreq')
            
            ax(4)=subplot(2,3,4);
            imagesc(CBVeqNoise(:,:,1));
            colormap(ax(4),gray)
            caxis(CBVAxisRange)
            title('CBV Index, posFreq')
            
            ax(5)=subplot(2,3,5);
            imagesc(CBVeqNoise(:,:,2));
            colormap(ax(5),gray)
            caxis(CBVAxisRange)
            title('CBV Index, negFreq')
            
            ax(6)=subplot(2,3,6);
            imagesc(CBVeqNoise(:,:,3));
            colormap(ax(6),gray)
            caxis(CBVAxisRange)
            title('CBV Index, allFreq')
            
            saveas(Fig,[FilePath, 'Norm_g1fUS',iFileName(3:end-4),'.tif'],'tif');
            saveas(Fig,[FilePath, 'Norm_g1fUS',iFileName(3:end-4),'.fig'],'fig');
            close(Fig)
            %% plot the unnormalized g1-based CBF indices
            CBVeqNoise=abs(unNorm_CBVind).^0.5./eqNoise.^0.5;
            %         speedAxisRange=[prctile(unNorm_CBF_speedInd(:),0.1) prctile(unNorm_CBF_speedInd(:),99.5)];
            %         CBVAxisRange=[prctile(unNorm_CBVind(:),8) prctile(unNorm_CBVind(:),98)];
            speedAxisRange=[0.05 0.5];
            CBVAxisRange=[1 3];
            Fig=figure('visible','on');
            set(Fig,'Position',[100 100 1600 700]);
            ax(1)=subplot(2,3,1);
            imagesc(unNorm_CBF_speedInd(:,:,1));
            colormap(ax(1),gray)
            caxis(speedAxisRange)
            title('CBF-speed Index, posFreq')
            
            ax(2)=subplot(2,3,2);
            imagesc(unNorm_CBF_speedInd(:,:,2));
            colormap(ax(2),gray)
            caxis(speedAxisRange)
            title('CBF-speed Index, negFreq')
            
            ax(3)=subplot(2,3,3);
            imagesc(unNorm_CBF_speedInd(:,:,3));
            caxis(speedAxisRange)
            colormap(ax(3),gray)
            title('CBF-speed Index,allFreq')
            
            ax(4)=subplot(2,3,4);
            imagesc(CBVeqNoise(:,:,1));
            colormap(ax(4),gray)
            %         colormap(ax(4),custColormap)
            caxis(CBVAxisRange)
            title('CBV Index, posFreq')
            
            ax(5)=subplot(2,3,5);
            imagesc(CBVeqNoise(:,:,2));
            colormap(ax(5),gray)
            %         colormap(ax(5),custColormap)
            caxis(CBVAxisRange)
            title('CBV Index, negFreq')
            
            ax(6)=subplot(2,3,6);
            imagesc(CBVeqNoise(:,:,3));
            colormap(ax(6),gray)
            %         colormap(ax(6),custColormap)
            caxis(CBVAxisRange)
            title('CBV Index, allFreq')
            
            saveas(Fig,[FilePath, 'noNorm_g1fUS',iFileName(3:end-4),'.tif'],'tif');
            saveas(Fig,[FilePath, 'noNorm_g1fUS',iFileName(3:end-4),'.fig'],'fig');
            disp(['Results are saved! - ', datestr(datetime('now'))]);
            close (Fig)
        end
    end
end
% %% figure plot
% [VzCmap]=Colormaps_fUS;
% if strcmp(useULMmsk, 'YES')
%     Coor.x=PRSSinfo.xCoor;
%     Coor.z=PRSSinfo.zCoor;
%     
%     Fig=figure;
%     set(Fig, 'Position',[300 400 1800 450]);
%     subplot(1,3,1)
%     Fuse2Images(V(:,:,1),V(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
%     title(['vUS, V [mm/s]']);
%     subplot(1,3,2)
%     Fuse2Images(Vz(:,:,1),Vz(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
%     title(['vUS, Vz [mm/s]']);
%     subplot(1,3,3)
%     Fuse2Images(Vcz(:,:,1),Vcz(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
%     title(['iCD, Vcz [mm/s]']);
% else
%     Coor.x=PRSSinfo.xCoor;
%     Coor.z=PRSSinfo.zCoor;
%     Fig=figure;
%     set(Fig, 'Position',[300 400 1800 450]);
%     subplot(1,3,1)
%     Fuse2Images(V(:,:,1),V(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
%     title(['vUS, V [mm/s]']);
%     subplot(1,3,2)
%     Fuse2Images(Vz(:,:,1),Vz(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
%     title(['vUS, Vz [mm/s]']);
%     subplot(1,3,3)
%     Fuse2Images(Vcz(:,:,1),Vcz(:,:,2),[-30 30],[-30 30],Coor.x,Coor.z,2.5);
%     title(['iCD, Vcz [mm/s]']);
%     %% plot PDI-based HSV velocity map
%     figure;
%     PLOTwtV(V,(PDIHHP).^0.4,Coor,[-30 30])
%     title('vUS-V [mm/s]');
% end
