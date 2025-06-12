 clear all; close all;
addpath D:\CODE\DataProcess\g1_based_fUS\A-US-fUS-g1\SubFunctions

%% I.load GG data
[FileName,FilePath]=uigetfile('V:\G1based fUS data\');
fileInfo=strsplit(FileName(1:end-4),'-');
myFile=matfile([FilePath,FileName]);
P=myFile.P;
prompt={'Start Repeat', 'Number of Repeats'};
name='File info'; 
defaultvalue={'1','325','1'};
numinput=inputdlg(prompt,name, 1, defaultvalue);
startRpt=str2num(numinput{1});
nRpt=str2num(numinput{2});          % number of repeat for each coronal plane

indSkipped=1; k = 1; 
for iRpt=startRpt:startRpt+nRpt-1
 iFileInfo=fileInfo;
 iFileInfo{8}=num2str(iRpt);
 iFileName=[strjoin(iFileInfo,'-'),'.mat'];
            %     load([FilePath,iFileName]);
 if exist([FilePath,iFileName],'file')
  myFile=matfile([FilePath,iFileName]);

    VV = myFile.V;
    VVz = myFile.Vz;
    VVx = myFile.Vx;
    VVcz = myFile.Vcz; 
    GG = myFile.GG;
    Ms = myFile.Ms;
    Mf = myFile.Mf;
    sumGG0 = myFile.sumGG0;% sumGG0_ means GG-Ms
    sumGGV0 = myFile.sumGGV0;
    sumGGFV0 = myFile.sumGGFV0;

  PDIHP = myFile.PDIHP;
  eqNoise = myFile.eqNoise;
  disp([iFileName,' was loaded!'])
  else
   disp([iFileName, ' skipped!'])
  SkipFile(indSkipped)=iRpt;
   indSkipped=indSkipped+1;  
 end

VV0(:,:,:,k) = VV;
VVz0(:,:,:,k)  = VVz;
VVx0(:,:,:,k)  = VVx;
VVcz0(:,:,:,k)  = VVcz;
GG0(:,:,:,k) = GG;
Ms0(:,:,:,k) = Ms;
Mf0(:,:,:,k) = Mf(:,:,:,1);%sqrt((Mf(:,:,:,1).^2+Mf(:,:,:,2).^2)/2); % MfR
sumGG(:,:,:,k) = sumGG0;
sumGGV(:,:,:,k) = sumGGV0;
sumGGFV(:,:,:,k) = sumGGFV0;
sumGGF(:,:,:,k) = sumGGV0.*sumGG0;

PDIHP0(:,:,:,k) = PDIHP;
eqNoise0(:,:,k) = eqNoise;
k = k+1;
end
eqNoise = mean(eqNoise0,3);
PDI = PDIHP0./eqNoise.^1.8;
[nz,nx,~,~] = size(VV0);
vmsk = ones([nz,nx]);
%% default trial defination 
trial.nBase = 25;
trial.n = 10;
trial.nStim = 5;
trial.nRecover = 25;
trial.nRest = 5; % included in nRecover
trial.nlength = trial.nStim + trial.nRecover;
% trial.parttern = Stim(trial.nBase-10+1:trial.nBase-10+(trial.nStim+trial.nRecover)); 
stim = zeros(nRpt,1);
for k = 1: trial.n % # of trials
stim(trial.nBase+(k-1)*trial.nlength+1: trial.nBase+trial.nStim+(k-1)*trial.nlength) = 1;
end
trial.stim = stim;
%% image background
ImgBG = mean(squeeze(PDI(:,:,3,:)),3).^0.35;
trialPDI = reshape(squeeze(PDI(:,:,3,trial.nBase-5+1: end-5)), [nz,nx,trial.nlength,trial.n]);
mtrialPDI = median(trialPDI,4);
%% Global HRF 
%hrf = hemodynamicResponse(1,[2 16 0.5 1 20 0]);
%  hrf = hemodynamicResponse(1,[1.5 10 0.5 1 20 0 16]);
hrf = hemodynamicResponse(0.1,[2 6 .5 1 10 0 12]);

Stim = zeros(trial.nlength*10,1);
Stim(5*10+1:5*10+trial.nStim*10,:)=1;
stimhrf0 = filter(hrf,1,Stim);
stimhrf = interp1(1:trial.nlength*10, stimhrf0, 1:10:trial.nlength*10, 'linear')';
% figure; subplot(211)
% plot(hrf);
% subplot(212)
% plot(0:0.1:trial.nlength-0.1, stimhrf0);
% hold on; plot(0: trial.nlength-1, stimhrf)
% % hold on;
% % plot(smooth(stimhrf,3));
% % hold on; plot(Stim);

%% GG
GG0_ = smoothdata(abs(GG0), 2, 'sgolay', 11);
%%
sumGG0 = squeeze(sqrt(abs(-log((GG0_(:,20,:,:))./(GG0_(:,2,:,:)))))); % (5,1);(5,2);(5,3);(5,4);(6,2);%(7,2);
sumGG = reshape(sumGG0, [nz, nx, 3, nRpt]);
%% GGV
sumGGV_ = sumGGV;
sumGGV = squeeze((Mf0)./(1-(Mf0)));

%% plot velocity map
 figure; subplot(121);imagesc(-squeeze(mean(VV0(:,:,1,:),4)));axis image; colorbar;title('ascending veins');axis off
subplot(122);imagesc(squeeze(mean(VV0(:,:,2,:),4)));axis image;colorbar;title('decending arteries');axis off 
%% data analysis
lib = {'G1(rCBFv)', 'G1(rCBV)','G1(rCBF)','G1(optimal)','PDI', 'G1(vUS)', 'ColorDoppler'} ;
%lib = {'G1(rCBFv)'};
for i = 1:numel(lib)
    datatype = lib{i};
    switch datatype
        case 'G1(rCBFv)'
            XX = sumGG;
        case 'G1(rCBV)'
            XX = sumGGV;
        case 'G1(rCBF)'
            XX = sumGGF;
        case 'G1(optimal)'

            XX = sumGGFV;
        case 'PDI'
            XX = PDI;
        case 'G1(vUS)'
            XX = VVx0;
        case 'ColorDoppler'
            XX = VVcz0;
    end
    %normalize data and preprocess
    trialXX = reshape(XX(:,:,1:2,trial.nBase-5+1: end-5), [nz,nx,2,trial.nlength,trial.n]);
    mtrialXX = median(trialXX,5);
    mtrialXX_bar = normalize(mtrialXX, 4, 'norm');
    mtrialXX_ratio = mtrialXX_bar./mean(mtrialXX_bar(:,:,:,1: trial.nRest),4);
    for iNP = 1:2
        %correlation map
        coefthld = 1.65; % 1.65 p<0.05; 2.35 p<0.01
        [actmap_XX(:,:,iNP), coefmap_XX(:,:,iNP)] = CoorCoeffMap(abs(squeeze(mtrialXX(:,:,iNP,:))), stimhrf', 0, coefthld);
%         actMsk(:,:,iNP) = calActivatedImage(abs(squeeze(mtrialXX(:,:,iNP,:))), trial);
%         actmap_XX(:,:,iNP) = coefmap_XX(:,:,iNP).*actMsk(:,:,iNP);
        %amplitude map
        window = [trial.nRest+1, trial.nRest+trial.nStim+1];            
        ampXX0(:,:,iNP) = max(mtrialXX_ratio(:,:,iNP,window(1):window(2)),[],4);
        ampXX(:,:,iNP) = ampXX0(:,:,iNP).*(actmap_XX(:,:,iNP)>0);
        %calculate CNR image
        [cnrXX(:,:,iNP), cnrXXnum(:,:,iNP), cnrXXden(:,:,iNP), cnrXX_final(:,:,iNP)] = calCNRImage(abs(squeeze(mtrialXX(:,:,iNP,:))), trial);
        %select roi
        exist MSKXX
        if ans ~= 1
            MSKXX = selectROI(actmap_XX(:,:,iNP)>0);
        end
       
        BW(:,:,iNP) = MSKXX;%.*(actmap_XX(:,:,iNP)>0);
        %count a/v pixels in actmap
        Npixel(:,iNP) = sum(BW(:,:,iNP),'all'); 
        %calculate roi time course
        roiXX(iNP,:) = squeeze(sum(sum(squeeze(XX(:,:,iNP,:)).*repmat(BW(:,:,iNP).*vmsk,[1, 1,nRpt]),1),2))/(Npixel(:,iNP)+1e-6);
    end
    roiXXall = abs(roiXX(1,:)).*Npixel(:,1)/sum(Npixel,2)+roiXX(2,:).*Npixel(:,2)/sum(Npixel,2);
    ROIXX = averagedTrials(roiXXall, trial);
    ROIXX.vein = averagedTrials(roiXX(1,:), trial);
    ROIXX.artery = averagedTrials(roiXX(2,:), trial);

    switch datatype
        case 'G1(rCBFv)'
            actmap_GG = actmap_XX; ampGG = ampXX0; CNR.GG = cnrXX_final;BWs.GG = BW;anpratio.GG = Npixel; ROIGG = ROIXX; coefmap_GG = coefmap_XX;
        case 'G1(rCBV)'
            actmap_GGV = actmap_XX; ampGGV = ampXX0; CNR.GGV = cnrXX_final;BWs.GGV = BW;anpratio.GGV = Npixel; ROIGGV = ROIXX; coefmap_GGV = coefmap_XX;
        case 'G1(rCBF)'
            actmap_GGF = actmap_XX; ampGGF = ampXX0; CNR.GGF = cnrXX_final;BWs.GGF = BW;anpratio.GGF = Npixel; ROIGGF = ROIXX; coefmap_GGF = coefmap_XX;
        case 'G1(optimal)'
            actmap_GGFV = actmap_XX; ampGGFV = ampXX0; CNR.GGFV = cnrXX_final;BWs.GGFV = BW;anpratio.GGFV = Npixel; ROIGGFV = ROIXX; coefmap_GGFV = coefmap_XX;
        case 'PDI'
            actmap_PDI = actmap_XX; ampPDI= ampXX0; CNR.PDI= cnrXX_final;BWs.PDI = BW;anpratio.PDI = Npixel; ROIPDI = ROIXX; coefmap_PDI = coefmap_XX;
        case 'G1(vUS)'
            actmap_VV = actmap_XX; ampVV= ampXX0; CNR.VV= cnrXX_final;BWs.VV = BW;anpratio.VV = Npixel; ROIVV = ROIXX; coefmap_VV = coefmap_XX;
        case 'ColorDoppler'
            actmap_VVcz = actmap_XX; ampVVcz= ampXX0; CNR.VVcz = cnrXX_final;BWs.VVcz = BW;anpratio.VVcz = Npixel; ROIVVcz = ROIXX; coefmap_VVcz = coefmap_XX;
    end  
end

%%
Fig = figure;
set(Fig, 'Position', [500 500 1200 300])
subplot(241);
imagesc(max(CNR.GG(:,:,:),[],3)); title('G1(rCBFv)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(245)
imagesc(max(CNR.GGV(:,:,:),[],3));title('G1(rCBV)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(246)
imagesc(max(CNR.GGF(:,:,:),[],3));title('G1(rCBF)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(247)
imagesc(max(CNR.GGFV(:,:,:),[],3));title('G1(optimal)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(244)
imagesc(max(CNR.PDI(:,:,:),[],3));title('PDI'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(242)
imagesc(max(CNR.VV(:,:,:),[],3));title('G1(vUS)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(243)
imagesc(max(CNR.VVcz(:,:,:),[],3));title('ColorDoppler'); axis image; caxis([0.5,2]);colorbar;axis off;
%%
BWall = sum((BWs.GG.*BWs.GGV),3)>0;%(sum((BWs.VV + BWs.PDI),3))>0;% + BWs.GGF + BWs.GGFV + BWs.VV + BWs.VVcz + BWs.PDI; sum((BWs.GG.*BWs.GGV),3)>0;
% BWall = (BWs.GG .* BWs.GGV .* BWs.GGF .* BWs.GGFV .* BWs.VV .* BWs.VVcz.*BWs.PDI);
ampGG2 = reshape(ampGG,[nx*nz,2]); ampGGV2 = reshape(ampGGV,[nx*nz,2]);
GGvsGGV(:,:,1) = [ampGG2(find(BWall>0),1), ampGGV2(find(BWall>0),1)];
GGvsGGV(:,:,2) = [ampGG2(find(BWall>0),2),ampGGV2(find(BWall>0),2)];
y1 = polyfit(GGvsGGV(:,1,1),GGvsGGV(:,2,1),1);
p1 = polyfit(GGvsGGV(:,1,1),GGvsGGV(:,2,1),1);
y1 = polyval(p1,GGvsGGV(:,1,1));
y2 = polyfit(GGvsGGV(:,1,2),GGvsGGV(:,2,2),1);
p2 = polyfit(GGvsGGV(:,1,2),GGvsGGV(:,2,2),1);
y2 = polyval(p2,GGvsGGV(:,1,2));
figure;
scatter(GGvsGGV(:,1,1),GGvsGGV(:,2,1),'b'); hold on; plot(GGvsGGV(:,1,1),y1,'b-','LineWidth',2); hold on;
scatter(GGvsGGV(:,1,2),GGvsGGV(:,2,2),'r'); hold on; plot(GGvsGGV(:,1,2),y2,'r-','LineWidth',2); 
axis equal; ylim([0.8,2.2]);xlim([0.8,2.2]);
xlabel('Velocity change');ylabel('Volume change')
clear GGvsGGV
%%
disp(['Activated pixels veins percentage (g1-CBFv) = ', num2str(anpratio.GG(1)/sum(anpratio.GG)*100),'%']);
disp(['Activated pixels veins percentage (g1-vUS) = ', num2str(anpratio.VV(1)/sum(anpratio.VV)*100),'%']);
disp(['Activated pixels veins percentage (ColorDoppler) = ', num2str(anpratio.VVcz(1)/sum(anpratio.VVcz)*100),'%']);
disp(['Activated pixels veins percentage (g1-CBV) = ', num2str(anpratio.GGV(1)/sum(anpratio.GGV)*100),'%']);
disp(['Activated pixels veins percentage (pdi) = ', num2str(anpratio.PDI(1)/sum(anpratio.PDI)*100),'%']);
disp(['Activated pixels veins percentage (g1-CBF) = ', num2str(anpratio.GGF(1)/sum(anpratio.GGF)*100),'%']);
disp(['Activated pixels veins percentage (g1-CBFV) = ', num2str(anpratio.GGFV(1)/sum(anpratio.GGFV)*100),'%']);
%% averaged time course compare

R1 = max(movmean(ROIGG.m, 1)/100.*movmean(ROIGGV.m, 1)/100-1)./(max(movmean(ROIGGV.m, 1)/100)-1);
R2 = max(movmean(ROIGG.m, 1)/100.*movmean(ROIPDI.m, 1)/100-1)./(max(movmean(ROIPDI.m, 1)/100)-1);
R3 = max(movmean(ROIVV.m, 1)/100.*movmean(ROIGGV.m, 1)/100-1)./(max(movmean(ROIGGV.m, 1)/100)-1);
R4 = max(movmean(ROIVV.m, 1)/100.*movmean(ROIPDI.m, 1)/100-1)./(max(movmean(ROIPDI.m, 1)/100)-1);
disp(['Ratio of Peak rCBF / Peak rCBV: (g1-CBFv*g1-CBV) = ', num2str(roundn(R1,-2))])
disp(['Ratio of Peak rCBF / Peak rCBV: (g1-CBFv*pdi) = ', num2str(roundn(R2,-2))])
disp(['Ratio of Peak rCBF / Peak rCBV: (g1-vUS*g1-CBV) = ', num2str(roundn(R3,-2))])
disp(['Ratio of Peak rCBF / Peak rCBV: (g1-vUS*pdi) = ', num2str(roundn(R4,-2))])

%%
Fig = figure; set(Fig,'Position',[400 600 1200 300])
subplot(131)
plot(ROIGG.m,'r'); 
hold on; plot(ROIVV.m,'b');
hold on; plot(ROIVVcz.m,'c');
plot(trial.nRest+1:trial.nRest+trial.nStim,ones(trial.nStim,1)*min(ROIPDI.m(:)),'-k','LineWidth',2);
legend({'g1-CBFv','g1-vUS','ColorDoppler-Vz'})
title(['Velocity change'])
xlabel('Time[s]');
ylabel('%');
ylim([90 150])
subplot(132)
hold on; plot(ROIGGV.m,'k');
hold on; plot(ROIPDI.m, 'k-.');
plot(trial.nRest+1:trial.nRest+trial.nStim,ones(trial.nStim,1)*min(ROIPDI.m(:)),'-k','LineWidth',2);
legend({'g1-CBV','pdi'})
title(['Volume change'])
xlabel('Time[s]');
ylabel('%');
ylim([90 150])
subplot(133)
hold on; plot(ROIGGV.m.*ROIGG.m/100,'r');hold on; 
hold on; plot(ROIPDI.m.*ROIGG.m/100,'r-.');
hold on; plot(ROIGGV.m.*ROIVV.m/100,'b');
hold on; plot(ROIPDI.m.*ROIVV.m/100,'b-.');
plot(trial.nRest+1:trial.nRest+trial.nStim,ones(trial.nStim,1)*min(ROIPDI.m(:)),'-k','LineWidth',2);
legend({'g1-CBFv*g1-CBV','g1-CBFv*pdi', 'g1-vUS*g1-CBV', 'g1-vUS*pdi'})
title(['Flow change'])
xlabel('Time[s]');
ylabel('%');
ylim([90 160])

%% a/v

