addpath D:\g1_based_fUS\A-US-fUS-g1\SubFunctions\
newSlect = 1;
[FileName,FilePath]=uigetfile('V:\G1based fUS data');
load([FilePath,'GG_invivo_results.mat'],'sumGG', 'sumGGV', 'sumGGFV', 'PDI', 'eqNoise', 'trial', 'dir', 'vmsk', 'nx', 'nz', 'ntau');
[~,~,nRpt] = size(sumGG);
tCoor = [1:trial.nlength]-trial.nRest;

bmsk = ones(size(dir));
bmsk(1:5,:) = 0;
bmsk(120:nz,:) = 0; % rough brain mask



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
%% 
%hrf = hemodynamicResponse(1,[2 16 0.5 1 20 0]);
 hrf = hemodynamicResponse(1,[1.5 10 0.5 1 20 0 16]);
% hrf = hemodynamicResponse(0.1,[2 6 .5 1 10 0 12]);

Stim = zeros(trial.nlength,1);
Stim(5+1:5+trial.nStim,:)=1;
stimhrf = filter(hrf,1,Stim);
% figure; subplot(211)
% plot(hrf);
% subplot(212)
% plot(0: trial.nlength-1, stimhrf)
% % hold on;
% % plot(smooth(stimhrf,3));
% hold on; plot(Stim);


%							defaults
%							(seconds)
%	p(1) - delay of response (relative to onset)	   6
%	p(2) - delay of undershoot (relative to onset)    16
%	p(3) - dispersion of response			   1
%	p(4) - dispersion of undershoot			   1
%	p(5) - ratio of response to undershoot		   6
%	p(6) - onset (seconds)				   0
%	p(7) - length of kernel (seconds)		  32
%% preprocess data

% get trial averaged data
trialsumGG = reshape(sumGG(:,:,trial.nBase-5+1: end-5), [nz,nx,trial.nlength,trial.n]);
mtrialsumGG = median(trialsumGG,4);

trialsumGGV = reshape(sumGGV(:,:,trial.nBase-5+1: end-5), [nz,nx,trial.nlength,trial.n]);
mtrialsumGGV = median(trialsumGGV,4);

trialsumGGFV = reshape(sumGGFV(:,:,trial.nBase-5+1: end-5), [nz,nx,trial.nlength,trial.n]);
mtrialsumGGFV = median(trialsumGGFV,4);

trialPDI = reshape(squeeze(PDI(:,:,3,trial.nBase-5+1: end-5)), [nz,nx,trial.nlength,trial.n]);
mtrialPDI = median(trialPDI,4);

%mtrialsumGGF = mtrialsumGG.*mtrialsumGGV;
sumGGF = sumGG.*sumGGV;
trialsumGGF = reshape(sumGGF(:,:,trial.nBase-5+1: end-5), [nz,nx,trial.nlength,trial.n]);
mtrialsumGGF = median(trialsumGGF,4);

% get normalized data
mtrialsumGG_bar = normalize(mtrialsumGG, 3, 'norm');
mtrialsumGGV_bar = normalize(mtrialsumGGV, 3, 'norm');
mtrialsumGGFV_bar = normalize(mtrialsumGGFV, 3, 'norm');
mtrialPDI_bar = normalize(mtrialPDI, 3, 'norm');
mtrialsumGGF_bar = normalize(mtrialsumGGF, 3, 'norm');

% get ralative change data
mtrialsumGG_ratio = mtrialsumGG_bar./mean(mtrialsumGG_bar(:,:,1: trial.nRest),3);
mtrialsumGGV_ratio = mtrialsumGGV_bar./mean(mtrialsumGGV_bar(:,:,1: trial.nRest),3);
mtrialsumGGFV_ratio = mtrialsumGGFV_bar./mean(mtrialsumGGFV_bar(:,:,1: trial.nRest),3);
mtrialPDI_ratio = mtrialPDI_bar./mean(mtrialPDI_bar(:,:,1: trial.nRest),3);
mtrialsumGGF_ratio = mtrialsumGGF_bar./mean(mtrialsumGGF_bar(:,:,1: trial.nRest),3);

%% Correlation map
% Stim1 = repmat(stimhrf', [1, trial.n]);
% [actmap_GG0, coefmap_GG0] = CoorCoeffMap(sumGG(:,:,trial.nBase-5+1: end-5).*vmsk.*bmsk, Stim1, 1);

coefthld = 1.65; % 1.65 p<0.05; 2.35 p<0.01

[actmap_GG0, coefmap_GG0] = CoorCoeffMap(mtrialsumGG, stimhrf', 0, coefthld);

[actmap_GGV0, coefmap_GGV0] = CoorCoeffMap(mtrialsumGGV_bar, stimhrf', 0, coefthld);

[actmap_GGFV0, coefmap_GGFV0] = CoorCoeffMap(mtrialsumGGFV_bar, stimhrf', 0, coefthld);

[actmap_PDI0, coefmap_PDI0] = CoorCoeffMap(mtrialPDI, stimhrf', 0, coefthld);

[actmap_GGF0, coefmap_GGF0] = CoorCoeffMap(mtrialsumGGF_bar, stimhrf', 0, coefthld);

actmap_GG0 = actmap_GG0.*vmsk.*bmsk;
actmap_GGV0 = actmap_GGV0.*vmsk.*bmsk;
actmap_GGFV0 = actmap_GGFV0.*vmsk.*bmsk;
actmap_GGF0 = actmap_GGF0.*vmsk.*bmsk;

%% ROIGG ROIGGV defined 
newSelect = 0;
ROItype = 'userdefined';%'actmap';
differentROI = 1;

if newSelect == 1
if differentROI == 1;
 switch ROItype
    case 'userdefined'
        if newSelect == 1
            MSKGG = selectROI(coefmap_GG0);
            MSKGGV = selectROI(coefmap_GGV0);
            MSKGGF = selectROI(coefmap_GGF0);          
        end
 end
end
else 
    [FileName,FilePath]=uigetfile('V:\G1based fUS data');
    load([FilePath,'GG_invivo_results.mat'],'MSKGG', 'MSKGGV', 'MSKGGF');
end

%% ROI averaged time course

newSelect = 1;
ROItype = 'actmap';%'actmap';%'rCBFv_based';%'userdefined';%
differentROI = 1;

if differentROI == 1;
 switch ROItype
    case 'actmap'
        BWGG0 = (actmap_GG0>0).*MSKGG;
        BWGGV0 = (actmap_GGV0>0).*MSKGGV;
        BWGGFV0 = (actmap_GGFV0>0).*MSKGGF;
        BWPDI0 = (actmap_PDI0>0).*MSKGGF;
        BWGGF0 = (actmap_GGF0>0).*MSKGGF;
    case 'userdefined'
        if newSelect == 1
            BWGG0 = selectROI(coefmap_GG0);
            BWGGV0 = selectROI(coefmap_GGV0);
            BWGGFV0 = selectROI(coefmap_GGFV0);
            BWPDI0 = selectROI(coefmap_PDI0);
            BWGGF0 = selectROI(coefmap_GGF0);
%             ROI0  = selectROI(coefmap_GG0);
%             BWGG0 = ROI0;
%             BWGGV0 = ROI0;
%             BWGGFV0 = ROI0;
%             BWPDI0 = ROI0;
%             BWGGF0 = ROI0;           
        end
 end
else 
    switch ROItype
        case "rCBFv_based"
        BWGG0 = (actmap_GG0>0).*MSKGG;
        BWGGV0 = (actmap_GG0>0).*MSKGG;
        BWGGFV0 = (actmap_GG0>0).*MSKGG;
        BWPDI0 = (actmap_GG0>0).*MSKGG;
        BWGGF0 = (actmap_GG0>0).*MSKGG;
        case "rCBV_based"
        BWGG0 = (actmap_GGV0>0).*MSKGGV;;
        BWGGV0 = (actmap_GGV0>0).*MSKGGV;;
        BWGGFV0 = (actmap_GGV0>0).*MSKGGV;;
        BWPDI0 = (actmap_GGV0>0).*MSKGGV;;
        BWGGF0 = (actmap_GGV0>0).*MSKGGV;;  
        case "rCBF_based"
        BWGG0 = (actmap_GGF0>0).*MSKGGF;
        BWGGV0 = (actmap_GGF0>0).*MSKGGF;
        BWGGFV0 = (actmap_GGF0>0).*MSKGGF;
        BWPDI0 = (actmap_GGF0>0).*MSKGGF;
        BWGGF0 = (actmap_GGF0>0).*MSKGGF;  
        case "userdefined"
            ROI0  = selectROI(coefmap_GG0.*MSKGG);
            BWGG0 = ROI0;
            BWGGV0 = ROI0;
            BWGGFV0 = ROI0;
            BWPDI0 = ROI0;
            BWGGF0 = ROI0;  
    end
end

roiGG = squeeze(sum(sum(sumGG.*repmat(BWGG0.*vmsk,[1, 1,nRpt]),1),2));
% Fig = figure;set(Fig,'Position',[400 400 800 170])                     
% yyaxis left;plot(-roiGG, 'b', 'LineWidth', 1);ylabel('G1-fUS(CBF)'); xlabel('Time[s]');                         
% yyaxis right; plot(stim);grid on;
roiGGV = squeeze(sum(sum(sumGGV.*repmat(BWGGF0.*vmsk,[1, 1,nRpt]),1),2));
% Fig = figure;set(Fig,'Position',[400 400 800 170]) 
% yyaxis left;plot(-roiGGV, 'b', 'LineWidth', 1);ylabel('G1-fUS(CBV)'); xlabel('Time[s]'); 
% yyaxis right; plot(stim);grid on;
roiGGFV = squeeze(sum(sum(sumGGFV.*repmat(BWGGFV0.*vmsk,[1, 1,nRpt]),1),2));
% Fig = figure;set(Fig,'Position',[400 400 800 170]) 
% yyaxis left;plot(-roiGGFV, 'b', 'LineWidth', 1);ylabel('G1-fUS(novel)'); xlabel('Time[s]'); 
% yyaxis right; plot(stim);grid on;
roiPDI = squeeze(sum(sum(squeeze(PDI(:,:,3,:)).*repmat(BWPDI0,[1, 1,nRpt]),1),2));
% yyaxis left;plot(roiPDI, 'b', 'LineWidth', 1);ylabel('PDI'); xlabel('Time[s]'); 
% yyaxis right; plot(stim);grid on;
roiGGF = squeeze(sum(sum(sumGGF.*repmat(BWGGF0.*vmsk,[1, 1,nRpt]),1),2));
% Fig = figure;set(Fig,'Position',[400 400 800 170]) 
% yyaxis left;plot(-roiGGFV, 'b', 'LineWidth', 1);ylabel('G1-fUS(novel)'); xlabel('Time[s]'); 
% yyaxis right; plot(stim);grid on;

% %% coefmap use ROI based HRF
% 
% % figure(2); imagesc(coefmap); axis image;colormap(jet);caxis([0,1]);
% HRF = calROIHRF(sumGG,trial,BWGG0);
% [actmap_GG, coefmap_GG] = CoorCoeffMap(mtrialsumGG.*vmsk, HRF, 0);
% 
% HRF = calROIHRF(sumGGV,trial,BWGGV0.*vmsk);
% [actmap_GGV, coefmap_GGV] = CoorCoeffMap(mtrialsumGGV.*vmsk, HRF, 0);
% 
% HRF = calROIHRF(sumGGFV,trial,BWGGFV0.*vmsk);
% [actmap_GGFV, coefmap_GGFV] = CoorCoeffMap(mtrialsumGGFV.*vmsk, HRF, 0);
% 
% HRF = calROIHRF(squeeze(PDI(:,:,3,:)),trial,BWPDI0.*vmsk);
% [actmap_PDI, coefmap_PDI] = CoorCoeffMap(mtrialPDI, HRF, 0); 
% 
% HRF = calROIHRF(sumGGF,trial,BWGGF0.*vmsk);
% [actmap_GGF, coefmap_GGF] = CoorCoeffMap(mtrialsumGGF.*vmsk, HRF, 0);

%% trial averaged time course GG vs. PDI

ROIGG = averagedTrials(roiGG, trial);

ROIGGV = averagedTrials(roiGGV, trial);

ROIGGFV = averagedTrials(roiGGFV, trial);

ROIPDI = averagedTrials(roiPDI, trial);

ROIGGF = averagedTrials(roiGGF, trial);

% %% plots
% % stimhrf correlation map 
% Fig = figure;set(Fig,'Position',[200 200 1000 800])
% subplot(421)
% imagesc(coefmap_GG0); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (g1-CBFv)');
% subplot(422)
% imagesc(actmap_GG0); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Activation map (Zscore>1.6)');
% subplot(423)
% imagesc(coefmap_GGV0); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (g1-CBV)');
% subplot(424)
% imagesc(actmap_GGV0); axis image; caxis([-1,1]); colormap('jet');colorbar;
% subplot(425)
% imagesc(coefmap_PDI0); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (PDI)');
% subplot(426)
% imagesc(actmap_PDI0); axis image; caxis([-1,1]); colormap('jet');colorbar;
% subplot(427)
% imagesc(coefmap_GGFV0); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (novel)');
% subplot(428)
% imagesc(actmap_GGFV0); axis image; caxis([-1,1]); colormap('jet');colorbar;
% 
% 
% % roihrf correlation map
% Fig = figure;set(Fig,'Position',[200 200 1000 800])
% subplot(421)
% imagesc(coefmap_GG); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (g1-CBFv)');
% subplot(422)
% imagesc(actmap_GG); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Activation map (Zscore>1.6)');
% subplot(423)
% imagesc(coefmap_GGV); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (g1-CBV)');
% subplot(424)
% imagesc(actmap_GGV); axis image; caxis([-1,1]); colormap('jet');colorbar;
% subplot(425)
% imagesc(coefmap_PDI); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (g1-PDI)');
% subplot(426)
% imagesc(actmap_PDI); axis image; caxis([-1,1]); colormap('jet');colorbar;
% subplot(427)
% imagesc(coefmap_GGFV); axis image; caxis([-1,1]); colormap('jet');colorbar;title('Correlation map (novel)');
% subplot(428)
% imagesc(actmap_GGFV); axis image; caxis([-1,1]); colormap('jet');colorbar;
% 
%% time course
stim = trial.stim;
Fig = figure;set(Fig,'Position',[400 400 700 700])
subplot(411)
yyaxis left;plot(roiGG, 'b', 'LineWidth', 1);ylabel('g1-fUS(CBFv)'); xlabel('Time[s]'); 
yyaxis right; plot(stim);grid on;
subplot(412)
yyaxis left;plot(roiGGV, 'b', 'LineWidth', 1);ylabel('g1-fUS(CBV)'); xlabel('Time[s]'); 
yyaxis right; plot(stim);grid on;
subplot(413)
yyaxis left;plot(roiPDI, 'b', 'LineWidth', 1);ylabel('PDI'); xlabel('Time[s]'); 
yyaxis right; plot(stim);grid on;
subplot(414)
yyaxis left;plot(roiGGFV, 'b', 'LineWidth', 1);ylabel('g1-fUS(novel)'); xlabel('Time[s]'); 
yyaxis right; plot(stim);grid on;

%% averaged time course compare
disp('Ratio of Peak rCBF / Peak rCBV: ')
R = max(movmean(ROIGGF.m, 1)/100-1)./(max(movmean(ROIGGV.m, 1)/100)-1)

Fig = figure; set(Fig,'Position',[400 600 450 300])
plot(ROIGG.m); 
hold on; plot(ROIGGV.m);
hold on; plot(ROIPDI.m);
hold on; plot(ROIGGFV.m);
hold on; plot(ROIGGF.m);
hold on; 
plot(trial.nRest+1:trial.nRest+trial.nStim,ones(trial.nStim,1)*min(ROIPDI.m(:)),'-k','LineWidth',2);
legend({'g1-CBFv','g1-CBV','pdi','novel','g1-CBF'})
title(['rCBFlow/rCBV = ', num2str(R)])
xlabel('Time[s]');
ylabel('%');
ylim([90 130])

%% arteries/veins divide

avdir = zeros(size(dir));
avdir(dir>0) = 1;
figure; 
subplot(121); imagesc(avdir); axis image; colormap('gray'); colorbar; title('vein g1 mask');
subplot(122); imagesc((1-avdir)); axis image; colormap('gray'); colorbar; title('artery g1 mask');

%% Directional Doppler based a/v mask
vPDm = (squeeze(mean(PDI(:,:,1,:)./(PDI(:,:,1,:)+PDI(:,:,2,:)),4)));
aPDm = (squeeze(mean(PDI(:,:,2,:)./(PDI(:,:,1,:)+PDI(:,:,2,:)),4)));
% 
% vPDm = log(squeeze(mean(PDI(:,:,1,:),4)));
% aPDm = log(squeeze(mean(PDI(:,:,2,:),4)));
figure; 
subplot(121); imagesc(vPDm); axis image; colormap('gray'); colorbar; title('vein Directional PD');
subplot(122); imagesc(aPDm); axis image; colormap('gray'); colorbar; title('artery Directional PD');
vPDmsk = vPDm>1.01*mean(vPDm(:));%(0.99*max(vPDm(:))+min(vPDm(:)))/2;
aPDmsk = aPDm>1.01*mean(aPDm(:));%(0.9*max(aPDm(:))+min(aPDm(:)))/2;
figure; 
subplot(121); imagesc(vPDmsk); axis image; colormap('gray'); colorbar; title('vein Directional PD mask');
subplot(122); imagesc(aPDmsk); axis image; colormap('gray'); colorbar; title('artery Directional PD mask');

%% 
% vMSK = vPDmsk.*avdir;
% aMSK = aPDmsk.*(1-avdir);

% vMSK = avdir;
% aMSK = (1-avdir);
% 
vMSK = vPDmsk;
aMSK = aPDmsk;

%% count a/v pixels in actmap
anp_GG0 = sum((actmap_GG0>0).*MSKGG.*aMSK,'all') 
vnp_GG0 = sum((actmap_GG0>0).*MSKGG.*vMSK,'all')

anp_GGV0 = sum((actmap_GGV0>0).*MSKGGV.*aMSK,'all')
vnp_GGV0 = sum((actmap_GGV0>0).*MSKGGV.*vMSK,'all')

%%
veinCmap_GG0 = coefmap_GG0.*vMSK;
arteriesCmap_GG0 = coefmap_GG0.*aMSK;
veinAmap_GG0 = actmap_GG0.*vMSK.*bmsk;
arteriesAmap_GG0 = actmap_GG0.*aMSK.*bmsk;
% figure; subplot(121); imagesc(veinCmap_GG0);caxis([-1,1]);axis image;title('veins')
% subplot(122); imagesc(arteriesCmap_GG0);caxis([-1,1]);axis image;title('arteries')
% figure; subplot(121); imagesc(veinAmap_GG0);caxis([-1,1]);axis image;title('veins')
% subplot(122); imagesc(arteriesAmap_GG0);caxis([-1,1]);axis image;title('arteries')

veinCmap_GGV0 = coefmap_GGV0.*vMSK;
arteriesCmap_GGV0 = coefmap_GGV0.*aMSK;
veinAmap_GGV0 = actmap_GGV0.*vMSK.*bmsk;
arteriesAmap_GGV0 = actmap_GGV0.*aMSK.*bmsk;
% figure; subplot(121); imagesc(veinCmap_GGV0);caxis([-1,1]);axis image;title('veins')
% subplot(122); imagesc(arteriesCmap_GGV0);caxis([-1,1]);axis image;title('arteries')
% figure; subplot(121); imagesc(veinAmap_GGV0);caxis([-1,1]);axis image;title('veins')
% subplot(122); imagesc(arteriesAmap_GGV0);caxis([-1,1]);axis image;title('arteries')

veinCmap_GGFV0 = coefmap_GGFV0.*vMSK;
arteriesCmap_GGFV0 = coefmap_GGFV0.*aMSK;
veinAmap_GGFV0 = actmap_GGFV0.*vMSK.*bmsk;
arteriesAmap_GGFV0 = actmap_GGFV0.*aMSK.*bmsk;
% figure; subplot(121); imagesc(veinCmap_GGV0);caxis([-1,1]);axis image;title('veins')
% subplot(122); imagesc(arteriesCmap_GGV0);caxis([-1,1]);axis image;title('arteries')
% figure; subplot(121); imagesc(veinAmap_GGFV0);caxis([-1,1]);axis image;title('veins')
% subplot(122); imagesc(arteriesAmap_GGFV0);caxis([-1,1]);axis image;title('arteries')

veinCmap_PDI0 = coefmap_PDI0.*vMSK;
arteriesCmap_PDI0 = coefmap_PDI0.*aMSK;
veinAmap_PDI0 = actmap_PDI0.*vMSK.*bmsk;
arteriesAmap_PDI0 = actmap_PDI0.*aMSK.*bmsk;
% figure; subplot(121); imagesc(veinCmap_GGV0);caxis([-1,1]);axis image;title('veins')
% subplot(122); imagesc(arteriesCmap_GGV0);caxis([-1,1]);axis image;title('arteries')
% figure; subplot(121); imagesc(veinAmap_PDI0);caxis([-1,1]);axis image;title('veins')
% subplot(122); imagesc(arteriesAmap_PDI0);caxis([-1,1]);axis image;title('arteries')

veinCmap_GGF0 = coefmap_GGF0.*vMSK;
arteriesCmap_GGF0 = coefmap_GGF0.*aMSK;
veinAmap_GGF0 = actmap_GGF0.*vMSK.*bmsk;
arteriesAmap_GGF0 = actmap_GGF0.*aMSK.*bmsk;
% figure; subplot(121); imagesc(veinCmap_GGV0);caxis([-1,1]);axis image;title('veins')
% subplot(122); imagesc(arteriesCmap_GGV0);caxis([-1,1]);axis image;title('arteries')
% figure; subplot(121); imagesc(veinAmap_GGF0);caxis([-1,1]);axis image;title('veins')
% subplot(122); imagesc(arteriesAmap_GGF0);caxis([-1,1]);axis image;title('arteries')

%% arteries/veins ROI averaged time course

newSelect = 1;
ROItype = 'actmap';%'rCBFv_based';%'rCBV_based';%'userdefined';%'actmap';%;%
differentROI = 1;

if differentROI == 1;
switch ROItype
    case 'actmap'
        vBWGG0 = (veinAmap_GG0>0).*MSKGG;
        aBWGG0 = (arteriesAmap_GG0>0).*MSKGG;
        vBWGGV0 = (veinAmap_GGV0>0).*MSKGGV;
        aBWGGV0 = (arteriesAmap_GGV0>0).*MSKGGV;
        vBWGGFV0 = (veinAmap_GGFV0>0).*MSKGGF;
        aBWGGFV0 = (arteriesAmap_GGFV0>0).*MSKGGF;
        vBWPDI0 = (veinAmap_PDI0>0).*MSKGGF;
        aBWPDI0 = (arteriesAmap_PDI0>0).*MSKGGF;
        vBWGGF0 = (veinAmap_GGF0>0).*MSKGGF;
        aBWGGF0 = (arteriesAmap_GGF0>0).*MSKGGF;
    case 'userdefined'
        if newSelect == 1
            vBWGG0 = selectROI(veinCmap_GG0);
            aBWGG0 = selectROI(arteriesCmap_GG0);
            vBWGGV0 = selectROI(veinCmap_GGV0);
            aBWGGV0 = selectROI(arteriesCmap_GGV0);
            vBWGGFV0 = selectROI(veinCmap_GGFV0);
            aBWGGFV0 = selectROI(arteriesCmap_GGFV0);
            vBWPDI0 = selectROI(veinCmap_PDI0);
            aBWPDI0 = selectROI(arteriesCmap_PDI0);
            vBWGGF0 = selectROI(veinCmap_GGF0);
            aBWGGF0 = selectROI(arteriesCmap_GGF0);

%             vROI = selectROI(veinCmap_GGV0);
%             aROI = selectROI(arteriesCmap_GGV0);
% 
%             vBWGG0 = vROI;
%             aBWGG0 = aROI;
%             vBWGGV0 = vROI;
%             aBWGGV0 = aROI;
%             vBWGGFV0 = vROI;
%             aBWGGFV0 = aROI;
%             vBWPDI0 = vROI;
%             aBWPDI0 = aROI;
%             vBWGGF0 = vROI;
%             aBWGGF0 = aROI;
        end
end
else
    switch ROItype
        case "rCBFv_based"
        vBWGG0 = (veinAmap_GG0>0).*MSKGG;
        aBWGG0 = (arteriesAmap_GG0>0).*MSKGG;
        vBWGGV0 = (veinAmap_GG0>0).*MSKGG;
        aBWGGV0 = (arteriesAmap_GG0>0).*MSKGG;
        vBWGGFV0 = (veinAmap_GG0>0).*MSKGG;
        aBWGGFV0 = (arteriesAmap_GG0>0).*MSKGG;
        vBWPDI0 = (veinAmap_GG0>0).*MSKGG;
        aBWPDI0 = (arteriesAmap_GG0>0).*MSKGG;
        vBWGGF0 = (veinAmap_GG0>0).*MSKGG;
        aBWGGF0 = (arteriesAmap_GG0>0).*MSKGG; 
        case "rCBV_based"
        vBWGG0 = (veinAmap_GGV0>0).*MSKGGV;
        aBWGG0 = (arteriesAmap_GGV0>0).*MSKGGV;
        vBWGGV0 = (veinAmap_GGV0>0).*MSKGGV;
        aBWGGV0 = (arteriesAmap_GGV0>0).*MSKGGV;
        vBWGGFV0 = (veinAmap_GGV0>0).*MSKGGV;
        aBWGGFV0 = (arteriesAmap_GGV0>0).*MSKGGV;
        vBWPDI0 = (veinAmap_GGV0>0).*MSKGGV;
        aBWPDI0 = (arteriesAmap_GGV0>0).*MSKGGV;
        vBWGGF0 = (veinAmap_GGV0>0).*MSKGGV;
        aBWGGF0 = (arteriesAmap_GGV0>0).*MSKGGV; 
        case "rCBF_based"
        vBWGG0 = (veinAmap_GGF0>0).*MSKGGF;
        aBWGG0 = (arteriesAmap_GGF0>0).*MSKGGF;
        vBWGGV0 = (veinAmap_GGF0>0).*MSKGGF;
        aBWGGV0 = (arteriesAmap_GGF0>0).*MSKGGF;
        vBWGGFV0 = (veinAmap_GGF0>0).*MSKGGF;
        aBWGGFV0 = (arteriesAmap_GGF0>0).*MSKGGF;
        vBWPDI0 = (veinAmap_GGF0>0).*MSKGGF;
        aBWPDI0 = (arteriesAmap_GGF0>0).*MSKGGF;
        vBWGGF0 = (veinAmap_GGF0>0).*MSKGGF;
        aBWGGF0 = (arteriesAmap_GGF0>0).*MSKGGF; 
        case "userdefined"
            vROI = selectROI(veinCmap_GG0);
            aROI = selectROI(arteriesCmap_GG0);

            vBWGG0 = vROI;
            aBWGG0 = aROI;
            vBWGGV0 = vROI;
            aBWGGV0 = aROI;
            vBWGGFV0 = vROI;
            aBWGGFV0 = aROI;
            vBWPDI0 = vROI;
            aBWPDI0 = aROI;
            vBWGGF0 = vROI;
            aBWGGF0 = aROI;
    end
end
%         vBWGG0 = BWGG0.*avdir;%(veinCmap_GG0>0);
%         aBWGG0 = BWGG0.*(1-avdir);%(arteriesAmap_GG0>0);
%         vBWGGV0 = BWGGV0.*avdir;%(veinAmap_GGV0>0);
%         aBWGGV0 = BWGGV0.*(1-avdir);%(arteriesAmap_GGV0>0);
%         vBWGGFV0 = BWGGFV0.*avdir;%(veinAmap_GGFV0>0);
%         aBWGGFV0 = BWGGFV0.*(1-avdir);%(arteriesAmap_GGFV0>0);
%         vBWPDI0 = BWPDI0.*avdir;%(veinAmap_PDI0>0);
%         aBWPDI0 = BWPDI0.*(1-avdir);%(arteriesAmap_PDI0>0);
%         vBWGGF0 = BWGGF0.*avdir;%(veinAmap_GGF0>0);
%         aBWGGF0 = BWGGF0.*(1-avdir);%(arteriesAmap_GGF0>0);

vroiGG = squeeze(sum(sum(sumGG.*repmat(vBWGG0.*vmsk,[1, 1,nRpt]),1),2));
aroiGG = squeeze(sum(sum(sumGG.*repmat(aBWGG0.*vmsk,[1, 1,nRpt]),1),2));
% Fig = figure;set(Fig,'Position',[400 400 800 170])                     
% yyaxis left;plot(-roiGG, 'b', 'LineWidth', 1);ylabel('G1-fUS(CBF)'); xlabel('Time[s]');                         
% yyaxis right; plot(stim);grid on;
vroiGGV = squeeze(sum(sum(sumGGV.*repmat(vBWGGV0.*vmsk,[1, 1,nRpt]),1),2));
aroiGGV = squeeze(sum(sum(sumGGV.*repmat(aBWGGV0.*vmsk,[1, 1,nRpt]),1),2));
% Fig = figure;set(Fig,'Position',[400 400 800 170]) 
% yyaxis left;plot(-roiGGV, 'b', 'LineWidth', 1);ylabel('G1-fUS(CBV)'); xlabel('Time[s]'); 
% yyaxis right; plot(stim);grid on;
vroiGGFV = squeeze(sum(sum(sumGGFV.*repmat(vBWGGFV0.*vmsk,[1, 1,nRpt]),1),2));
aroiGGFV = squeeze(sum(sum(sumGGFV.*repmat(aBWGGFV0.*vmsk,[1, 1,nRpt]),1),2));
% Fig = figure;set(Fig,'Position',[400 400 800 170]) 
% yyaxis left;plot(-roiGGFV, 'b', 'LineWidth', 1);ylabel('G1-fUS(novel)'); xlabel('Time[s]'); 
% yyaxis right; plot(stim);grid on;
vroiPDI = squeeze(sum(sum(squeeze(PDI(:,:,3,:)).*repmat(vBWPDI0.*vmsk,[1, 1,nRpt]),1),2));
aroiPDI = squeeze(sum(sum(squeeze(PDI(:,:,3,:)).*repmat(aBWPDI0.*vmsk,[1, 1,nRpt]),1),2));
% yyaxis left;plot(roiPDI, 'b', 'LineWidth', 1);ylabel('PDI'); xlabel('Time[s]'); 
% yyaxis right; plot(stim);grid on;

vroiGGF = squeeze(sum(sum(sumGGF.*repmat(vBWGGF0.*vmsk,[1, 1,nRpt]),1),2));
aroiGGF = squeeze(sum(sum(sumGGF.*repmat(aBWGGF0.*vmsk,[1, 1,nRpt]),1),2));
% Fig = figure;set(Fig,'Position',[400 400 800 170]) 
% yyaxis left;plot(-roiGGFV, 'b', 'LineWidth', 1);ylabel('G1-fUS(novel)'); xlabel('Time[s]'); 
% yyaxis right; plot(stim);grid on;

%% arteries/veins ROI trial averaged time course GG vs. PDI

ROIGG.vein = averagedTrials(vroiGG, trial);
ROIGG.artery = averagedTrials(aroiGG, trial);

ROIGGV.vein = averagedTrials(vroiGGV, trial);
ROIGGV.artery = averagedTrials(aroiGGV, trial);

ROIGGFV.vein = averagedTrials(vroiGGFV, trial);
ROIGGFV.artery = averagedTrials(aroiGGFV, trial);

ROIPDI.vein = averagedTrials(vroiPDI, trial);
ROIPDI.artery = averagedTrials(aroiPDI, trial);

ROIGGF.vein = averagedTrials(vroiGGF, trial);
ROIGGF.artery = averagedTrials(aroiGGF, trial);

%% GGFV contrast compared to GG GGV PDI

[snrGG, snrGGnum, snrGGden, snrGG_final] = calSNRImage(mtrialsumGG_bar, trial);
[snrGGV, snrGGVnum, snrGGVden, snrGGV_final] = calSNRImage(mtrialsumGGV_bar, trial);
[snrGGFV, snrGGFVnum, snrGGFVden, snrGGFV_final] = calSNRImage(mtrialsumGGFV_bar, trial);
[snrPDI, snrPDInum, snrPDIden, snrPDI_final] = calSNRImage(mtrialPDI_bar, trial);
[snrGGF, snrGGFnum, snrGGFden, snrGGF_final] = calSNRImage(mtrialsumGGF_bar, trial);
% figure; 
% subplot(221);
% imagesc(snrGG);colorbar; title('G1(rCBFv)')
% subplot(222)
% imagesc(snrGGV);colorbar;title('G1(rCBV)')
% subplot(223)
% imagesc(snrGGFV);colorbar;title('G1(optimal)')
% subplot(224) 
% imagesc(snrPDI);colorbar;title('PDI')

%%
Fig = figure;
set(Fig, 'Position', [500 500 900 450])
subplot(231);
imagesc(snrGG_final); title('G1(rCBFv)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(232)
imagesc(snrGGV_final);title('G1(rCBV)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(233)
imagesc(snrGGF_final);title('G1(rCBF)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(234)
imagesc(snrGGFV_final);title('G1(optimal)'); axis image; caxis([0.5,2]);colorbar;axis off;
subplot(235)
imagesc(snrPDI_final);title('PDI'); axis image; caxis([0.5,2]);colorbar;axis off;

%%
Fig = figure;
set(Fig, 'Position', [500 500 900 450])
subplot(231);
imagesc(snrGG_final); title('G1(rCBFv)'); axis image; 
subplot(232)
imagesc(snrGGV_final);title('G1(rCBV)'); axis image; 
subplot(233)
imagesc(snrGGF_final);title('G1(rCBF)'); axis image;
subplot(234)
imagesc(snrGGFV_final);title('G1(optimal)'); axis image;
subplot(235)
imagesc(snrPDI_final);title('PDI'); axis image;


% %%
% cmax = 100;
% figure; 
% subplot(221);
% imagesc(snrGGnum);colorbar;title('G1(rCBFv)'); %caxis([0,0.2]);
% subplot(222)
% imagesc(snrGGVnum);colorbar;title('G1(rCBV)'); %caxis([0,0.2]);
% subplot(223)
% imagesc(abs(snrGGFVnum));colorbar;title('G1(optimal)');% caxis([0,0.2]);
% subplot(224) 
% imagesc(snrPDInum);colorbar;title('PDI'); %caxis([0,0.2]);

% %%
% cmax = 100;
% figure; 
% subplot(221);
% imagesc(snrGGden);colorbar;title('G1(rCBFv)')
% subplot(222)
% imagesc(snrGGVden);colorbar;title('G1(rCBV)')
% subplot(223)
% imagesc(snrGGFVden);colorbar;title('G1(optimal)')
% subplot(224) 
% imagesc(snrPDIden);colorbar;title('PDI')
% %%
% Fig = figure;
% set(Fig, 'Position', [500 500 900 450])
% subplot(231);
% imagesc(snrGGdB);caxis([0,40]); title('G1(rCBFv)'); axis image
% subplot(232)
% imagesc(snrGGVdB);caxis([0,40]);title('G1(rCBV)'); axis image
% subplot(233)
% imagesc(snrGGFdB);caxis([0,40]);title('G1(rCBF)'); axis image
% subplot(234)
% imagesc(snrGGFVdB);caxis([0,40]);title('G1(optimal)'); axis image
% subplot(235) 
% imagesc(snrPDIdB);caxis([0,40]);title('PDI'); axis image
% %%
% slt = [18,154];%[53,131];%[95,91];
% figure;
% subplot(221);
% plot(squeeze(mtrialsumGG_bar(slt(1),slt(2),:)));title(['(',num2str(slt),') ', num2str(snrGG(slt(1),slt(2)))]);
% subplot(222);
% plot(squeeze(mtrialsumGGV_bar(slt(1),slt(2),:)));title([num2str(snrGGV(slt(1),slt(2)))]);
% subplot(223);
% plot(squeeze(mtrialsumGGFV_bar(slt(1),slt(2),:)));title([num2str(snrGGFV(slt(1),slt(2)))]);
% subplot(224);
% plot(squeeze(mtrialPDI_bar(slt(1),slt(2),:)));title([num2str(snrPDI(slt(1),slt(2)))]);

%% amplitude of time course image GG GGV GGF

window = [trial.nRest+1, trial.nRest+trial.nStim+1];

ampGG0 = max(mtrialsumGG_ratio(:,:,window(1):window(2)),[],3);
ampGGV0 = max(mtrialsumGGV_ratio(:,:,window(1):window(2)),[],3);
ampGGF0 = max(mtrialsumGGF_ratio(:,:,window(1):window(2)),[],3);
% ampGGFV = max(mtrialsumGGFV_ratio(:,:,window(1):window(2)),[],3);
% ampPDI = max(mtrialPDI_ratio(:,:,window(1):window(2)),[],3);

% ampGG0 = convn(ampGG0,BB,'same');
% ampGGV0 = convn(ampGGV0,BB,'same');
% ampGGF0 = convn(ampGGF0,BB,'same');

ampGG = ampGG0.*(actmap_GG0>0);
ampGGV = ampGGV0.*(actmap_GGV0>0);
ampGGF = ampGGF0.*(actmap_GG0>0);

R  = (ampGGF0-1)./(ampGGV0-1);
actR = R.*(actmap_GGF0>0);
% figure; imagesc(R.*(actmap_GGF0>0));axis image; caxis([1.5 3]);

%% arteries/vein amplitude of time course image GG GGV GGF

window = [trial.nRest+1, trial.nRest+trial.nStim+1];

ampGG0_v = ampGG0.*(veinCmap_GG0>0);
ampGG0_a = ampGG0.*(arteriesCmap_GG0>0);

ampGGV0_v = ampGGV0.*(veinCmap_GG0>0);
ampGGV0_a = ampGGV0.*(arteriesCmap_GG0>0);

ampGGF0_v = ampGGF0.*(veinCmap_GG0>0);
ampGGF0_a = ampGGF0.*(arteriesCmap_GG0>0);

ampGG_v = ampGG0_v.*(actmap_GG0>0).*MSKGG;
ampGG_a = ampGG0_a.*(actmap_GG0>0).*MSKGG;
ampGGV_v = ampGGV0_v.*(actmap_GGV0>0).*MSKGGV;
ampGGV_a = ampGGV0_a.*(actmap_GGV0>0).*MSKGGV;
ampGGF_v = ampGGF0_v.*(actmap_GGF0>0).*MSKGGF;
ampGGF_a = ampGGF0_a.*(actmap_GGF0>0).*MSKGGF;

% R = (ampGGF-1)./(ampGGV-1);
% figure; imagesc(ampGG_v);axis image; caxis([1 1.3]);
% figure; imagesc(ampGG_a);axis image; caxis([1 1.3]);




%% save analzed data
save([FilePath, 'ROIdata.mat'], 'ROIGG', 'ROIGGV', 'ROIGGF', 'ROIGGFV', 'ROIPDI');
%% save manually selected roi
save([FilePath, 'BW.mat'], 'BWGG0', 'BWGGV0', 'BWGGFV0', 'BWPDI0', 'aBWGG0', 'vBWGG0', 'aBWGGV0', 'vBWGGV0');