function [ROIGG, ROIGGV, ROIGGF, ROIGGFV, ROIPDI, Ratio, Peak, T2Peak, CNR, anpratio] = g1fUS_analysis_main
addpath D:\g1_based_fUS\A-US-fUS-g1\SubFunctions\
newSlect = 1;
[FileName,FilePath]=uigetfile('V:\G1based fUS data');
load([FilePath,'GG_invivo_results.mat']);
dir1 = load([FilePath,'GG_invivo_results.mat'],'dir');
dir0 = dir1.dir;

[~,~,nRpt] = size(sumGG);
tCoor = [1:trial.nlength]-trial.nRest;

bmsk = ones(nz,nx);
bmsk(1:5,:) = 0;
bmsk(120:nz,:) = 0; % rough brain mask

% %% Global HRF 
% %hrf = hemodynamicResponse(1,[2 16 0.5 1 20 0]);
% %  hrf = hemodynamicResponse(1,[1.5 10 0.5 1 20 0 16]);
% hrf = hemodynamicResponse(0.1,[2 6 .5 1 10 0 12]);
% 
% Stim = zeros(trial.nlength*10,1);
% Stim(5*10+1:5*10+trial.nStim*10,:)=1;
% stimhrf0 = filter(hrf,1,Stim);
% stimhrf = interp1(1:trial.nlength*10, stimhrf0, 1:10:trial.nlength*10, 'linear')';
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

[actmap_GG0, coefmap_GG0] = CoorCoeffMap(mtrialsumGG.*vmsk.*bmsk, stimhrf', 0);

[actmap_GGV0, coefmap_GGV0] = CoorCoeffMap(mtrialsumGGV_bar.*vmsk.*bmsk, stimhrf', 0);

[actmap_GGFV0, coefmap_GGFV0] = CoorCoeffMap(mtrialsumGGFV_bar.*vmsk.*bmsk, stimhrf', 0);

[actmap_PDI0, coefmap_PDI0] = CoorCoeffMap(mtrialPDI, stimhrf', 0);

[actmap_GGF0, coefmap_GGF0] = CoorCoeffMap(mtrialsumGGF_bar.*vmsk.*bmsk, stimhrf', 0);

%% arteries/veins divide

avdir = zeros(nz,nx);
avdir(dir0>0) = 1;
% figure; 
% subplot(121); imagesc(avdir); axis image; colormap('gray'); colorbar; title('vein g1 mask');
% subplot(122); imagesc((1-avdir)); axis image; colormap('gray'); colorbar; title('artery g1 mask');

%% Directional Doppler based a/v mask
vPDm = (squeeze(mean(PDI(:,:,1,:)./(PDI(:,:,1,:)+PDI(:,:,2,:)),4)));
aPDm = (squeeze(mean(PDI(:,:,2,:)./(PDI(:,:,1,:)+PDI(:,:,2,:)),4)));
% 
% vPDm = log(squeeze(mean(PDI(:,:,1,:),4)));
% aPDm = log(squeeze(mean(PDI(:,:,2,:),4)));
% figure; 
% subplot(121); imagesc(vPDm); axis image; colormap('gray'); colorbar; title('vein Directional PD');
% subplot(122); imagesc(aPDm); axis image; colormap('gray'); colorbar; title('artery Directional PD');
vPDmsk = vPDm>1.01*mean(vPDm(:));%(0.99*max(vPDm(:))+min(vPDm(:)))/2;
aPDmsk = aPDm>1.01*mean(aPDm(:));%(0.9*max(aPDm(:))+min(aPDm(:)))/2;

%% 
vMSK = vPDmsk;%.*avdir;
aMSK = aPDmsk;%.*(1-avdir);
figure; 
subplot(121); imagesc(vMSK); axis image; colormap('gray'); colorbar; title('vein mask');
subplot(122); imagesc(aMSK); axis image; colormap('gray'); colorbar; title('artery mask');

%% count a/v pixels in actmap
anp_GG0 = sum((actmap_GG0>0).*MSKGG.*aMSK,'all'); 
vnp_GG0 = sum((actmap_GG0>0).*MSKGG.*vMSK,'all');
anpratio.GG = anp_GG0./(anp_GG0+vnp_GG0);

anp_GGV0 = sum((actmap_GGV0>0).*MSKGGV.*aMSK,'all');
vnp_GGV0 = sum((actmap_GGV0>0).*MSKGGV.*vMSK,'all');
anpratio.GGV = anp_GGV0./(anp_GGV0+vnp_GGV0);

anp_GGF0 = sum((actmap_GGF0>0).*MSKGGF.*aMSK,'all');
vnp_GGF0 = sum((actmap_GGF0>0).*MSKGGF.*vMSK,'all');
anpratio.GGF = anp_GGF0./(anp_GGF0+vnp_GGF0);

anp_GGFV0 = sum((actmap_GGFV0>0).*MSKGGF.*aMSK,'all');
vnp_GGFV0 = sum((actmap_GGFV0>0).*MSKGGF.*vMSK,'all');
anpratio.GGFV = anp_GGFV0./(anp_GGFV0+vnp_GGFV0);

anp_PDI0 = sum((actmap_PDI0>0).*MSKGGF.*aMSK,'all');
vnp_PDI0 = sum((actmap_PDI0>0).*MSKGGF.*vMSK,'all');
anpratio.PDI = anp_PDI0./(anp_PDI0+vnp_PDI0);

%%
veinCmap_GG0 = coefmap_GG0.*vMSK;
arteriesCmap_GG0 = coefmap_GG0.*aMSK;
veinAmap_GG0 = actmap_GG0.*vMSK.*bmsk;
arteriesAmap_GG0 = actmap_GG0.*aMSK.*bmsk;

veinCmap_GGV0 = coefmap_GGV0.*vMSK;
arteriesCmap_GGV0 = coefmap_GGV0.*aMSK;
veinAmap_GGV0 = actmap_GGV0.*vMSK.*bmsk;
arteriesAmap_GGV0 = actmap_GGV0.*aMSK.*bmsk;

veinCmap_GGFV0 = coefmap_GGFV0.*vMSK;
arteriesCmap_GGFV0 = coefmap_GGFV0.*aMSK;
veinAmap_GGFV0 = actmap_GGFV0.*vMSK.*bmsk;
arteriesAmap_GGFV0 = actmap_GGFV0.*aMSK.*bmsk;

veinCmap_PDI0 = coefmap_PDI0.*vMSK;
arteriesCmap_PDI0 = coefmap_PDI0.*aMSK;
veinAmap_PDI0 = actmap_PDI0.*vMSK.*bmsk;
arteriesAmap_PDI0 = actmap_PDI0.*aMSK.*bmsk;

veinCmap_GGF0 = coefmap_GGF0.*vMSK;
arteriesCmap_GGF0 = coefmap_GGF0.*aMSK;
veinAmap_GGF0 = actmap_GGF0.*vMSK.*bmsk;
arteriesAmap_GGF0 = actmap_GGF0.*aMSK.*bmsk;
%% GGFV contrast compared to GG GGV PDI

[snrGG, snrGGnum, snrGGden, snrGG_final] = calSNRImage(mtrialsumGG_bar, trial);
[snrGGV, snrGGVnum, snrGGVden, snrGGV_final] = calSNRImage(mtrialsumGGV_bar, trial);
[snrGGFV, snrGGFVnum, snrGGFVden, snrGGFV_final] = calSNRImage(mtrialsumGGFV_bar, trial);
[snrPDI, snrPDInum, snrPDIden, snrPDI_final] = calSNRImage(mtrialPDI_bar, trial);
[snrGGF, snrGGFnum, snrGGFden, snrGGF_final] = calSNRImage(mtrialsumGGF_bar, trial);

CNR.GG = snrGG_final;
CNR.GGV = snrGGV_final;
CNR.GGFV = snrGGFV_final;
CNR.PDI = snrPDI_final;
CNR.GGF = snrGGF_final;
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

%% amplitude of time course image GG GGV GGF

window = [trial.nRest+1, trial.nRest+trial.nStim+1];

ampGG0 = max(mtrialsumGG_ratio(:,:,window(1):window(2)),[],3);
ampGGV0 = max(mtrialsumGGV_ratio(:,:,window(1):window(2)),[],3);
ampGGF0 = max(mtrialsumGGF_ratio(:,:,window(1):window(2)),[],3);

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

ampGG_v = ampGG0_v.*(actmap_GG0>0);
ampGG_a = ampGG0_a.*(actmap_GG0>0);
ampGGV_v = ampGGV0_v.*(actmap_GGV0>0);
ampGGV_a = ampGGV0_a.*(actmap_GGV0>0);
ampGGF_v = ampGGF0_v.*(actmap_GGF0>0);
ampGGF_a = ampGGF0_a.*(actmap_GGF0>0);

%% ROIGG ROIGGV defined 
 
if exist('MSKGG','var')==0
    MSKGG = selectROI(coefmap_GG0);
    MSKGGV = selectROI(coefmap_GGV0);
    MSKGGF = selectROI(coefmap_GGF0);          
end

%% ROI averaged time course

lib = {'rCBFv_based', 'rCBV_based','rCBF_based'} ;
for i = 1:numel(lib)
    ROItype = lib{i};

newSelect = 1;
%ROItype = 'rCBFv_based';%'userdefined';%'actmap';
differentROI = 0;

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
        BWGGFV0 = (actmap_GGFV0>0).*MSKGGF; % GGF
        BWPDI0 = (actmap_PDI0>0).*MSKGGF; % GGF
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

roiGGV = squeeze(sum(sum(sumGGV.*repmat(BWGGF0.*vmsk,[1, 1,nRpt]),1),2));

roiGGFV = squeeze(sum(sum(sumGGFV.*repmat(BWGGFV0.*vmsk,[1, 1,nRpt]),1),2));

roiPDI = squeeze(sum(sum(squeeze(PDI(:,:,3,:)).*repmat(BWPDI0,[1, 1,nRpt]),1),2));

roiGGF = squeeze(sum(sum(sumGGF.*repmat(BWGGF0.*vmsk,[1, 1,nRpt]),1),2));

%% trial averaged time course GG vs. PDI

ROIGG{i} = averagedTrials(roiGG, trial);

ROIGGV{i} = averagedTrials(roiGGV, trial);

ROIGGFV{i} = averagedTrials(roiGGFV, trial);

ROIPDI{i} = averagedTrials(roiPDI, trial);

ROIGGF{i} = averagedTrials(roiGGF, trial);


%% arteries/veins ROI averaged time course

newSelect = 1;
% ROItype = 'rCBV_based';%'rCBV_based';%'userdefined';%'actmap';%;%
differentROI = 0;

if differentROI == 1
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
            vROI = selectROI(veinCmap_GGV0);
            aROI = selectROI(arteriesCmap_GGV0);

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

vroiGG = squeeze(sum(sum(sumGG.*repmat(vBWGG0.*vmsk,[1, 1,nRpt]),1),2));
aroiGG = squeeze(sum(sum(sumGG.*repmat(aBWGG0.*vmsk,[1, 1,nRpt]),1),2));

vroiGGV = squeeze(sum(sum(sumGGV.*repmat(vBWGGV0.*vmsk,[1, 1,nRpt]),1),2));
aroiGGV = squeeze(sum(sum(sumGGV.*repmat(aBWGGV0.*vmsk,[1, 1,nRpt]),1),2));

vroiGGFV = squeeze(sum(sum(sumGGFV.*repmat(vBWGGFV0.*vmsk,[1, 1,nRpt]),1),2));
aroiGGFV = squeeze(sum(sum(sumGGFV.*repmat(aBWGGFV0.*vmsk,[1, 1,nRpt]),1),2));

vroiPDI = squeeze(sum(sum(squeeze(PDI(:,:,3,:)).*repmat(vBWPDI0.*vmsk,[1, 1,nRpt]),1),2));
aroiPDI = squeeze(sum(sum(squeeze(PDI(:,:,3,:)).*repmat(aBWPDI0.*vmsk,[1, 1,nRpt]),1),2));

vroiGGF = squeeze(sum(sum(sumGGF.*repmat(vBWGGF0.*vmsk,[1, 1,nRpt]),1),2));
aroiGGF = squeeze(sum(sum(sumGGF.*repmat(aBWGGF0.*vmsk,[1, 1,nRpt]),1),2));

%% arteries/veins ROI trial averaged time course GG vs. PDI

ROIGG{i}.vein = averagedTrials(vroiGG, trial);
ROIGG{i}.artery = averagedTrials(aroiGG, trial);

ROIGGV{i}.vein = averagedTrials(vroiGGV, trial);
ROIGGV{i}.artery = averagedTrials(aroiGGV, trial);

ROIGGFV{i}.vein = averagedTrials(vroiGGFV, trial);
ROIGGFV{i}.artery = averagedTrials(aroiGGFV, trial);

ROIPDI{i}.vein = averagedTrials(vroiPDI, trial);
ROIPDI{i}.artery = averagedTrials(aroiPDI, trial);

ROIGGF{i}.vein = averagedTrials(vroiGGF, trial);
ROIGGF{i}.artery = averagedTrials(aroiGGF, trial);

%% measure peakratio peak time2peak 
Ratio{i}.all = (max(movmean(ROIGGF{i}.m, 3)/100)-1)/(max(movmean(ROIGGV{i}.m, 3)/100)-1);
Ratio{i}.vein = (max(movmean(ROIGGF{i}.vein.m, 3)/100)-1)/(max(movmean(ROIGGV{i}.vein.m, 3)/100)-1); 
Ratio{i}.artery = (max(movmean(ROIGGF{i}.artery.m, 3)/100)-1)/(max(movmean(ROIGGV{i}.artery.m, 3)/100)-1); 

[vpGG,~] = max(movmean(ROIGG{i}.vein.m(window(1):window(2)), 3)/100);
[vpGGV,~] = max(movmean(ROIGGV{i}.vein.m(window(1):window(2)), 3)/100);
[vpGGF,~] = max(movmean(ROIGGF{i}.vein.m(window(1):window(2)), 3)/100);

[apGG,~] = max(movmean(ROIGG{i}.artery.m(window(1):window(2)), 3)/100);
[apGGV,~] = max(movmean(ROIGGV{i}.artery.m(window(1):window(2)), 3)/100);
[apGGF,~] = max(movmean(ROIGGF{i}.artery.m(window(1):window(2)), 3)/100);

interpx = interp1(1:trial.nlength, ROIGG{i}.vein.m, 1:0.2:trial.nlength, "spline");
[~,vtpGG0] = max(interpx(window(1)*5:window(2)*5));
vtpGG = vtpGG0/5;
interpx = interp1(1:trial.nlength, ROIGGV{i}.vein.m, 1:0.2:trial.nlength, "spline");
[~,vtpGGV0] = max(interpx(window(1)*5:window(2)*5));
vtpGGV = vtpGGV0/5;
interpx = interp1(1:trial.nlength, ROIGGF{i}.vein.m, 1:0.2:trial.nlength, "spline");
[~,vtpGGF0] = max(interpx(window(1)*5:window(2)*5));
vtpGGF = vtpGGF0/5;

interpx = interp1(1:trial.nlength, ROIGG{i}.artery.m, 1:0.2:trial.nlength, "spline");
[~,atpGG0] = max(interpx(window(1)*5:window(2)*5));
atpGG = atpGG0/5;
interpx = interp1(1:trial.nlength, ROIGGV{i}.artery.m, 1:0.2:trial.nlength, "spline");
[~,atpGGV0] = max(interpx(window(1)*5:window(2)*5));
atpGGV = atpGGV0/5;
interpx = interp1(1:trial.nlength, ROIGGF{i}.artery.m, 1:0.2:trial.nlength, "spline");
[~,atpGGF0] = max(interpx(window(1)*5:window(2)*5));
atpGGF = atpGGF0/5;

Peak{i}.GG.vein = vpGG;
Peak{i}.GG.artery = apGG;
Peak{i}.GGV.vein = vpGGV;
Peak{i}.GGV.artery = apGGV;
Peak{i}.GGF.vein = vpGGF;
Peak{i}.GGF.artery = apGGF;

T2Peak{i}.GG.vein = vtpGG;
T2Peak{i}.GG.artery = atpGG;
T2Peak{i}.GGV.vein = vtpGGV;
T2Peak{i}.GGV.artery = atpGGV;
T2Peak{i}.GGF.vein = vtpGGF;
T2Peak{i}.GGF.artery = atpGGF;
end
save([FilePath, 'GG_invivo_results.mat'],'sumGG', 'sumGGV', 'sumGGFV', 'PDI', 'eqNoise', 'trial', 'dir','dir0', 'vmsk', 'nx', 'nz', 'ntau', 'MSKGG', 'MSKGGV', 'MSKGGF');
save([FilePath, 'ROIdata.mat'],'ROIGG', 'ROIGGV', 'ROIGGF', 'ROIGGFV', 'ROIPDI');
end



