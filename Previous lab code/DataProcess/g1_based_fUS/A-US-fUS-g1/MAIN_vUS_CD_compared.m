addpath D:\g1_based_fUS\A-US-fUS-g1\SubFunctions\

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
    sumGG0 = myFile.sumGG0;
    sumGGV0 = myFile.sumGGV0;
    sumGGFV0 = myFile.sumGGFV0;

  PDIHP = myFile.PDIHHP;
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
ImgBG = mean(squeeze(PDI(:,:,3,:)),3).^0.35;
trialPDI = reshape(squeeze(PDI(:,:,3,trial.nBase-5+1: end-5)), [nz,nx,trial.nlength,trial.n]);
mtrialPDI = median(trialPDI,4);
[nz,nx,~,~] = size(VV0);

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
%xCoor = (1:ntau)*1e3/P.CCFR; % tau axis[ms]

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
 figure; subplot(121);imagesc(-squeeze(mean(VV0(:,:,1,:),4)));axis image; colorbar;title('ascending veins');
subplot(122);imagesc(squeeze(mean(VV0(:,:,2,:),4)));axis image;colorbar;title('decending arteries');
%%

trialV0 = reshape(sumGGV(:,:,1:2,trial.nBase-5+1: end-5), [nz,nx,2,trial.nlength,trial.n]);
mtrialV0 = median(trialV0,5);
mtrialV0_bar = normalize(mtrialV0, 4, 'norm');
mtrialV0_ratio = mtrialV0_bar./mean(mtrialV0_bar(:,:,:,1: trial.nRest),4);
coefthld = 1.65; % 1.65 p<0.05; 2.35 p<0.01

[actmap_VV0(:,:,1), coefmap_VV0(:,:,1)] = CoorCoeffMap(squeeze(mtrialV0(:,:,1,:)), stimhrf', 1, coefthld);
[actmap_VV0(:,:,2), coefmap_VV0(:,:,2)] = CoorCoeffMap(squeeze(mtrialV0(:,:,2,:)), stimhrf', 1, coefthld);

% BWVV0 = logical((actmap_VV0(:,:,1)>0)+(actmap_VV0(:,:,2)>0)).*MSKGGF;
% roiV0 = squeeze(sum(sum(VVcz0.*repmat(BWVV0.*vmsk,[1, 1,1,nRpt]),1),2))/sum(sum(BWVV0.*vmsk));

%%
vBWVV0 = (actmap_VV0(:,:,1)>0).*MSKGG;
aBWVV0 = (actmap_VV0(:,:,2)>0).*MSKGG;
figure; subplot(121); imagesc(vBWVV0); axis image;title('veinMsk'); subplot(122); imagesc(aBWVV0); axis image;title('arteryMsk');

roiV0(1,:) = squeeze(sum(sum(squeeze(sumGGV(:,:,1,:)).*repmat(vBWVV0.*vmsk,[1, 1,nRpt]),1),2))/sum(vBWVV0.*vmsk,'all');
roiV0(2,:) = squeeze(sum(sum(squeeze(sumGGV(:,:,2,:)).*repmat(aBWVV0.*vmsk,[1, 1,nRpt]),1),2))/sum(aBWVV0.*vmsk,'all');

ROIVV.vein = averagedTrials(roiV0(1,:), trial);
ROIVV.artery = averagedTrials(roiV0(2,:), trial);
        
window = [trial.nRest+1, trial.nRest+trial.nStim+1];            
ampVV0 = max(mtrialV0_ratio(:,:,:,window(1):window(2)),[],4);
ampVV = ampVV0.*(actmap_VV0>0);

stim = trial.stim;
Fig = figure;set(Fig,'Position',[400 400 700 700])
subplot(411)
yyaxis left;plot(roiV0(1,:), 'b', 'LineWidth', 1);ylabel('g1-vUS(CBFv) vein'); xlabel('Time[s]'); 
yyaxis right; plot(stim);grid on;
subplot(412)
yyaxis left;plot(roiV0(2,:), 'b', 'LineWidth', 1);ylabel('g1-vUS(CBFv) artery'); xlabel('Time[s]'); 
yyaxis right; plot(stim);grid on;

        arteriesAmap = actmap_VV0(:,:,2);
        veinAmap = actmap_VV0(:,:,1);
        ROI = ROIVV;

        amp0_a = ampVV(:,:,2);
        amp0_v = ampVV(:,:,1);