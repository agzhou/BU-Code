% JCF variation Testing on bio data
%% Initialize data
clearvars
close all

% Extract Current Path
currentDir = matlab.desktop.editor.getActiveFilename; 
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

dataFile{1} = dataFilePath + "Boas Lab Data\Calf45AngMultiPower.mat";
filetype = 14;

% 9LD Calf Data multipower
% dataFile{1} = dataFilePath + "Calf Data\Calf45AngleMultiPowerGE9LD.mat";
% filetype = 10;

%% Load Data
[p,RFData] = initParams(dataFile,filetype);
nFrames = size(RFData,3);

%% Define grid
pL = computeNewGrid(p,[1,p.szX],[1,p.szZ],p.szX*4,p.szZ*4);

% pS = computeNewGrid(pL,[1,400],[100,600]);
pS = computeNewGrid(pL,[pL.szX/2,pL.szX],[100,800]);

%% Generate K matrix
cRF = computeCRF(RFData(:,:,4),pS);
beamform = reconraw.DASBModeOffline(pS);

tic; idxtMTX = beamform.computeDASFullKSpace(cRF); toc
tic; idxtMTX = reshape(idxtMTX,[pS.numEl,pS.na,pS.nPoints]); toc

%% Compute JCF variations

% No exponential factor
tic; JCF1 = JCFtestAlpha(idxtMTX,pS,1); toc
% Squared
tic; JCF2 = JCFtestAlpha(idxtMTX,pS,2); toc
% Cubed
tic; JCF3 = JCFtestAlpha(idxtMTX,pS,3); toc
% 4th power
tic; JCF4 = JCFtestAlpha(idxtMTX,pS,4); toc

% Repeat all of the above but with no conj factor in numerator
tic; JCF5 = JCFtestAlpha2(idxtMTX,pS,1); toc
tic; JCF6 = JCFtestAlpha2(idxtMTX,pS,2); toc
tic; JCF7 = JCFtestAlpha2(idxtMTX,pS,3); toc
tic; JCF8 = JCFtestAlpha2(idxtMTX,pS,4); toc

% Current version of JCF
tic; JCF9 = JCFCPP(idxtMTX,pS); toc

%% Plot results
figure
tiledlayout(2,4)
nexttile(1,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF1,0.1)
axis image
title('Raised to 1')
nexttile(2,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF2,0.1)
axis image
title('Raised to 2')
nexttile(3,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF3,0.1)
axis image
title('Raised to 3')
nexttile(4,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF4,0.1)
axis image
title('Raised to 4')

figure
tiledlayout(2,4)
nexttile(1,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF5,0.1)
axis image
title('No Conj Raised to 1')
nexttile(2,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF6,0.1)
axis image
title('No Conj Raised to 2')
nexttile(3,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF7,0.1)
axis image
title('No Conj Raised to 3')
nexttile(4,[2,1])
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF8,0.1)
axis image
title('No Conj Raised to 4')
figure
plotGammaScaleImage(pS.xCoord*1e3,pS.zCoord*1e3,JCF9,0.1)
axis image
title('Default CPP')

%% Helper Functions
function [weightedBFM] = JCFtestAlpha(idxtMTX,p,alpha)

weightedBFM = zeros(size(idxtMTX));
denomFactor = double(p.numEl*p.na);

for i = 1:p.nPoints
    kMTX = idxtMTX(:,:,i);
    kMTX2 = abs(kMTX).^alpha;
    
    weightMTX = abs(sum(kMTX,1).*conj(sum(kMTX,2))).^alpha ./ ...
      (denomFactor.^(alpha-1)*sum(kMTX2,1).*sum(kMTX2,2));
      

    weightMTX(isnan(weightMTX)) = 0;
    weightedBFM(:,:,i) = weightMTX;
end

weightedBFM = idxtMTX.*(weightedBFM);
weightedBFM = reshape(squeeze(sum(weightedBFM,[1,2])/denomFactor),[p.szZ,p.szX]);

end

function [weightedBFM] = JCFtestAlpha2(idxtMTX,p,alpha)

weightedBFM = zeros(size(idxtMTX));
denomFactor = double(p.numEl*p.na);

for i = 1:p.nPoints
    kMTX = idxtMTX(:,:,i);
    kMTX2 = abs(kMTX).^alpha;
    
    weightMTX = abs(sum(kMTX,1).*(sum(kMTX,2))).^alpha ./ ...
      (denomFactor.^(alpha-1)*sum(kMTX2,1).*sum(kMTX2,2));
      

    weightMTX(isnan(weightMTX)) = 0;
    weightedBFM(:,:,i) = weightMTX;
end

weightedBFM = idxtMTX.*(weightedBFM);
weightedBFM = reshape(squeeze(sum(weightedBFM,[1,2])/denomFactor),[p.szZ,p.szX]);

end

