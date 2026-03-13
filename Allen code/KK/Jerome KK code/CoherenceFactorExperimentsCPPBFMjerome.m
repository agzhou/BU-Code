% Coherence Factor Experiments
clearvars
close all

%% Select Data

dataFilePath = "G:\Shared drives\Biomicroscopy Lab\Active Projects\Ultrasound\Datasets\";

% CPP Test Data Options

% CPP Test Data Spherical target
% dataFile{1} = dataFilePath + "CPP Test Data\Type4Sphere_20mm.mat";
% filetype = 2;

%  CPP Test Data Tall Calibration point targets
% dataFile{1} = dataFilePath + "CPP Test Data\TallVerticalPoints.mat";
% filetype = 2;

% CPP Test Data Tall Elastic Targets + lateral resolution targets
% dataFile{1} = dataFilePath + "CPP Test Data\TallElasticTargets.mat";
% filetype = 2;


% Human Data Options
% dataFile{1} = dataFilePath + "Calf Data\RFdataAngle1unStrainedSingleFile.mat";
% filetype = 2;

%Boas lab data
dataFile{1} = dataFilePath + "Boas Lab Data\Calf15AngMultiPower.mat";
filetype = 14;

%% Load Data
[p,RFData] = initParams(dataFile,filetype);
p.tShift = 0.0;
p.txPL = zeros(length(p.TXangle));

% Only if not relying on default grid

% Description: There are 5 parameters within the main parameter structure
% that need to be updated if we want a new beamforming grid to operate on:
%       xCoord -    x Coordinates to be beamformed at (units in mm)
%       zCoord -    z Coordinates to be beamformed at (units in mm)
%       szX -       length(xCoord)
%       szZ -       length(zCoord)-1. Note: Do not forget the -1 here!!!
%       nPoints -   total number of beamforming points = szX*szZ;
%
% The xCoord and zCoord can be set at any points that the user wants. Only
% restrictions are that zCoord must be positive and the szX,szZ and nPoints
% values must be consistent with the definitions given.

% Here we have an example where we take the default grid and recompute a
% finer spacing in the x and z directions over a smaller region. Then we
% initialize the beamformer class and proceed as normal. This example is
% set up to zoom in on the resolution targets in the Tall Blue GSE phantom.

% p.xCoord = interp1(((p.szX/4+1):(3*p.szX/4)),p.xCoord((p.szX/4+1):(3*p.szX/4)),((p.szX/4+1):0.25:(3*p.szX/4)));
%p.xCoord = interp1((50:140),p.xCoord(50:140),(50:0.25:140));

% p.zCoord = interp1((400:450),p.zCoord(400:450),(400:0.5:450));
% p.zCoord = p.zCoord(85:165);

% zoom for lateral resolution targets
% p.zCoord = interp1((128:256),p.zCoord(128:256),(128:0.5:256));
% p.zCoord = p.zCoord(85:165);


% p.szX = length(p.xCoord);
% p.szZ = length(p.zCoord)-1;
% p.nPoints = p.szX*p.szZ;

% Initialize full k-space beamformer
beamform = reconraw.DASBModeOffline(p);

%% DAS + CF computation
% tic
cRF = computeCRF(RFData(:,:,1),p);

idxtMTX = beamform.computeDASFullKSpace(cRF);
idxtMTX = reshape(idxtMTX,[p.numEl,p.na,p.nPoints]);
BMode = reshape(squeeze(sum(idxtMTX,[1,2])),[p.szZ,p.szX]);
% toc

%% United Coherence Factor
unitedCF = (sum(idxtMTX,[1,2]).^2)./(double(p.numEl*p.na)*sum(idxtMTX.^2,[1,2]));
unitedCF = reshape(squeeze(sum((unitedCF.*sum(idxtMTX,2))/p.na,1)),p.szZ,p.szX);


%% Evaluate Jerome CF
weightMTX = zeros(p.szX,p.na);
weightedBFM = zeros(size(idxtMTX));

for i = 1:p.nPoints
    kMTX = idxtMTX(:,:,i);
    kMTX2 = abs(kMTX).^2;
    
    weightMTX = abs(sum(kMTX,1).*(sum(kMTX,2))).^2 ./ ...
        (p.szX*p.na*sum(kMTX2,1).*sum(kMTX2,2));
    weightMTX(isnan(weightMTX)) = 0;

%     weightedBFM(i,:,:) = permute(weightMTX,[3,1,2]);
%    weightedBFM(:,:,i) = weightMTX./max(weightMTX,[],'all');
    weightedBFM(:,:,i) = weightMTX;
    
end

weightedBFM = idxtMTX.*(weightedBFM);
weightedBFM = reshape(squeeze(sum(weightedBFM,[1,2])),[p.szZ,p.szX]);
% CFBmode = reshape(squeeze(sum(weightedBFM,[1,2])),[p.szZ,p.szX]);

%% Test coherence factor 1
weightMTX1 = zeros(p.szX,p.na);
weightedBFM1 = zeros(size(idxtMTX));

for i = 1:p.nPoints
    kMTX = idxtMTX(:,:,i);
    kMTX2 = abs(kMTX).^2;
    
    
%     weightMTX1 = abs(sum(kMTX,1).*conj(sum(kMTX,2))) ./ ...
%       (sum(abs(kMTX),1).*sum(abs(kMTX),2));
% 

%     weightMTX1 = abs(sum(kMTX,1).*conj(sum(kMTX,2))).^2 ./ ...
%         (p.szX*p.na*sum(kMTX2.*conj(kMTX2),[1,2]));

%      weightMTX1 = abs(sum(kMTX,1).*(sum(kMTX,2))).^2 ./ ...
%          (p.szX*p.na*sum(abs(kMTX*conj(kMTX')).^2,[1,2]));

%      weightMTX1 = abs(p.na*sum(kMTX,1)+p.szX*sum(kMTX,2)).^2 ./ ...
%          (p.szX*p.na*(p.na*sum(kMTX2,1)+p.szX*sum(kMTX2,2)));

        % weightMTX1 = abs(sum(kMTX,1)/p.szX+sum(kMTX,2)/p.na).^2 ./ ...
        %  (p.szX*p.na*(p.na*sum(kMTX2,1)+p.szX*sum(kMTX2,2)));

      % weightMTX1 = abs(sum(kMTX*conj(kMTX'),[1,2])) ./ ...
      %     (p.szX*p.na*sum(abs(kMTX*conj(kMTX')),[1,2]));
      % 
      % weightMTX1 = abs(sum(kMTX,1).*(sum(kMTX,2))) ./ ...
      %    (p.numEl*p.na*sum(abs(kMTX),1).*sum(abs(kMTX),2));
    
      % weightMTX1 = abs(sum(kMTX,1).*conj(sum(kMTX,2))).^2 ./ ...
      %     (p.numEl*p.na*sum(kMTX2,1).*sum(kMTX2,2));

      kMTX2 = abs(kMTX).^4;
      weightMTX1 = abs(sum(kMTX,1).*conj(sum(kMTX,2))).^4 ./ ...
          (p.numEl*p.na*sum(kMTX2,1).*sum(kMTX2,2));

      % weightMTX1 = real(sum(kMTX,1).*conj(sum(kMTX,2))) ./ ...
      %    (p.numEl*p.na*sum(abs(kMTX),1).*sum(abs(kMTX),2));

      
      

    weightMTX1(isnan(weightMTX1)) = 0;
    
    weightedBFM1(:,:,i) = weightMTX1;
    
end

weightedBFM1 = idxtMTX.*(weightedBFM1);
weightedBFMtest = reshape(squeeze(sum(weightedBFM1,[1,2])),[p.szZ,p.szX]);


%% Test coherence factor 2
% weightMTX2 = zeros(p.szX,p.na);
% weightedBFM2 = zeros(size(idxtMTX));
% 
% r=0.4;
% h=fspecial('average',[2*round(r*p.szX/2)+1,2*round(r*p.na/2)+1]);
% 
% for i = 1:p.nPoints
%     kMTX = idxtMTX(:,:,i);
%     kMTX2 = abs(kMTX).^2;
%     
%     akMTX = imfilter(kMTX, h, "replicate");
%     akMTX2 = imfilter(kMTX2, h, "replicate");
%     weightMTX2 = (kMTX2-abs(akMTX).^2.)./abs(akMTX).^2;
%     weightMTX2(isnan(weightMTX2)) = 0;
% 
%      
%     weightedBFM2(:,:,i) = weightMTX2;
%     
% end
% 
% weightedBFM2 = idxtMTX.*weightedBFM2;
% weightedBFMtest = reshape(squeeze(sum(weightedBFM2,[1,2])),[p.szZ,p.szX]);
% 

%% LogCampress
DynRange=60;
BModeLog=10*log10(abs(BMode).^2);
BModeLog(BModeLog<(max(BModeLog,[],'all')-DynRange))=(max(BModeLog,[],'all')-DynRange);
weightedBFMLog=10*log10(abs(weightedBFM));
weightedBFMLog(weightedBFMLog<(max(weightedBFMLog,[],'all')-DynRange))=(max(weightedBFMLog,[],'all')-DynRange);
weightedBFMtestLog=10*log10(abs(weightedBFMtest));
weightedBFMtestLog(weightedBFMtestLog<(max(weightedBFMtestLog,[],'all')-DynRange))=(max(weightedBFMtestLog,[],'all')-DynRange);
unitedCFLog=10*log10(abs(unitedCF));
unitedCFLog(unitedCFLog<(max(unitedCFLog,[],'all')-DynRange))=(max(unitedCFLog,[],'all')-DynRange);


% Produce images
% figure(100)
% colormap(gray)
% subplot(1,4,1)
% %plotLogScaleImage(BMode(450:480,80:110))
% %imagesc(abs(BMode(450:480,80:110)).^0.1)
% %plotLogScaleImage(BMode)
% imagesc(BModeLog)
% %imagesc(abs(BMode).^0.15)
% title('Log BMode')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% subplot(1,4,2)
% %plotLogScaleImage(weightedBFM(450:480,80:110))
% %imagesc(abs(weightedBFM(450:480,80:110)).^0.1)
% %plotLogScaleImage(weightedBFM)
% imagesc(weightedBFMLog)
% %imagesc(abs(weightedBFM).^0.15)
% title('Log JCF')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% subplot(1,4,3)
% %plotLogScaleImage(weightedBFM(450:480,80:110))
% %imagesc(abs(weightedBFM(450:480,80:110)).^0.1)
% %plotLogScaleImage(weightedBFMtest)
% imagesc(weightedBFMtestLog)
% %imagesc(abs(weightedBFMtest).^0.15)
% title('Log Test')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% subplot(1,4,4)
% %plotLogScaleImage(unitedCF(450:480,80:110))
% %imagesc(abs(unitedCF(450:480,80:110)).^0.1)
% %plotLogScaleImage(unitedCF)
% imagesc(unitedCFLog)
% %imagesc(abs(unitedCF).^0.15)
% title('Log UCF')
% set(gca,'xtick',[])
% set(gca,'ytick',[])

figure
gammaFactor=0.15;
colormap(gray)
subplot(1,4,1)
%plotLogScaleImage(BMode(450:480,80:110))
%imagesc(abs(BMode(450:480,80:110)).^0.1)
%plotLogScaleImage(BMode)
imagesc(abs(BMode).^2.^gammaFactor)
title('Gamma BMode')
set(gca,'xtick',[])
set(gca,'ytick',[])
subplot(1,4,2)
%plotLogScaleImage(weightedBFM(450:480,80:110))
%imagesc(abs(weightedBFM(450:480,80:110)).^0.1)
%plotLogScaleImage(weightedBFM)
imagesc(abs(weightedBFM).^gammaFactor)
title('Gamma JCF')
set(gca,'xtick',[])
set(gca,'ytick',[])
subplot(1,4,3)
%plotLogScaleImage(weightedBFM(450:480,80:110))
%imagesc(abs(weightedBFM(450:480,80:110)).^0.1)
%plotLogScaleImage(weightedBFMtest)
imagesc(abs(weightedBFMtest).^gammaFactor)
title('Gamma Test')
set(gca,'xtick',[])
set(gca,'ytick',[])
subplot(1,4,4)
%plotLogScaleImage(unitedCF(450:480,80:110))
%imagesc(abs(unitedCF(450:480,80:110)).^0.1)
%plotLogScaleImage(unitedCF)
imagesc(abs(unitedCF).^gammaFactor)
title('Gamma UCF')
set(gca,'xtick',[])
set(gca,'ytick',[])


% figure
% hold off
% semilogy(abs(BMode(465,80:110))/max(abs(BMode(465,80:110)),[],'all'))
% hold on
% semilogy(abs(weightedBFM(465,80:110))/max(abs(weightedBFM(465,80:110)),[],'all'))
% 
% % zcoord=198;
% %zcoord=264;
% zcoord=202;
% % zcoord=465;
% zcoord=64;
% figure (10)
% hold on
% plot(10*log10(abs(BMode(zcoord,:)/max(BMode(zcoord,:)))))
% plot(10*log10(abs(weightedBFM(zcoord,:)/max(weightedBFM(zcoord,:)))))
% plot(10*log10(abs(weightedBFMtest(zcoord,:)/max(weightedBFMtest(zcoord,:)))))
% plot(10*log10(abs(unitedCF(zcoord,:)/max(unitedCF(zcoord,:)))))
% hold off
% legend("BMode","Jerome CF",'Test',"uCF")
% xlabel('Element Number')
% ylabel('Norm DB scale')
% 
% % 
% % xcoord=97;
% % figure
% % imagesc(abs(idxtMTX(:,:,p.szZ*(xcoord-1)+zcoord)))

% figure (20)
% plot(sum(abs(BMode),2)./sum(abs(BMode(2,:))))
% hold on
% plot(sum(abs(weightedBFM),2)./sum(abs(weightedBFM(2,:))))




