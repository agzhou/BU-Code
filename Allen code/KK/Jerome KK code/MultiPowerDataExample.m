% Variable Power Data
%% Initialize file location
clearvars
% close all

% Extract Current Path
currentDir = matlab.desktop.editor.getActiveFilename; 
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

% Boas Lab data multiframe
% dataFile{1} = dataFilePath + "Boas Lab Data\DeepTargetsMultiPower.mat";
dataFile{1} = dataFilePath + "Boas Lab Data\Calf15AngMultiPower.mat";
% dataFile{1} = dataFilePath + "Boas Lab Data\Calf45AngMultiPower.mat";
% dataFile{1} = dataFilePath + "Boas Lab Data\AbdomenA15AngMultiPower.mat";
% dataFile{1} = dataFilePath + "Boas Lab Data\AbdomenB15AngMultiPower.mat";
% dataFile{1} = dataFilePath + "Boas Lab Data\carotidCrossSec15AngMultiPower.mat";
% dataFile{1} = dataFilePath + "Boas Lab Data\carotidLongitudinalSec15AngMultiPower.mat";

% dataFile{1} = dataFilePath + "Boas Lab Data\cardiacParasternalShortAxis15Ang.mat";
filetype = 14;

%% Load Data
[p,RFData] = initParams(dataFile,filetype);
p.tShift = 0.0;
p.txPL = zeros(length(p.TXangle),1);
nFrames = size(RFData,3);

%% Init Coords and BFM
p2 = p;
% p2 = computeNewGrid(p,[p.xCoord(1),p.xCoord(end)],[p.zCoord(1),0.0754],p.szX,875);
% p2 = computeNewGrid(p,[p.xCoord(1),p.xCoord(end)],[p.zCoord(1),p.zCoord(end)*2],p.szX,p.szZ*2);
beamform = reconraw.DASBModeOffline(p2);

%% Compute BMode of multiframe results
Recon = zeros(p2.szZ,p2.szX,nFrames);
for n = 1:nFrames
    cRF = computeCRF(RFData(:,:,n),p2);
    Recon(:,:,n) = beamform.computeDAScrfBMode(cRF);
end

genCustomSlider(@plotBioSlider,10*log10(abs(Recon)),p2.xCoord,p2.zCoord)

%% Compute JCF of low and high power

% Low Power Data
cRF = computeCRF(RFData(:,:,1),p2);
idxtMTX = beamform.computeDASFullKSpace(cRF);
idxtMTX = permute(reshape(idxtMTX.',[p2.numEl,p2.nPoints,p2.na]),[1,3,2]);

JCF = JCFcombinedDenom(idxtMTX,p2);

figure
plotLogScaleImage(p.xCoord,p.zCoord,JCF)
title("JCF")
axis image

% High power data
cRF = computeCRF(RFData(:,:,4),p2);
idxtMTX = beamform.computeDASFullKSpace(cRF);
idxtMTX = permute(reshape(idxtMTX.',[p2.numEl,p2.nPoints,p2.na]),[1,3,2]);

JCF = JCFcombinedDenom(idxtMTX,p2);

figure
plotLogScaleImage(p.xCoord,p.zCoord,JCF)
title("JCF")
axis image

%% Helper Functions
function plotBioSlider(vars)
    % Plots the image
    val = round(get(vars.slider1_handle,'Value'));
    imagesc(vars.Axes,vars.field1*1e3,vars.field2*1e3,vars.Data(:,:,val));
    colormap(vars.Axes,'gray');
    set(vars.Label,'String',num2str(val));
    axis(vars.Axes,'image');
end