clearvars
close all

%% Load Params

% Extract Path to current file location. These three lines must be run from
% inside this script and the script must be saved in one of the Ultrasound
% subfolders.
currentDir = matlab.desktop.editor.getActiveFilename; 
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

% Now define the filename and path. The dataset does not have to be saved
% on the drive. If it is not, simply replace the right hand side of the
% dataFile{1} definition with the complete file path and file name.

% Simulation Data
dataFile{1} = dataFilePath + "Jerome Data\Sim15_Point.mat";
filetype = 2;   % use this filetype if using the C++ beamformer. See initParams comments/documentation for other options for filetype

[p,RFData] = initParams(dataFile,filetype);


%% Beamforming Demonstrations

% This line defines the beamformer class. This beamformer class contains
% methods for doing traditional Delay and Sum Beamforming. Other beamformer
% classes will contain methods for doing other types of beamforming in C++, e.g.
% DPC, DSI, DMAS etc.
beamform = reconraw.DASBModeOffline(p);

% Compute the complex version of the RFData
cRF = computeCRF(RFData,p);

% These are the implementations of various methods you can call from the
% DAS class

tic; Recon1 = beamform.computeDASBMode(RFData); toc  % Using raw, real RF data. Ideally should implement fourier transforms in C++, but not yet working/implemented fully.
tic; ReconC = beamform.computeDAScrfBMode(cRF); toc  % Using hilbert transformed RF data.
tic; idxtMTX = beamform.computeDASFullKSpace(cRF); toc  % Generates uncompounded data.
idxtMTX = reshape(idxtMTX,[p.numEl,p.na,p.nPoints]);    % Need this line to conver to final dimensions of [num receive elements, num transmits, num pixels]

%% Plotting Demonstrations

% Plotting using log scale plotting function
figure
plotLogScaleImage(ReconC)   % If you want your axes to be labelled with indices, use this syntax

figure
plotLogScaleImage(p.xCoord*1e3,p.zCoord*1e3,ReconC) % If you want your axes to be labelled with coordinates, use this syntax


% Plotting using gamma scale plotting function
figure
plotGammaScaleImage(ReconC,0.5) % Follows similar syntax as log scale. The extra parameter defines the power the normalized magnitude would be raised to.

figure
plotGammaScaleImage(p.xCoord*1e3,p.zCoord*1e3,ReconC,0.1) % another gamma compression factor for example


%% Slider demonstrations
% Basic slider
ReconAng = reshape(squeeze(sum(idxtMTX,1)).',[p.szZ,p.szX,p.na]); % First, we make a 3D matrix of the beamformed data for each plane wave angle for this demonstration
genSliderV2(log10(abs(ReconAng))) % Then we pass the 3D matrix to genSliderV2. The matrix must be 3D and must not be complex, hence the abs().

% Custom slider (if you want multiple sub plots, Titles, axes labels,
% special plot styles, some other function rather than imagesc, etc.)

% TODO: forthcoming example from Nikunj

%% Export Nice Figure Demonstrations

% TODO: forthcoming examples from Nikunj
