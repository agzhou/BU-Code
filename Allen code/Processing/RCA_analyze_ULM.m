%% Description
% ULM data analysis:
% Take reconstructed data from individual buffers/batches/superframes, each with
% some number of frames/subframes, each containing N acquisitions/angles.

% Re-sample with a rolling method to get more effective frames.

% SVD to separate the bubble signals from the tissue signal and other
% clutter

% ...
%% Use parallel processing for speed
% https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox

pp = gcp('nocreate');
if isempty(pp)
    % There is no parallel pool
    parpool LocalProfile1

end

%% Load parameters and make folder for saving the processed data

% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', 'G:\Allen\Data\');
    load([params_pathname, params_filename])
end

% Get data path of the reconstructed IQ data
datapath = uigetdir('G:\Allen\Data\', 'Select the IQ data path');
datapath = [datapath, '\'];

% Load Verasonics reconstruction parameters: datapath\PData.mat
if ~exist('PData', 'var')
    load([datapath, 'PData.mat'])
end

% Prompt for parameter user input
parameterPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound', 'Image refinement factor - x', 'Image refinement factor - y', 'Image refinement factor - z', 'XC Adaptive Threshold Factor', 'x pixel spacing [um]', 'y pixel spacing [um]', 'z pixel spacing [um]'};
parameterDefaults = {'1', '', '10', '150', '2', '2', '2', '0.2', num2str(PData.PDelta(1) * P.wl * 1e6), num2str(PData.PDelta(2) * P.wl * 1e6), num2str(PData.PDelta(3) * P.wl * 1e6)};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});
numFiles = endFile - startFile + 1;
sv_threshold_lower = str2double(parameterUserInput{3});
sv_threshold_upper = str2double(parameterUserInput{4});
imgRefinementFactor = [str2double(parameterUserInput{5}), str2double(parameterUserInput{6}), str2double(parameterUserInput{7})];
if any(floor(imgRefinementFactor) ~= imgRefinementFactor) || any(imgRefinementFactor < 1)
    error('Image refinement factors must be whole numbers')
end
XCThresholdFactor = str2double(parameterUserInput{8});
xpix_spacing = str2double(parameterUserInput{9});
ypix_spacing = str2double(parameterUserInput{10});
zpix_spacing = str2double(parameterUserInput{11});

% IQfolderName = 'IQ Data - Verasonics Recon\'; % 'IQ data\'
saveFolderName = 'Processed Data 03-17-2025\';
% savepath = [datapath, saveFolderName];
% mkdir([datapath, saveFolderName])
savepath = ['G:\Allen\Data\03-17-2025 AZ02 ULM\RC15gV\run 2 right eye\', saveFolderName];
mkdir(savepath)
% savepath = 'F:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\FMAS Processed Data\';

filename_structure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

addpath('C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Allen code\Processing\normxcorr3.m')
addpath('C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Allen code\Processing\toolbox_nlmeans_version2')
%% Parameters for processing the data
% Define various processing parameters
% Singular value thresholds
% sv_threshold_lower = 10;
% sv_threshold_upper = 150;

% % Region of interest
xrange = int16(1:80);
yrange = int16(1:80);
zrange = int16(36:142);

% framerange = 1:200;
% framerange = 1:size(IQf, 3);
% range = {xrange, yrange, zrange, framerange};
range = {xrange, yrange, zrange};

% % Image refinement and localization parameters
% irfc = 2;
% imgRefinementFactor = [2, 2, 2]; % z, x pixel refinement factor
% imgRefinementFactor = ones(1, 3) .* irfc;

% xpix_spacing = P.Trans.spacingMm / 1e3;
% ypix_spacing = P.Trans.spacingMm / 1e3;
% zpix_spacing = P.wl / 2;
% imgRefinementFactor = [irfc * xpix_spacing/zpix_spacing, irfc * ypix_spacing/zpix_spacing, irfc];

% XCThresholdFactor = 0.2;
 
% Load and refine simulated PSF
if ~exist('PSF', 'var')
    load('G:\Allen\Data\RC15gV PSF sim\PSF.mat', 'PSF')
    % figure; imagesc(squeeze(abs(PSF(40, :, :)))')
end

% PSFs = PSF(190:210, 118:138, :); % PSF section
PSFs = PSF(30:50, 30:50, 92:110); % PSF section
refPSF = imresize3(PSFs, [size(PSFs, 1) * imgRefinementFactor(1), size(PSFs, 2) * imgRefinementFactor(2), size(PSFs, 3) * imgRefinementFactor(3)]);
% volumeViewer(abs(refPSF))

% allCenters = {};
%% resampling size test
% s = size(IQ);
% test = zeros([s(1:3) .* 4, s(4)]);
%% Process the data
% tic
% for filenum = startFile:endFile
for filenum = 201
    tic
%     load([datapath, IQfolderName, filename_structure, num2str(filenum), '.mat'])  % load each reconstructed buffer/batch/superframe
    load([datapath, filename_structure, num2str(filenum), '.mat'])  % load each reconstructed buffer/batch/superframe
%     IQr = LA_rollingFrames(IQ);                                                 % rolling method to get more effective frames
    
    IQ = squeeze(IData + 1i .* QData);   % Combine I and Q, which are saved separately. It's easier to save the big reconstructed data with savefast, which doesn't support complex values. The data is already a coherent sum.
    clear IData QData
    
%     if filenum == 1
        [xp, yp, zp, nf] = size(IQ);
        range{4} = int16(1:nf); % set frame range after rolling on the first file
%     end

    % SVD proc part 1
%     tic
    [PP, EVs, V_sort] = getSVs2D(IQ);
    disp('SVs decomposed')
%     toc
    % SVD proc part 2
%     tic
    [IQf] = applySVs2D(IQ, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')

    clear PP EVs V_sort

%     save([savepath, 'Filtered-Data-', num2str(filenum)], 'IQr', 'PP', 'EVs', 'V_sort', 'IQf', "-v6")

%     [centers, refIQs, XC] = localizeBubbles3D(IQf, refPSF, range, imgRefinementFactor, XCThreshold);
%     [centers, ~, ~, XCThresholdAdaptive] = localizeBubbles3D(IQf, refPSF, range, imgRefinementFactor, XCThresholdFactor);
    [centers, refIQs, XC, XCThresholdAdaptive] = localizeBubbles3D(IQf, refPSF, range, imgRefinementFactor, XCThresholdFactor);

%     [coords, img_size, XCThresholdsAdaptive] = localizeBubbles3D_chunk(IQf, refPSF, range, imgRefinementFactor, XCThresholdFactor);

    %     save([savepath, 'IQf-', num2str(filenum)], 'IQf', "-v6")

%     save([savepath, 'dataproc-', num2str(filenum)], 'IQf', 'centroidCoordinates', "-v6")
    savefast([savepath, 'centers-', num2str(filenum)], 'centers', 'XCThresholdAdaptive')
%     savefast([savepath, 'coords-', num2str(filenum)], 'coords', 'img_size', 'XCThresholdsAdaptive')

%     allCenters = [allCenters; centers];
    disp(strcat("Centroid finding done: file ", num2str(filenum)))
    toc
end
save([savepath, 'proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'PSF', 'range', 'imgRefinementFactor', 'XCThresholdFactor')
% save([savepath, 'allCenters'], 'allCenters', "-v7.3")
% toc
