% registerULM

%% Description:
%       Load in 3D ULM speed maps, apply an initial transformation
%       manually, and then register all of them to the first input.
%       The registration is done on a downsampled version of the speed 
%       maps for speed, and then the transformation is used on the
%       full-sized data.

%% Get data path
% Get data path of the speed map data
datapath = uigetdir('D:\Allen\Data\', 'Select the ULM speed map data path');
datapath = [datapath, '\'];

% filenames = {dir(datapath).name};
% % filenames = {dir(fullfile(datapath, '*.m')).name};
% 
% % filenames_trim = {};
% % for i = 1:length(filenames)
% %     if ~(filenames{i} == '.' | filenames{i} == '..')
% %     end
% % end
% [ind, ~] = listdlg('PromptString', {'Select the files to register'}, ...
%     'SelectionMode', 'multiple', 'ListString', filenames);

%% Temporary data loading until I figure out the automatic stuff
load([datapath, 'speed_maps_baseline'])
load([datapath, 'speed_maps_hour1'])
load([datapath, 'speed_maps_day1'])

SMs_raw = {SMs_AZ03_baseline.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ03_hour1.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ03_day1.SM_SmoothedKFConstrained_LI_Rfn};

clearvars SMs_AZ03_baseline SMs_AZ03_hour1 SMs_AZ03_day1

%% rotate all speed maps manually so we can avoid smearing in our MIPs (think of some more robust way to do this)
rx = 0;
ry = 0;
rz = -8; % z rotation in degrees

Rx = [1 0 0; 0 cosd(rx) -sind(rx); 0 sind(rx) cosd(rx)];
Ry = [cosd(ry) 0 sind(ry); 0 1 0; -sind(ry) 0 cosd(ry)];
Rz = [cosd(rz) -sind(rz) 0; sind(rz) cosd(rz) 0; 0 0 1];
R = Rz*Ry*Rx;

baseline_tform = rigidtform3d([rx, ry, rz], [0, 0, 0]); % No translation

SMs_IT = cell(size(SMs_raw)); % Speed maps with an Initial Transformation

for mn = 1:length(SMs_raw) % Go through each map (number) 
    SMs_IT{mn} = imwarp(SMs_raw{mn}, baseline_tform, "OutputView", imref3d(size(SMs_raw{mn})));
end

%% Plot the rotated base speed maps

vcmap = colormap_ULM;

for mn = 1:length(SMs_IT)
    % volumeSegmenter(SMs_IT{mn})
    figure; imagesc(squeeze(max(SMs_IT{mn}(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
end

%% Downsample before registering
dsf = 5; % Downsampling factor
SMs_IT_ds = cell(size(SMs_IT)); % Cell array of the downsampled speed maps

% Downsample the speed maps
for mn = 1:length(SMs_IT)
    SMs_IT_ds{mn} = imresize3(SMs_IT{mn}, 1/dsf);
end

%% Show two of the downsampled base speed map volumes on top of each other

compareVolumes(SMs_IT_ds{1}, SMs_IT_ds{2})

%% Register the bubble maps
SMs_reg_ds = cell(size(SMs_IT_ds)); % Initialize the registered downsampled speed maps
tforms_ds = cell(size(SMs_IT_ds)); % Initialize the transformation info for the registered, downsampled speed maps

if length(SMs_IT_ds) < 2
    error('There needs to be at least 2 volumes to run the registration')
end
SMs_reg_ds{1} = SMs_IT_ds{1}; % The first volume (baseline) does not need to be registered to itself

for mn = 2:length(SMs_IT_ds)
    [SMs_reg_ds{mn}, tforms_ds{mn}] = rigidRegTF(SMs_IT_ds{mn}, SMs_IT_ds{1});
end

%% Compare the downsampled registration to the fixed image
% figure; imagesc(squeeze(max(SM_baseline_rotated(400:600, :, :), [], 1))'); colormap(cmap); clim([0, 40])
% figure; imagesc(squeeze(max(SM_hour1_registered(400:600, :, :), [], 1))'); colormap(cmap); clim([0, 40])
figure; imagesc(squeeze(max(SM_baseline_rotated_ds(80:120, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
figure; imagesc(squeeze(max(SM_hour1_ds_registered(80:120, :, :), [], 1))'); colormap(vcmap); clim([0, 40])

% compareVolumes(SM_baseline_rotated_ds, SM_hour1_ds_registered)

%% Use the downsampled transformation information and transform the full-sized data
tforms_reg = cell(size(tforms_ds)); % Initialize the transformation info for the registered, full-sized speed maps
SMs_reg = cell(size(SMs_reg_ds)); % Initialize the registered, full-sized speed maps
SMs_reg{1} = SMs_IT{1}; % The first volume (baseline) does not need to be registered to itself

for mn = 2:length(SMs_reg_ds)
    tforms_reg{mn} = tforms_ds{mn};
    tforms_reg{mn}.Translation = tforms_reg{mn}.Translation .* dsf; % Adjust the transformation to account for the prior downsampling

    SMs_reg{mn} = imwarp(SMs_IT{mn}, tforms_reg{mn}, "OutputView", imref3d(size(SMs_IT{mn}))); % Apply the full-sized transformation
end

%% MIPs of the final transformed data
vcmap = colormap_ULM;
% figure; imagesc(squeeze(max(SM_baseline_rotated(400:600, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
% figure; imagesc(squeeze(max(SM_hour1_registered(400:600, :, :), [], 1))'); colormap(vcmap); clim([0, 40])

for mn = 1:length(SMs_reg)
    % volumeSegmenter(SMs_reg{mn})
    figure; imagesc(squeeze(max(SMs_reg{mn}(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
end

%% Save MIPs

vcmap = colormap_ULM;
speedRange = [0, 40];
% speedRange = [0, 5];
MIP_windowsize = 50;
region_size = [8.8, 8.8, 8];

% Bubble density map
% BDMs_registered = {baselineBDM_registered, hour1BDM_registered, day3BDM_registered, day7BDM};
% BDMs_registered_gamma = {baselineBDM_registered .^ 0.4, hour1BDM_registered .^ 0.5, day3BDM_registered .^ 0.5, day7BDM .^ 0.5};
% generateTiffStack_multi(BDMs_registered_gamma, region_size, 'hot', MIP_windowsize)

% Speed map
SMs_registered = {SM_baseline_rotated, SM_hour1_registered};
generateTiffStack_multi(SMs_registered, region_size, vcmap, MIP_windowsize, speedRange)

% is = size(hour1SM); % image size
% uf = 2; % upsampling factor
% volumeDataUpsampled = {imresize3(baselineBDM_registered, is .* uf), imresize3(hour1SM, is .* uf), imresize3(day3BDM_registered, is .* uf)};


%% Helper functions
% Get the registered image and the translation transformation for registering
% 'img' to 'fixed'
function [img_reg, tform] = translRegTF(img, fixed)
    [optimizer, metric] = imregconfig('monomodal');
    % optimizer.GradientMagnitudeTolerance = 1e-7;
    optimizer.MaximumIterations = 1000; %%%%%%%%%%
    % optimizer.MaximumIterations = 10;
    optimizer.MinimumStepLength = 1e-5;
    % optimizer.MaximumStepLength = 1;
    
    % Inputs: moving, fixed, transform type, optimizer, metric

    tic
    [tform] = imregtform(img, fixed, 'translation', optimizer, metric, 'DisplayOptimization', true, 'PyramidLevels', 3);
    img_reg = imwarp(img, tform, "OutputView", imref3d(size(fixed)));
    disp('Registration done')
    toc
end

% Get the registered image and the rigid transformation for registering
% 'img' to 'fixed'
function [img_reg, tform] = rigidRegTF(img, fixed)
    [optimizer, metric] = imregconfig('monomodal');
    % optimizer.GradientMagnitudeTolerance = 1e-7;
    optimizer.MaximumIterations = 1000; %%%%%%%%%%
    % optimizer.MaximumIterations = 10;
    optimizer.MinimumStepLength = 1e-5;
    % optimizer.MaximumStepLength = 10;
    
    % Inputs: moving, fixed, transform type, optimizer, metric

    tic
    [tform] = imregtform(img, fixed, 'rigid', optimizer, metric, 'DisplayOptimization', true, 'PyramidLevels', 3);
    img_reg = imwarp(img, tform, "OutputView", imref3d(size(fixed)));
    disp('Registration done')
    toc
end

function compareVolumes(vol1, vol2) % Can change this so it has a cell array input and goes through more than 2 volumes

    viewerThresholded = viewer3d(BackgroundColor = "black", BackgroundGradient="off");
    volshow(vol1 .^ 1, Parent=viewerThresholded, RenderingStyle = "Isosurface", ...
        Colormap=[1 0 1],Alphamap=1);
    volshow(vol2 .^ 1, Parent=viewerThresholded, RenderingStyle = "Isosurface", ...
        Colormap=[0 1 0],Alphamap=1);

end

% Apply previously obtained transforms 'tforms' to register each volume in the cell
% array 'data' to the fixed volume 'fixed'. Output is a cell array of the
% registered data.
function [registeredData] = applyTransforms(data, tforms, fixed)
    registeredData = cell(size(data));
    for i = 1:length(data)
        registeredData{i} = imwarp(data{i}, tforms{i}, "OutputView", imref3d(size(fixed)));
    end
end