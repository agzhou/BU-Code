% registerULM

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

%%
load([datapath, 'speed_maps_baseline_SmoothedKF'])
load([datapath, 'speed_maps_hour1_SmoothedKF'])
load([datapath, 'speed_maps_day3_SmoothedKF'])
load([datapath, 'speed_maps_day7_SmoothedKF'])

SMs_raw = {SMs_AZ02_day7.SM_SmoothedKF_LI_Rfn, SMs_AZ02_baseline.SM_SmoothedKF_LI_Rfn, SMs_AZ02_hour1.SM_SmoothedKF_LI_Rfn, SMs_AZ02_day3.SM_SmoothedKF_LI_Rfn};

% clearvars SMs_AZ02_baseline SMs_AZ02_hour1 SMs_AZ02_day3 SMs_AZ02_day7

%%
% load([datapath, 'speed_maps_baseline'])
% load([datapath, 'speed_maps_hour1'])
% 
% BDM_baseline = SMs_AZ03_baseline.SM_SmoothedKFConstrained_counter;
% SM_baseline = SMs_AZ03_baseline.SM_SmoothedKFConstrained_LI_Rfn;
% 
% BDM_hour1 = SMs_AZ03_hour1.SM_SmoothedKFConstrained_counter;
% SM_hour1 = SMs_AZ03_hour1.SM_SmoothedKFConstrained_LI_Rfn;
% 
% clearvars SMs_AZ03_baseline SMs_AZ03_hour1
%% rotate baseline manually
rx = 0;
ry = 0;
rz = -8; % z rotation in degrees

Rx = [1 0 0; 0 cosd(rx) -sind(rx); 0 sind(rx) cosd(rx)];
Ry = [cosd(ry) 0 sind(ry); 0 1 0; -sind(ry) 0 cosd(ry)];
Rz = [cosd(rz) -sind(rz) 0; sind(rz) cosd(rz) 0; 0 0 1];
R = Rz*Ry*Rx;

% baseline_tform = rigidtform3d(R, [0, 0, 0])
baseline_tform = rigidtform3d([rx, ry, rz], [0, 0, 0]);
SM_baseline_rotated = imwarp(SM_baseline, baseline_tform, "OutputView", imref3d(size(SM_baseline)));
SM_hour1_rotated = imwarp(SM_hour1, baseline_tform, "OutputView", imref3d(size(SM_hour1)));

%%
volumeSegmenter(SM_baseline)
volumeSegmenter(SM_baseline_rotated)
volumeSegmenter(SM_hour1)
volumeSegmenter(SM_hour1_rotated)
%%
cmap = colormap_ULM;
figure; imagesc(squeeze(max(SM_baseline_rotated(400:600, :, :), [], 1))'); colormap(cmap); clim([0, 40])
figure; imagesc(squeeze(max(SM_hour1_rotated(400:600, :, :), [], 1))'); colormap(cmap); clim([0, 40])

%% Downsample before registering
dsf = 5; % Downsampling factor
SM_baseline_rotated_ds = imresize3(SM_baseline_rotated, 1/dsf);
SM_hour1_rotated_ds = imresize3(SM_hour1_rotated, 1/dsf);

%% Show two volumes on top of each other

compareVolumes(SM_baseline_rotated_ds, SM_hour1_rotated_ds)

%% Register the bubble maps

% [SM_hour1_registered, tf_hour1_SM] = rigidRegTF(SM_hour1_rotated, SM_baseline_rotated);
% [SM_hour1_registered, tf_hour1_SM] = translRegTF(SM_hour1_rotated, SM_baseline_rotated);

[SM_hour1_ds_registered, tf_hour1_ds_SM] = rigidRegTF(SM_hour1_rotated_ds, SM_baseline_rotated_ds);


%% Compare the downsampled registration to the fixed image
% figure; imagesc(squeeze(max(SM_baseline_rotated(400:600, :, :), [], 1))'); colormap(cmap); clim([0, 40])
% figure; imagesc(squeeze(max(SM_hour1_registered(400:600, :, :), [], 1))'); colormap(cmap); clim([0, 40])
figure; imagesc(squeeze(max(SM_baseline_rotated_ds(80:120, :, :), [], 1))'); colormap(cmap); clim([0, 40])
figure; imagesc(squeeze(max(SM_hour1_ds_registered(80:120, :, :), [], 1))'); colormap(cmap); clim([0, 40])

% compareVolumes(SM_baseline_rotated_ds, SM_hour1_ds_registered)

%% Use the downsampled transformation information
tf_hour1_SM = tf_hour1_ds_SM;
tf_hour1_SM.Translation = tf_hour1_SM.Translation .* dsf;
SM_hour1_registered = imwarp(SM_hour1_rotated, tf_hour1_SM, "OutputView", imref3d(size(SM_baseline_rotated)));

%% MIPs of the final transformed data
cmap = colormap_ULM;
figure; imagesc(squeeze(max(SM_baseline_rotated(400:600, :, :), [], 1))'); colormap(cmap); clim([0, 40])
figure; imagesc(squeeze(max(SM_hour1_registered(400:600, :, :), [], 1))'); colormap(cmap); clim([0, 40])


%% Save MIPs

cmap = colormap_ULM;
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
generateTiffStack_multi(SMs_registered, region_size, cmap, MIP_windowsize, speedRange)

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