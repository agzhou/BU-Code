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
% load([datapath, 'speed_maps_baseline'])
% load([datapath, 'speed_maps_hour1'])
% load([datapath, 'speed_maps_day1'])
% load([datapath, 'speed_maps_day4'])
% load([datapath, 'speed_maps_day7'])
load([datapath, 'speed_maps_baseline_SmoothedKF'])
load([datapath, 'speed_maps_hour1_SmoothedKF'])
load([datapath, 'speed_maps_day1_SmoothedKF'])
load([datapath, 'speed_maps_day4_SmoothedKF'])
load([datapath, 'speed_maps_day7_SmoothedKF'])
% SMs_raw = {SMs_AZ03_baseline.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ03_hour1.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ03_day1.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ03_day4.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ03_day7.SM_SmoothedKFConstrained_LI_Rfn};
SMs_raw = {SMs_AZ03_baseline.SM_SmoothedKF_LI_Rfn, SMs_AZ03_hour1.SM_SmoothedKF_LI_Rfn, SMs_AZ03_day1.SM_SmoothedKF_LI_Rfn, SMs_AZ03_day4.SM_SmoothedKF_LI_Rfn, SMs_AZ03_day7.SM_SmoothedKF_LI_Rfn};
% SMs_raw = {SMs_AZ03_baseline.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ03_hour1.SM_SmoothedKFConstrained_LI_Rfn};
% 
clearvars SMs_AZ03_baseline SMs_AZ03_hour1 SMs_AZ03_day1 SMs_AZ03_day4 SMs_AZ03_day7

%% Temporary data loading until I figure out the automatic stuff
% speed_map_name = 'SM_SmoothedKFConstrained_LI_Rfn';
load([datapath, 'speed_maps_baseline'])
load([datapath, 'speed_maps_hour1'])
load([datapath, 'speed_maps_day3'])
load([datapath, 'speed_maps_day7'])

% SMs_raw = {SMs_AZ02_baseline.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ02_hour1.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ02_day3.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ02_day7.SM_SmoothedKFConstrained_LI_Rfn};
SMs_raw = {SMs_AZ02_day7.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ02_baseline.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ02_hour1.SM_SmoothedKFConstrained_LI_Rfn, SMs_AZ02_day3.SM_SmoothedKFConstrained_LI_Rfn};

clearvars SMs_AZ02_baseline SMs_AZ02_hour1 SMs_AZ02_day3 SMs_AZ02_day7

%% Plot the raw speed maps

vcmap = colormap_ULM;

for mn = 1:length(SMs_raw)
    % figure; imagesc(squeeze(max(SMs_raw{mn}(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
    figure; imagesc(squeeze(max(SMs_raw{mn}(:, :, :), [], 3))'); colormap(vcmap); clim([0, 40])
end

%% rotate all speed maps manually so we can avoid smearing in our MIPs (think of some more robust way to do this)
% figure; imagesc(squeeze(max(SMs_raw{1}(:, :, :), [], 3))); colormap(vcmap); clim([0, 40])

rx = 0;
ry = 0;
rz = -7; % z rotation in degrees
% rz = -4; % z rotation in degrees

Rx = [1 0 0; 0 cosd(rx) -sind(rx); 0 sind(rx) cosd(rx)];
Ry = [cosd(ry) 0 sind(ry); 0 1 0; -sind(ry) 0 cosd(ry)];
Rz = [cosd(rz) -sind(rz) 0; sind(rz) cosd(rz) 0; 0 0 1];
R = Rz*Ry*Rx;

baseline_tform = rigidtform3d([rx, ry, rz], [0, 0, 0]); % No translation

SMs_IT = cell(size(SMs_raw)); % Speed maps with an Initial Transformation

% for mn = 1:length(SMs_raw) % Go through each map (number) 
% for mn = 2:length(SMs_raw)
for mn = 1
    SMs_IT{mn} = imwarp(SMs_raw{mn}, baseline_tform, 'cubic', "OutputView", imref3d(size(SMs_raw{mn})));
end

SMs_IT(2:end) = SMs_raw(2:end);
% figure; imagesc(squeeze(max(SMs_IT{1}(:, :, :), [], 3))); colormap(vcmap); clim([0, 40])

%% Plot the rotated base speed maps

vcmap = colormap_ULM;

% Coronal
% for mn = 1:length(SMs_IT)
%     % volumeSegmenter(SMs_IT{mn})
%     figure; imagesc(squeeze(max(SMs_IT{mn}(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
% end

% Axial
for mn = 1:length(SMs_IT)
    % volumeSegmenter(SMs_IT{mn})
    figure; imagesc(squeeze(max(SMs_IT{mn}(:, :, :), [], 3))); colormap(vcmap); clim([0, 40])
end

% generateTiffStack_multi({SMs_IT{2}}, region_size, vcmap, MIP_windowsize, speedRange)

%% Downsample before registering
dsf = 5; % Downsampling factor
SMs_IT_ds = cell(size(SMs_IT)); % Cell array of the downsampled speed maps

% Downsample the speed maps
for mn = 1:length(SMs_IT)
    SMs_IT_ds{mn} = imresize3(SMs_IT{mn}, 1/dsf, "Method", "cubic");
end

%% Plot the downsampled base speed maps

vcmap = colormap_ULM;

for mn = 1:length(SMs_IT_ds)
    % volumeSegmenter(SMs_IT{mn})
    figure; imagesc(squeeze(max(SMs_IT_ds{mn}(30:50, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
end

%%

%% Show two of the downsampled base speed map volumes on top of each other

compareVolumes(SMs_IT_ds{1}, SMs_IT_ds{2})
% compareVolumes(SMs_IT_ds{2}, SMs_IT_ds{3})

%% Register the bubble maps
SMs_reg_ds = cell(size(SMs_IT_ds)); % Initialize the registered downsampled speed maps
tforms_ds = cell(size(SMs_IT_ds)); % Initialize the transformation info for the registered, downsampled speed maps

if length(SMs_IT_ds) < 2
    error('There needs to be at least 2 volumes to run the registration')
end
SMs_reg_ds{1} = SMs_IT_ds{1}; % The first volume (baseline) does not need to be registered to itself

%%
% for mn = 2:length(SMs_IT_ds)
for mn = [2]
    [SMs_reg_ds{mn}, tforms_ds{mn}] = rigidRegTF(SMs_IT_ds{mn}, SMs_IT_ds{1});
    disp("Downsampled speed map #" + num2str(mn) + " registered.")
end

%% manual cheating
% test_tform = tforms_ds{mn};
% rxt = 0; ryt = 0; rzt = 0;
% Rxt = [1 0 0; 0 cosd(rxt) -sind(rxt); 0 sind(rxt) cosd(rxt)];
% Ryt = [cosd(ryt) 0 sind(ryt); 0 1 0; -sind(ryt) 0 cosd(ryt)];
% Rzt = [cosd(rzt) -sind(rzt) 0; sind(rzt) cosd(rzt) 0; 0 0 1];
% Rt = Rzt*Ryt*Rxt;
% 
% test_tform.R = Rt;
% test_tform.Translation = [0 0 0];
% % test_tform.Translation(2) = -tforms_ds{mn}.Translation(2);
% test_tform.Translation(1) = 30; 
% test = imwarp(SMs_reg_ds{2}, test_tform, 'cubic', "OutputView", imref3d(size(SMs_IT_ds{2})));
% 
% SMs_reg_ds_old = SMs_reg_ds;
% SMs_reg_ds{2} = test;
% 
% tforms_ds_old = tforms_ds;
% tforms_ds{2}.Translation = tforms_ds_old{2}.Translation + test_tform.Translation;

% Manually register visually

fixed = SMs_reg_ds{1};
moving1 = SMs_reg_ds{2};
moving2 = SMs_reg_ds{3};
moving3 = SMs_reg_ds{4};
moving4 = SMs_reg_ds{5};

medicalRegistrationEstimator
SMs_reg_ds_old = SMs_reg_ds;
tforms_ds_old = tforms_ds;

moving1_manual_reg = imwarp(moving1, movingReg1.tform, 'cubic', "OutputView", imref3d(size(moving1)));
SMs_reg_ds{2} = moving1_manual_reg;
tforms_ds{2}.A = movingReg1.tform.A * tforms_ds_old{2}.A;

moving2_manual_reg = imwarp(moving2, movingReg2.tform, 'cubic', "OutputView", imref3d(size(moving2)));
SMs_reg_ds{3} = moving2_manual_reg;
tforms_ds{3}.A = movingReg2.tform.A * tforms_ds_old{3}.A;

moving3_manual_reg = imwarp(moving3, movingReg3.tform, 'cubic', "OutputView", imref3d(size(moving3)));
SMs_reg_ds{4} = moving3_manual_reg;
tforms_ds{4}.A = movingReg3.tform.A * tforms_ds_old{4}.A;

moving4_manual_reg = imwarp(moving4, movingReg4.tform, 'cubic', "OutputView", imref3d(size(moving4)));
SMs_reg_ds{5} = moving4_manual_reg;
tforms_ds{5}.A = movingReg4.tform.A * tforms_ds_old{5}.A;

%% Compare the downsampled registration(s) to the fixed image

vcmap = colormap_ULM;

for mn = 1:length(SMs_reg_ds)
    % volumeSegmenter(SMs_IT{mn})
    figure; imagesc(squeeze(max(SMs_reg_ds{mn}(30:50, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
end
%%
compareVolumes(SMs_reg_ds{1}, SMs_reg_ds{2})
% compareVolumes(SMs_reg_ds{1}, test)

%% Use the downsampled transformation information and transform the full-sized data
tforms_reg = cell(size(tforms_ds)); % Initialize the transformation info for the registered, full-sized speed maps
SMs_reg = cell(size(SMs_reg_ds)); % Initialize the registered, full-sized speed maps
SMs_reg{1} = SMs_IT{1}; % The first volume (baseline) does not need to be registered to itself

for mn = 2:length(SMs_reg_ds)
% for mn = 2
    tforms_reg{mn} = tforms_ds{mn};
    tforms_reg{mn}.Translation = tforms_reg{mn}.Translation .* dsf; % Adjust the transformation to account for the prior downsampling

    SMs_reg{mn} = imwarp(SMs_IT{mn}, tforms_reg{mn}, 'cubic', "OutputView", imref3d(size(SMs_IT{mn}))); % Apply the full-sized transformation
end

%% MIPs of the final transformed data
vcmap = colormap_ULM;

for mn = 1:length(SMs_reg)
    % volumeSegmenter(SMs_reg{mn})
    figure; imagesc(squeeze(max(SMs_reg{mn}(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
end

%% Save MIPs

vcmap = colormap_ULM;
speedRange = [0, 40];
% speedRange = [0, 5];
MIP_windowsize = 50;
% region_size = [8.8, 8.8, 8];
total_padding = [80, 80, 80];
region_size = (1 + total_padding ./ (size(SMs_reg{1}) - total_padding)) .* [8.8, 8.8, 8];


% Bubble density map
% BDMs_registered = {baselineBDM_registered, hour1BDM_registered, day3BDM_registered, day7BDM};
% BDMs_registered_gamma = {baselineBDM_registered .^ 0.4, hour1BDM_registered .^ 0.5, day3BDM_registered .^ 0.5, day7BDM .^ 0.5};
% generateTiffStack_multi(BDMs_registered_gamma, region_size, 'hot', MIP_windowsize)

% Speed map
% SMs_registered = {SM_baseline_rotated, SM_hour1_registered};
% SMs_reg_outoforder = SMs_reg;
% SMs_reg = [SMs_reg(2:4), SMs_reg(1)];

for mn = 1:length(SMs_reg)
    % volumeSegmenter(SMs_reg{mn})
    figure; imagesc(squeeze(max(SMs_reg{mn}(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
end

generateTiffStack_multi(SMs_reg, region_size, vcmap, MIP_windowsize, speedRange)

generateTiffStack_multi({SMs_reg{1}, SMs_reg{2}}, region_size, vcmap, MIP_windowsize, speedRange)
% generateTiffStack_multi({SMs_reg{2}}, region_size, vcmap, MIP_windowsize, speedRange)
% generateTiffStack_multi({SMs_reg{2}, SMs_reg{4}}, region_size, vcmap, MIP_windowsize, speedRange)

% is = size(hour1SM); % image size
% uf = 2; % upsampling factor
% volumeDataUpsampled = {imresize3(baselineBDM_registered, is .* uf), imresize3(hour1SM, is .* uf), imresize3(day3BDM_registered, is .* uf)};

% Plot xy MIPs
figure; imagesc(squeeze(max(SMs_reg{4}(:, :, :), [], 3))'); colormap(vcmap)
figure; imagesc(squeeze(max(SMs_reg{3}(:, :, :), [], 3))'); colormap(vcmap)
figure; imagesc(squeeze(max(SMs_reg{2}(:, :, :), [], 3))'); colormap(vcmap)
figure; imagesc(squeeze(max(SMs_reg{1}(:, :, :), [], 3))'); colormap(vcmap)

%% Subtraction to accentuate the stroke core?
% % figure; imagesc(squeeze(max(SMs_reg{1}(300:500, :, :), [], 1))'); colormap(vcmap)
figure; imagesc(squeeze(max(SMs_reg{4}(300:500, :, :), [], 1))'); colormap(vcmap)
figure; imagesc(squeeze(max(SMs_reg{2}(300:500, :, :), [], 1))'); colormap(vcmap)


% baseline_minus_hour1 = SMs_reg{1} - SMs_reg{2};
day7_minus_hour1 = SMs_reg{4} - SMs_reg{2};

% volumeViewer(baseline_minus_hour1)
% figure; imagesc(squeeze(max(baseline_minus_hour1(300:500, :, :), [], 1))'); colormap(vcmap)
figure; imagesc(squeeze(max(day7_minus_hour1(300:500, :, :), [], 1))'); colormap(vcmap)

% Could try binarizing and subtracting...

%% Hard coded stuff for figure making of AZ02
savepath = "D:\Allen\Official Stuff\June 2025 ultrasound grant\Figures\02 ULM slices pre and post stroke (AZ02)\";

actualSize = region_size;
hwRatio_xz = actualSize(3) / actualSize(2); % height to width ratio for the full image
xz_slice_yrange = 222:272;
xz_slice_xrange = 480:960;
for mn = 1:length(SMs_reg)
    figure; imagesc(squeeze(max(SMs_reg{mn}(xz_slice_yrange, xz_slice_xrange, :), [], 1))'); colormap(vcmap); clim([0, 40])
    axis tight
    ax = gca;
    ax.PlotBoxAspectRatio = [1, hwRatio_xz ./ (length(xz_slice_xrange)/size(SMs_reg{mn}, 2)), 1];

end

%% Save the images

for mn = 1:length(SMs_reg)
    % Get the axes as an image
    figure(mn)
    axh = gca;
    fr = getframe(axh);

    % Save the image
    imwrite(fr.cdata, savepath + num2str(mn) + ".png");
end
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
    % [tform] = imregtform(img, fixed, 'translation', optimizer, metric, 'DisplayOptimization', true, 'PyramidLevels', 3);
    [tform] = imregtform(img, fixed, 'translation', optimizer, metric, 'PyramidLevels', 3);
    img_reg = imwarp(img, tform, 'cubic', "OutputView", imref3d(size(fixed)));
    disp('Registration done')
    toc
end

% Get the registered image and the rigid transformation for registering
% 'img' to 'fixed'
function [img_reg, tform] = rigidRegTF(img, fixed)
    [optimizer, metric] = imregconfig('monomodal');
    % optimizer.GradientMagnitudeTolerance = 1e-5;
    optimizer.GradientMagnitudeTolerance = 1e-12;
    % optimizer.MaximumIterations = 10000;
    % optimizer.MaximumIterations = 20000;
    optimizer.MaximumIterations = 50000;
    % optimizer.MaximumIterations = 10;
    optimizer.MinimumStepLength = 1e-9;
    % optimizer.MinimumStepLength = 1e-4;
    % optimizer.MaximumStepLength = .5;
    optimizer.MaximumStepLength = .002;
    
    % Inputs: moving, fixed, transform type, optimizer, metric

    tic
    [tform] = imregtform(img, fixed, 'rigid', optimizer, metric, 'DisplayOptimization', true, 'PyramidLevels', 3);
    % [tform] = imregtform(img, fixed, 'rigid', optimizer, metric, 'PyramidLevels', 3);
    img_reg = imwarp(img, tform, 'cubic', "OutputView", imref3d(size(fixed)));
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
        registeredData{i} = imwarp(data{i}, tforms{i}, 'cubic', "OutputView", imref3d(size(fixed)));
    end
end