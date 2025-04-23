% registerULM

load('D:\Allen\Data\AZ02 Stroke ULM RC15gV\03-17-2025 baseline right eye\ULM processing results\Pairing and tracking all files\bubble_density_maps.mat')
load('D:\Allen\Data\AZ02 Stroke ULM RC15gV\04-15-25 1h left eye\ULM processing results\Pairing and tracking results\bubble_density_maps.mat', 'BDMs_AZ02_hour1')
load('D:\Allen\Data\AZ02 Stroke ULM RC15gV\04-18-2025 3d right eye\ULM processing results\Pairing and tracking results\bubble_density_maps.mat', 'BDMs_AZ02_day3')
%%
hour1BDM = BDMs_AZ02_hour1.BDM_LI_RSC;
day3BDM = BDMs_AZ02_day3.BDM_LI_RSC;
baselineBDM = BDMs_AZ02_baseline.BDM_LI_RSC;

% fixed = fixed ./ max(fixed, [], 'all');
% moving_orig = moving_orig ./ max(moving_orig, [], 'all');

% View histograms to see where we can set the cutoff. It's to remove those
% points where some false positive bubble is present for every frame, so
% it's really bright
% figure; hf = histogram(fixed(:), 10, 'BinWidth', 10)
% figure; hm = histogram(moving_orig, 10, 'BinWidth', 10)
%
cutoff = 500;
hour1BDM(hour1BDM > cutoff) = 0;
day3BDM(day3BDM > cutoff) = 0;
baselineBDM(baselineBDM > cutoff) = 0;

%% Adjust the histogram
% moving = imhistmatchn(moving_orig, fixed);
day3BDM = day3BDM;

%%
sv1 = size(hour1BDM);
sv2 = size(day3BDM);
figure
imshowpair(max(hour1BDM(:, :, sv1(3)/2 - 5 : sv1(3)/2 + 5) .^ 0.3, [], 3), max(day3BDM(:, :, sv2(3)/2 - 5 : sv2(3)/2 + 5) .^ 0.3, [], 3))

%% Base MIPs for the original data
mipPower = 1;
yr = 70:90; % y range for MIP
figure; imagesc(squeeze(max(hour1BDM(yr, :, :), [], 1))' .^ mipPower); colormap hot
figure; imagesc(squeeze(max(day3BDM(yr, :, :), [], 1))' .^ mipPower); colormap hot

%% Show the two volumes on top of each other
viewerUnregistered = viewer3d(BackgroundColor="black",BackgroundGradient="off");
volshow(hour1BDM .^ mipPower, Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
    Colormap=[1 0 1],Alphamap=1);
volshow(day3BDM .^ mipPower, Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
    Colormap=[0 1 0],Alphamap=1);

%% imregdemons test

[D, ird] = imregdemons(day3BDM, hour1BDM, 100, 'DisplayWaitbar', true, 'PyramidLevels', 3, 'AccumulatedFieldSmoothing', 3);

%% display results of imregdemons
% volumeViewer(ird .^ 0.3)
figure; imagesc(squeeze(max(ird(yr, :, :), [], 1))' .^ mipPower); colormap hot

% viewerRegistered = viewer3d(BackgroundColor="black",BackgroundGradient="off");
% volshow(fixed .^ 0.3, Parent=viewerRegistered,RenderingStyle="Isosurface", ...
%     Colormap=[1 0 1],Alphamap=1);
% volshow(moving .^ 0.3, Parent=viewerRegistered,RenderingStyle="Isosurface", ...
%     Colormap=[0 1 0],Alphamap=1);

%% Trying imregister

[optimizer, metric] = imregconfig('monomodal');
% optimizer.GradientMagnitudeTolerance = 1e-7;
optimizer.MaximumIterations = 10000;
optimizer.MinimumStepLength = 1e-6;

% Inputs: moving, fixed, transform type, optimizer, metric
% tic
% [day3BDM_registered, R_reg_day3BDM] = imregister(moving, fixed, 'rigid', optimizer, metric, 'DisplayOptimization', true);
% toc

tic
[baselineBDM_registered, R_reg_baselineBDM] = imregister(baselineBDM, hour1BDM, 'rigid', optimizer, metric, 'DisplayOptimization', true);
toc

% figure; imagesc(squeeze(max(moving_registered(yr, :, :), [], 1))' .^ 1); colormap hot
figure; imagesc(squeeze(max(baselineBDM_registered(yr, :, :), [], 1))' .^ 0.6); colormap hot


%% Trying imregtform

% test = imregtform(moving, fixed, 'rigid', optimizer, metric, 'DisplayOptimization', true);
movingRegistered = imwarp(day3BDM, moving_registered, "OutputView",imref3d(size(hour1BDM))); % apply the transformation
% movingRegistered = imwarp(vol2, test);
%%
% figure; imagesc(squeeze(max(test(yr, :, :), [], 1))' .^ 1); colormap hot
figure; imagesc(squeeze(max(hour1BDM(yr, :, :), [], 1))' .^ mipPower); colormap hot
% figure; imagesc(squeeze(max(movingReg.regVol.Voxels(yr, :, :), [], 1))'); colormap hot
figure; imagesc(squeeze(max(moving_registered(yr, :, :), [], 1))' .^ 0.5); colormap hot

% figure; imagesc(squeeze(max(movingReg1.regVol.Voxels(yr, :, :), [], 1))' .^ 0.3); colormap hot
% figure; imagesc(squeeze(max(movingRegistered(yr, :, :), [], 1))' .^ 0.3); colormap hot

%% Check the results of imregister

viewerThresholded = viewer3d(BackgroundColor = "black", BackgroundGradient="off");
volshow(hour1BDM .^ 1, Parent=viewerThresholded, RenderingStyle = "Isosurface", ...
    Colormap=[1 0 1],Alphamap=1);
volshow(moving_registered .^ 1, Parent=viewerThresholded, RenderingStyle = "Isosurface", ...
    Colormap=[0 1 0],Alphamap=1);

%%

volumeData = {baselineBDM_registered, hour1BDM, day3BDM_registered};

is = size(hour1BDM); % image size
uf = 2; % upsampling factor
volumeDataUpsampled = {imresize3(baselineBDM_registered, is .* uf), imresize3(hour1BDM, is .* uf), imresize3(day3BDM_registered, is .* uf)};

volumeDataGammaCompressed = {baselineBDM_registered .^ 0.7, hour1BDM .^ 1, day3BDM_registered .^ 0.9};
