%% Add paths for reading nrrd files
% addpath('C:\Users\agzhou\Downloads\nrrdread')
addpath('J:\Allen Atlas CCFv3\nrrdread')
%% Load the data and atlas
% Load fUS data (template map)
% load('H:\Ultrasound data from 7-21-2025\10-21-2025 AZ01 FC RCA\Test proc 10-23-25\PDI average across sf 100-500\fUSmap_50um.mat')
load('J:\Ultrasound data from 09-x-2025\10-21-2025 RCA US FC\Test proc 10-23-25\fUSmap_50um.mat')

% Load Allen atlas (10 um voxel size)
% [atlas, atlas_metadata] = nrrdread('C:\Users\agzhou\Downloads\average_template_10.nrrd');
% [atlas] = nrrdread('C:\Users\agzhou\Downloads\average_template_10.nrrd');
[atlas, metadata] = nrrdread('J:\Allen Atlas CCFv3\average_template_10.nrrd');
% [atlas, atlas_metadata] = nrrdread('D:\Allen\Data\Allen Atlas CCFv3\Annotated 25\annotation_25.nrrd');
% atlas_double = permute(double(atlas), [1, 3, 2]);
atlas_double = permute(double(atlas), [2, 3, 1]);
% Note: atlas dimensions are [dorsal-ventral, anterior-posterior, lateral]
% when read

% atlas = permute(atlas, [2, 3, 1]);
%%
% fUSmap = permute(PDI_allSF_avg, [3, 1, 2]);
% fUSmap_rs = permute(PDI_allSF_avg_rs, [3, 1, 2]); % resampled fUS map
fUSmap_rs = PDI_allSF_avg_rs; % resampled fUS map

%%
atlas_resized = imresize3(atlas_double, 1/5, 'Method', 'cubic'); % 10 um --> 50 um voxel size
% atlas_resized = imresize3(atlas_double, 1/2.5, 'Method', 'cubic'); % 10 um --> 25 um voxel size
% fUSmap_resized = imresize3(fUSmap, 2, 'Method', 'cubic'); % ~110 um --> ~50 um voxel size

%% Possible steps
% 1. Center the fUS volume (rotationally) and crop to a known size
% 2. Upsample/resize fUS volume to have isotropic voxel size equal to what
%    we use for the atlas
% 3. Manually ___ 

%% Temp testing for scatteredInterpolant
% Might have to do this for every time point???????????????
tempPDI = PDIallSF(:, :, :, 1);
tempRigidWarpedPDI = imwarp(tempPDI, fUSmap_50um_rigid_reg.tform);

UP = []; % Ultrasound points
AP = []; % Atlas points
for m = 1:size(controlPoints, 1) % Go through each set of control points (rows of the cell array)
    UP(m, :) = str2num(controlPoints{m, 2});
    AP(m, :) = str2num(controlPoints{m, 3});
end

UV = NaN(size(UP, 1), 1); % Ultrasound values (at each control point)
for m = 1:size(controlPoints, 1) % Go through each set of control points (rows of the cell array)
    UV(m) = tempPDI(UP(m, 1), UP(m, 2), UP(m, 3));
end

F = scatteredInterpolant(UP, UV, 'linear');


%% random
test_tform = fitgeotrans(UP, AP, "pwl"); % fitgeotrans(movingPoints, fixedPoints, "pwl")