%%
addpath('C:\Users\agzhou\Downloads\nrrdread')

%%
load('H:\Ultrasound data from 7-21-2025\08-28-2025 AZ01 RCA fUS air puff\run 1\fUS processing 09-02-2025 50 to 496\09-11-2025 processing  with 50-496 SVs and time lags 2-3 for CBFspeed with no extra g1 constraints\Trials 1 to 13 diff proc\activation_maps.mat')
% [atlas, atlas_metadata] = nrrdread('C:\Users\agzhou\Downloads\average_template_10.nrrd');
[atlas, atlas_metadata] = nrrdread('D:\Allen\Data\Allen Atlas CCFv3\Annotated 25\annotation_25.nrrd');

% Note: atlas dimensions are [dorsal-ventral, anterior-posterior, lateral]
% when read

% atlas = permute(atlas, [2, 3, 1]);
%%
fUSmap = permute(CBVi_allSF_avg, [3, 1, 2]);
%%
atlas_resized = imresize3(atlas, 1/5, 'Method', 'cubic'); % 10 um --> 50 um voxel size
fUSmap_resized = imresize3(fUSmap, 2, 'Method', 'cubic'); % ~110 um --> ~50 um voxel size

%% Possible steps
% 1. Center the fUS volume (rotationally) and crop to a known size
% 2. Upsample/resize fUS volume to have isotropic voxel size equal to what
%    we use for the atlas
% 3. Manually ___