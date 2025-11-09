%% Description
% Apply deformable registration and do correlation analysis for functional
% connectivity (in 3D space + time)

% REQUIREMENTS:
% - Clone the allenCCF repository: https://github.com/cortex-lab/allenCCF.git
% - Clone the npy-matlab repository: https://github.com/kwikteam/npy-matlab.git

%% Add paths for reading nrrd files + Load the template atlas

% Load the Allen Atlas CCFv3 path
AACCFv3dirpath = uigetdir('H:\Allen Atlas CCFv3', 'Select the Allen Atlas CCFv3 directory');
AACCFv3dirpath = [AACCFv3dirpath, '\'];
addpath(genpath(AACCFv3dirpath))
% addpath([AACCFv3dirpath, 'nrrdread']) % Add the nrrdread path

% % Load the template (10 um voxel size)
% [AA_template_FilePathFN, AA_template_FilePath] = uigetfile([AACCFv3dirpath, 'average_template_10.nrrd'], 'Select the average template atlas (10 um)');
% AA_template_FilePath = [AA_template_FilePath, AA_template_FilePathFN];
% [AA_template_10um_up, AA_template_10um_up_metadata] = nrrdread(AA_template_FilePath); % Unpermuted annotated atlas
% AA_template_10um = double(permute(AA_template_10um_up, [2, 3, 1])); % Permuted annotated atlas
% % Note: atlas dimensions are [dorsal-ventral, anterior-posterior, lateral] when read
% AA_template_50um = imresize3(AA_template_10um, 1/5, 'Method', 'cubic');

% Load the template (50 um voxel size)
[AA_template_FilePathFN, AA_template_FilePath] = uigetfile([AACCFv3dirpath, 'average_template_50.nrrd'], 'Select the average template atlas (50 um)');
AA_template_FilePath = [AA_template_FilePath, AA_template_FilePathFN];
[AA_template_50um_up, AA_template_50um_up_metadata] = nrrdread(AA_template_FilePath); % Unpermuted annotated atlas
AA_template_50um = double(permute(AA_template_50um_up, [2, 3, 1])); % Permuted annotated atlas
% Note: atlas dimensions are [dorsal-ventral, anterior-posterior, lateral] when read


%% Load the ultrasound map of interest
% Load fUS data (template map)
% [US_template_FilePathFN, US_template_FilePath] = uigetfile(['PDI_template_10um.mat'], 'Select the ultrasound template (10 um)');
% US_template_FilePath = [US_template_FilePath, US_template_FilePathFN];
% load(US_template_FilePath)
% ... resize to 50 um for ease of manual registration ...

% 50 um version
[US_template_FilePathFN, US_template_FilePath] = uigetfile(['PDI_template_50um.mat'], 'Select the ultrasound template (50 um)');
US_template_FilePath = [US_template_FilePath, US_template_FilePathFN];
load(US_template_FilePath)


%% (Optional) Open GUI for slightly "easier" manual registration
US_Atlas_Reg

%% Store the manual registration transformation info and 
% rigid_tform_50um = fUSmap_50um_rigid_reg.tform;
rigid_tform_50um = rigidtform3d([-16, 0, -6], [20, 40, 44.5]); % HARD CODING TEMPORARILY FOR 8/28/25 fUS DATA
% Not sure if we need this: redefine the tform for 10 um voxel size
% rigid_tform_10um = fUSmap_50um_rigid_reg.tform;
% rigid_tform_10um.A(1:3, 4) = rigid_tform_50um.A(1:3, 4) .* (50/10);

%% Apply the registration/warping to the ultrasound data !!!!!!!!!!!!

Rout = imref3d(size(AA_template_50um)); % Reference for the output of the transformation
% temp_US_template = imwarp(PDI_allSF_avg_rs, rigid_tform_50um, 'cubic', 'OutputView', Rout);

% fUS template maps
CBVi_allSF_avg_rs_reg = imwarp(CBVi_allSF_avg_rs, rigid_tform_50um, 'cubic', 'OutputView', Rout);
CBFsi_allSF_avg_rs_reg = imwarp(CBFsi_allSF_avg_rs, rigid_tform_50um, 'cubic', 'OutputView', Rout);
PDI_allSF_avg_rs_reg = imwarp(PDI_allSF_avg_rs, rigid_tform_50um, 'cubic', 'OutputView', Rout);

% fUS activation maps
am_rCBV_rs_reg = imwarp(am_rCBV_rs, rigid_tform_50um, 'cubic', 'OutputView', Rout);
am_rCBFspeed_rs_reg = imwarp(am_rCBFspeed_rs, rigid_tform_50um, 'cubic', 'OutputView', Rout);
am_rPDI_rs_reg = imwarp(am_rPDI_rs, rigid_tform_50um, 'cubic', 'OutputView', Rout);

%% Load the NPY path
NPYdirpath = uigetdir('C:\Users\Allen\Documents\GitHub\npy-matlab', 'Select the NPY-Matlab directory');
NPYdirpath = [NPYdirpath, '\'];

addpath(genpath(NPYdirpath)) % Add the NPY-Matlab path + all subfolders
% addpath([NPYdirpath, 'npy-matlab\']) % Add the NPY-Matlab path

%% Load the allenCCF path (use the Cortex Lab's functions)
allenCCFdirpath = uigetdir('C:\Users\Allen\Documents\GitHub\allenCCF', 'Select the allenCCF directory');
allenCCFdirpath = [allenCCFdirpath, '\'];

addpath(genpath(allenCCFdirpath))
% addpath([allenCCFdirpath, 'npy-matlab\']) % Add the NPY-Matlab path

%% Load the 2017 modified templates from the Cortex Lab
% Load the annotated volume
% Note: the number at each voxel corresponds to the region/area index:
% "The original volume has numbers that correspond to the "id" field in the structure tree, 
%  but since I wanted to make a colormap for these, I re-indexed the annotation volume by 
%  the row number of the structure tree. So in this version the values correspond to "index"+1. 
%  This also allows using uint16 datatype, cutting file size in half. See setup_utils.m."

AA_up = readNPY([AACCFv3dirpath, '2017 modified by Cortex Lab UCL\annotation_volume_10um_by_index.npy']); % Unpermuted annotated atlas (10 um voxel size)
AA = permute(AA_up, [1, 3, 2]); % Permuted annotated atlas (10 um voxel size)

% Load the structure tree for the annotated volume (should be in the allenCCF Github directory)
ST = loadStructureTree('structure_tree_safe_2017.csv'); % Structure tree: a table of what all the labels in the annotated volume mean

%% (Checking) view the annotation
% volumeViewer(AA)

figure; imagesc(squeeze(AA(600, :, :))')

%% Define which regions to look at
% use the "index + 1" part of the table and combine subregions inf needed
region_indices = {};


% % MO_ind = [13:18]; % Somatomotor areas (MO)
% MOp_ind = [19:24]; % Primary motor area (MOp)
% MOs_ind = [25:30]; % Secondary motor area (MOs)
% % SSp_ind = [38:44]; % Primary somatosensory area (SSp)
% SSPn_ind = [45:51]; % Primary somatosensory area, nose (SSp-n)
% 
% VISp_ind = [186:192]; % Primary visual area (VISp)

% ==== Regions are defined in a spreadsheet: auto read them ==== %

[regions_FilePathFN, regions_FilePath] = uigetfile('C:\Users\Allen\Documents\GitHub\BU-Code\Allen code\Processing\Allen_Atlas_CCFv3_regions_of_interest_whiskerstim.csv', 'Select the .csv for the brain regions of interest (Allen Atlas CCFv3)');
regions_FilePath = [regions_FilePath, regions_FilePathFN];

% Read the .csv
regions_opts = detectImportOptions(regions_FilePath);
regionsTable = readtable(regions_FilePath, regions_opts);
regionsCell = table2cell(regionsTable);

% Store the regions of interest's names, abbreviations, and start/end
% "index + 1"s into cell arrays
region_names = regionsCell(:, 1);
region_acronyms = regionsCell(:, 2);
region_inds = regionsCell(:, 3:4);
num_regions = size(regionsCell, 1); % # of regions of interest

% region_masks_10um = cell(num_regions, 1);
region_masks_50um = cell(num_regions, 1);
for rn = 1:num_regions % region number
    region_mask_10um_temp = ( AA >= region_inds{rn, 1} & AA <= region_inds{rn, 2} );
    % region_masks_10um{rn} = region_mask_10um_temp;

    % Resize the 10um mask to 50um voxel size
    region_masks_50um{rn} = imresize3(region_mask_10um_temp, 10/50, 'Method', 'cubic');
end

clearvars region_mask_10um_temp rn


