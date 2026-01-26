%% Description
% - Apply deformable registration and do correlation analysis for functional
%   connectivity (in 3D space + time)
% - The ultrasound data should already be processed

% Inputs:
% - Power Doppler (PDI) template volume, usually an average across all
%   superframes in an experiment, at [50, 50, 50] um voxel size
% - PDI volumes at each superframe (PDIallSF)

% REQUIREMENTS:
% - Matlab 2024 and newer (for the registration, which uses the medical registration package)
% - Clone the allenCCF repository: https://github.com/cortex-lab/allenCCF.git
% - Clone the npy-matlab repository: https://github.com/kwikteam/npy-matlab.git
% - Add the Allen Atlas CCFv3 .nrrd files (average and annotated atlases) to a folder

%% Add paths for the preprocessed data (and any processed data to re-load)

data_dirpath = uigetdir('H:\Ultrasound data from 7-21-2025\10-21-2025 AZ01 FC RCA\Test proc 10-23-25\FC', 'Select the data directory');
data_dirpath = [data_dirpath, '\'];
% addpath(genpath(data_dirpath))

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
[US_template_FilePathFN, US_template_FilePath] = uigetfile([data_dirpath, 'PDI_template_50um.mat'], 'Select the ultrasound template (50 um)');
US_template_FilePath = [US_template_FilePath, US_template_FilePathFN];
load(US_template_FilePath)


%% (Optional) Open GUI for slightly "easier" manual registration
vs = prereg_params.target_voxel_size; % Voxel size in y, x, z
Rout = imref3d(size(AA_template_50um)); % Reference for the output of the transformation
% Rout = imref3d(size(AA_template_50um), vs(1), vs(2), vs(3)); % Reference for the output of the transformation

addpath([cd, '\..\'])
US_Atlas_Reg

%% (Optional) If loading already-registered data, load here
[saved_reg_data_FilePathFN, saved_reg_data_FilePath] = uigetfile([data_dirpath, 'fUSmap_50um.mat'], 'Select the registered ultrasound data (50 um)');
saved_reg_data_FilePath = [saved_reg_data_FilePath, saved_reg_data_FilePathFN];
load(saved_reg_data_FilePath)

%% Store the manual registration transformation info and 
rigid_tform_50um = fUSmap_50um_rigid_reg.tform;

% % **** manually define the rigid tform ****
% rigidRegParameterPrompt = {'x angle [deg]', 'y angle [deg]', 'z angle [deg]', 'x translation [voxels]', 'y translation [voxels]', 'z translation [voxels]'};
% rigidRegParameterDefaults = {'', '', '', '', '', ''};
% rigidRegParameterUserInput = inputdlg(rigidRegParameterPrompt, 'Input Parameters', 1, rigidRegParameterDefaults);
% 
% % define # of files manually for now
% % str2double(parameterUserInput{});
% rr.angles(2) = str2double(rigidRegParameterUserInput{1});
% rr.angles(1) = str2double(rigidRegParameterUserInput{2});
% rr.angles(3) = str2double(rigidRegParameterUserInput{3});
% rr.transl(2) = str2double(rigidRegParameterUserInput{4});
% rr.transl(1) = str2double(rigidRegParameterUserInput{5});
% rr.transl(3) = str2double(rigidRegParameterUserInput{6});
% 
% % rr.angles(1) = str2double(rigidRegParameterUserInput{1});
% % rr.angles(2) = str2double(rigidRegParameterUserInput{2});
% % rr.angles(3) = str2double(rigidRegParameterUserInput{3});
% % rr.transl(1) = str2double(rigidRegParameterUserInput{4});
% % rr.transl(2) = str2double(rigidRegParameterUserInput{5});
% % rr.transl(3) = str2double(rigidRegParameterUserInput{6});
% 
% rigid_tform_50um = rigidtform3d(rr.angles, rr.transl);
% 
% Rtest = imref3d(size(PDI_allSF_avg_rs))

% test = single(imwarp(PDI_allSF_avg_rs, rigid_tform_50um, 'cubic', 'OutputView', Rout));
% volumeViewer(test)

% Not sure if we need this: redefine the tform for 10 um voxel size
% rigid_tform_10um = fUSmap_50um_rigid_reg.tform;
% rigid_tform_10um.A(1:3, 4) = rigid_tform_50um.A(1:3, 4) .* (50/10);

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

[regions_FilePathFN, regions_FilePath] = uigetfile('C:\Users\Allen\Documents\GitHub\BU-Code\Allen code\Processing\Allen_Atlas_CCFv3_regions_of_interest.csv', 'Select the .csv for the brain regions of interest (Allen Atlas CCFv3)');
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

%% Create the ROI masks (at 50 um voxel size)
% region_masks_10um = cell(num_regions, 1);
region_masks_50um = cell(num_regions, 1);
for rn = 1:num_regions % region number
    region_mask_10um_temp = ( AA >= region_inds{rn, 1} & AA <= region_inds{rn, 2} );
    % region_masks_10um{rn} = region_mask_10um_temp;

    % Resize the 10um mask to 50um voxel size
    region_masks_50um{rn} = imresize3(region_mask_10um_temp, 10/50, 'Method', 'cubic');
end

clearvars region_mask_10um_temp rn

%% Separate the ROI masks by hemisphere (this all only works if the hemisphere axis is #2 in the data)
region_masks_50um_hemis = cell(num_regions, 2); % second dimension corresponds to 1 = left, 2 = right

% Define the indices for separating the hemispheres
hemis_dim = 2; % Dimension of the left-right axis
if mod(size(region_masks_50um{1}, hemis_dim), 2) ~= 0
    error("Code is not set up to deal with an odd # of voxels in the hemisphere axis direction")
end

% volumeViewer(squeeze(fUSmap_50um_rigid_reg.regVol.Voxels(:, 1:size(region_masks_50um{1}, hemis_dim)/2, :)))
% volumeViewer(squeeze(fUSmap_50um_rigid_reg.regVol.Voxels(:, size(region_masks_50um{1}, hemis_dim)/2 + 1:end, :)))
% 
% figure; imagesc(squeeze(max(fUSmap_50um_rigid_reg.regVol.Voxels(:, size(region_masks_50um{1}, hemis_dim)/2 + 1:end, :), [], 1)))

% Left should start from index = 1
hemis_inds.left = 1:size(region_masks_50um{1}, hemis_dim)/2;
hemis_inds.right = size(region_masks_50um{1}, hemis_dim)/2 + 1:size(region_masks_50um{1}, hemis_dim);

for rn = 1:num_regions % region number
    
    % Left
    % region_masks_50um_hemis{rn, 1} = region_masks_50um{rn}(:, hemis_inds.left, :);
    region_masks_50um_hemis{rn, 1} = region_masks_50um{rn};
    region_masks_50um_hemis{rn, 1}(:, hemis_inds.right, :) = false;

    % Right
    % region_masks_50um_hemis{rn, 2} = region_masks_50um{rn}(:, hemis_inds.right, :);
    region_masks_50um_hemis{rn, 2} = region_masks_50um{rn};
    region_masks_50um_hemis{rn, 2}(:, hemis_inds.left, :) = false;
end

% Add ROI name and acronyms for the hemisphere separation
region_names_hemis = cell(size(region_names, 1), 2);
region_acronyms_hemis = cell(size(region_acronyms, 1), 2);
for rn = 1:num_regions % region number
    region_name_temp = region_names{rn};
    region_acronym_temp = region_acronyms{rn};

    region_names_hemis{rn, 1} = region_name_temp + " - left";
    region_names_hemis{rn, 2} = region_name_temp + " - right";
    region_acronyms_hemis{rn, 1} = region_acronym_temp + " (L)";
    region_acronyms_hemis{rn, 2} = region_acronym_temp + " (R)";
end
clearvars region_mask_10um_temp rn region_name_temp region_acronym_temp

%% Add ROI info to a struct for saving
roi.names = region_names;
roi.acronyms = region_acronyms;
roi.inds = region_inds;
roi.num_regions = num_regions;
roi.masks_50um = region_masks_50um;
roi.names_hemis = region_names_hemis;
roi.acronyms_hemis = region_acronyms_hemis;
roi.acronyms_hemis_interleaved = roi.acronyms_hemis'; roi.acronyms_hemis_interleaved = roi.acronyms_hemis_interleaved(:);
roi.masks_50um_hemis = region_masks_50um_hemis;


% Save ROI info as a mat file
roi_savepath = uigetdir('J:\', 'Select the ROI info save path');
roi_savepath = [roi_savepath, '\'];

save([roi_savepath, 'roi_info.mat'], 'roi')

%% (Optional) Overlay the ROI masks onto a PDI template
% compareUStoAtlasROIs(PDI_template_reg_50um, region_masks_50um)
% compareUStoAtlasROIs(fUSmap_50um_rigid_reg.regVol.Voxels, region_masks_50um)
compareUStoAtlasROIs(fUSmap_50um_rigid_reg.regVol.Voxels, roi.masks_50um) % Non-hemisphere separated ROIs
compareUStoAtlasROIs(fUSmap_50um_rigid_reg.regVol.Voxels, roi.masks_50um_hemis(:)) % Hemisphere separated ROIs

%% Load the timing data (output of plotfUStiming_FC.m) and convert to actual time
[timingFilePathFN, timingFilePath] = uigetfile([data_dirpath, '..\..\Timing data\TD.mat'], 'Select the timing data');
timingFilePath = [timingFilePath, timingFilePathFN];
load(timingFilePath)

% Rename the actual time
t = TD.sfTimeTags;


% %% Upsample/interpolate over time to match the timing data...
% pupil_fr = 10; % Pupil data (behavioral camera) frame rate

%% Load the PDI across superframes data
[PDIallSF_FilePathFN, PDIallSF_FilePath] = uigetfile([data_dirpath, 'PDIallSF.mat'], 'Select the PDI across superframes data');
PDIallSF_FilePath = [PDIallSF_FilePath, PDIallSF_FilePathFN];
load(PDIallSF_FilePath)

%% Apply the registration/warping to the ultrasound data

% temp_US_template = imwarp(PDI_allSF_avg_rs, rigid_tform_50um, 'cubic', 'OutputView', Rout);

tic
% Go through every superframe and 1) resize to the 50 um voxel size, and 2) apply the registration
for sfi = 1:size(PDIallSF, 4) % superframe index
% for sfi = 1:2
    disp(sfi)
    PDI_sfi = squeeze(PDIallSF(:, :, :, sfi)); % PDI volume at superframe # sfi
    PDI_sfi_rs = imresize3(PDI_sfi, 'Scale', prereg_params.prereg_interp_factor, 'Method', 'cubic'); % Resize/resample

    if sfi == 1
        PDIallSF_reg = zeros([size(AA_template_50um), size(PDIallSF, 4)], 'single'); % Initialize the registered PDI across time
    end
    % Apply the transformation
    PDIallSF_reg(:, :, :, sfi) = single(imwarp(PDI_sfi_rs, rigid_tform_50um, 'cubic', 'OutputView', Rout));
end
toc

clearvars PDI_sfi sfi PDI_sfi_rs

PDIallSF_reg_avg = mean(PDIallSF_reg, 4);

%% (Optional) Save the registered PDI across superframes data, along with the ROI masks
% save([data_dirpath, 'PDIallSF_reg_ROI_masks_50um.mat'], "PDIallSF_reg", "region_masks_50um", '-v7.3')
save([data_dirpath, 'PDIallSF_reg_ROI_info_50um.mat'], "PDIallSF_reg", "roi", '-v7.3')

%% (Optional) Load the registered PDI across superframes data, along with the ROI masks
[PDIallSF_reg_ROI_masks_FilePathFN, PDIallSF_reg_ROI_masks_FilePath] = uigetfile([data_dirpath, 'PDIallSF_reg_ROI_info_50um.mat'], 'Select the registered PDI across superframes data');
PDIallSF_reg_ROI_masks_FilePath = [PDIallSF_reg_ROI_masks_FilePath, PDIallSF_reg_ROI_masks_FilePathFN];
load(PDIallSF_reg_ROI_masks_FilePath)

%% Global mean subtraction
PDI_reg_global_mean = squeeze(mean(PDIallSF_reg, [1, 2, 3])); % Global mean timecourse (voxelwise) of the registered PDI
figure; plot(t, PDI_reg_global_mean); title("Global mean (voxelwise)"); xlabel("Time [s]"); ylabel("PDI")

%% Get ROI averaged PDI timecourses
% (old incorrect title: Correlation without resampling PDI in time (with sliding window))
% figure; plot(sfTimeTags)
% figure; plot(diff(sfTimeTags))
num_sf = size(PDIallSF_reg, 4); % # of superframes

% PDI_ROI_timecourses = cell(num_regions, num_sf); % Store average ROI PDI timecourses in a cell array (each cell is an average timecourse)
PDI_ROI_timecourses = cell(num_regions, 1); % Store average ROI PDI timecourses in a cell array (each cell is an average timecourse)
PDI_ROI_hemis_timecourses = cell(num_regions, 2); % Store average ROI PDI (hemisphere-separated) timecourses in a cell array (each cell is an average timecourse)

% -- Calculate PDI [ROI average] timecourses -- %
tic
for ti = 1:num_sf % "time" index -- go through each superframe
% for ti = 11
    disp(ti)
    PDIallSF_reg_ti_temp = squeeze(PDIallSF_reg(:, :, :, ti)); % Registered volume at "time" index ti

    for ri = 1:num_regions % region/ROI index -- loop through each region

        % Normal
        ROI_mask_temp = region_masks_50um{ri}; % ROI #ri mask
        PDI_ri_masked_temp = PDIallSF_reg_ti_temp(ROI_mask_temp); % Vectorized voxels of the registered PDI at "time" index ti
        PDI_ROI_timecourses{ri}(ti) = mean(PDI_ri_masked_temp);
        
        % Hemisphere-separated
        ROI_mask_temp_left = region_masks_50um_hemis{ri, 1}; % ROI #ri mask (left)
        ROI_mask_temp_right = region_masks_50um_hemis{ri, 2}; % ROI #ri mask (right)
        PDI_ri_masked_temp_left = PDIallSF_reg_ti_temp(ROI_mask_temp_left); % Vectorized voxels of the registered PDI at "time" index ti
        PDI_ri_masked_temp_right = PDIallSF_reg_ti_temp(ROI_mask_temp_right); % Vectorized voxels of the registered PDI at "time" index ti
        PDI_ROI_hemis_timecourses{ri, 1}(ti) = mean(PDI_ri_masked_temp_left);
        PDI_ROI_hemis_timecourses{ri, 2}(ti) = mean(PDI_ri_masked_temp_right);

    end
end
% clearvars ti ri PDIallSF_reg_ti_temp ROI_mask_temp PDI_ri_masked_temp
toc

% Make any row timecourses into column vectors
for ri = 1:num_regions
    PDI_ROI_timecourses{ri} = squeeze(PDI_ROI_timecourses{ri}');

    PDI_ROI_hemis_timecourses{ri, 1} = squeeze(PDI_ROI_hemis_timecourses{ri, 1}');
    PDI_ROI_hemis_timecourses{ri, 2} = squeeze(PDI_ROI_hemis_timecourses{ri, 2}');

end

%% Add global mean subtracted versions
PDI_ROI_GMS_timecourses = cell(num_regions, 1);
PDI_ROI_hemis_GMS_timecourses = cell(num_regions, 2);
for ri = 1:num_regions
    PDI_ROI_GMS_timecourses{ri} = PDI_ROI_timecourses{ri} - PDI_reg_global_mean;
    PDI_ROI_hemis_GMS_timecourses{ri, 1} = PDI_ROI_hemis_timecourses{ri, 1} - PDI_reg_global_mean;
    PDI_ROI_hemis_GMS_timecourses{ri, 2} = PDI_ROI_hemis_timecourses{ri, 2} - PDI_reg_global_mean;
end

%% Store the ROI timecourses in matrix form
PDI_ROI_timecourses_mat = zeros(length(t), num_regions); % Still ROI-averaged PDI timecourses, but in matrix form (each column is a separate ROI timecourse). Dimensions: [# time points, # ROIs]
for ri = 1:num_regions % region/ROI index -- loop through each region
    PDI_ROI_timecourses_mat(:, ri) = PDI_ROI_timecourses{ri};
end
% figure; plot(PDI_ROI_timecourses{1})

% Hemisphere-separated version
PDI_ROI_hemis_timecourses_mat = zeros(length(t), num_regions*2); % Still ROI-averaged PDI timecourses, but in matrix form (each column is a separate ROI timecourse). Dimensions: [# time points, # ROIs]
for ri = 1:num_regions % region/ROI index -- loop through each region
    PDI_ROI_hemis_timecourses_mat(:, (ri - 1)*2 + 1) = PDI_ROI_hemis_timecourses{ri, 1};
    PDI_ROI_hemis_timecourses_mat(:, ri*2) = PDI_ROI_hemis_timecourses{ri, 2};
end
% figure; plot(PDI_ROI_timecourses{1})

%% Store the ROI timecourses in matrix form (GMS version)
PDI_ROI_GMS_timecourses_mat = zeros(length(t), num_regions); % Still ROI-averaged PDI timecourses, but in matrix form (each column is a separate ROI timecourse). Dimensions: [# time points, # ROIs]
for ri = 1:num_regions % region/ROI index -- loop through each region
    PDI_ROI_GMS_timecourses_mat(:, ri) = PDI_ROI_GMS_timecourses{ri};
end

% Hemisphere-separated version
PDI_ROI_hemis_GMS_timecourses_mat = zeros(length(t), num_regions*2); % Still ROI-averaged PDI timecourses, but in matrix form (each column is a separate ROI timecourse). Dimensions: [# time points, # ROIs]
for ri = 1:num_regions % region/ROI index -- loop through each region
    PDI_ROI_hemis_GMS_timecourses_mat(:, (ri - 1)*2 + 1) = PDI_ROI_hemis_GMS_timecourses{ri, 1};
    PDI_ROI_hemis_GMS_timecourses_mat(:, ri*2) = PDI_ROI_hemis_GMS_timecourses{ri, 2};
end
% figure; plot(PDI_ROI_timecourses{1})

%% Calculate the Global Variance of the Temporal Derivative (GVTD): data-driven motion quantification
% diff_PDIallSF_reg = diff(PDIallSF_reg, 1, length(size(PDIallSF_reg))); % 1st order diff of the registered PDI across superframes, across time (the last dimension)
size_PDIallSF_reg = size(PDIallSF_reg);
numVoxelsInVolume = size_PDIallSF_reg(1) * size_PDIallSF_reg(2) * size_PDIallSF_reg(3); % # of voxels in each volume
GVTD = squeeze( sum( diff(PDIallSF_reg, 1, length(size(PDIallSF_reg))) .^ 2, [1, 2, 3] ) ./ numVoxelsInVolume ) .^ 0.5;
GVTD(end + 1) = NaN; % pad the end with a NaN, since there is no forward point past the last time point

% figure; plot(t, GVTD); title("Global Variance of the Temporal Derivative (GVTD) - PDIallSF"); xlabel("Time [s]"); ylabel("GVTD")

figure
% yyaxis left

% Add accelerometer
% yyaxis right
accel = sqrt(sum(TD.inScanData.^2, 2));
accel_zm = accel - mean(accel);

plot(TD.daqTimeTags, accel_zm ./ max(accel_zm), '--')
hold on
plot(t, GVTD ./ max(GVTD), 'LineWidth', 2); title("Global Variance of the Temporal Derivative (GVTD) of PDIallSF vs. Accelerometer"); xlabel("Time [s]"); ylabel("GVTD")
plot(t, PDI_reg_global_mean ./ max(PDI_reg_global_mean))
ylabel("Accelerometer amplitude")
% legend("GVTD", "Accelerometer component 1", "Accelerometer component 2", "Accelerometer component 3")
legend("Accelerometer component 1", "GVTD", "PDI global mean")

%% Plot each ROI's PDI timecourse
% subplot(num_regions, 1, 1)
% figure;
% ROI_PDI_timecourse_tl = tiledlayout("vertical"); % Vertical tile layout
% for ri = 1:num_regions % region/ROI index -- loop through each region
%     nexttile
%     % subplot(num_regions, 1, ri)
%     plot(t, PDI_ROI_timecourses{ri})
% end
% ROI_PDI_timecourse_tl.TileSpacing = 'compact';
% ROI_PDI_timecourse_tl.Padding = 'compact';
% title(ROI_PDI_timecourse_tl, "ROI average PDI timecourses")
% xlabel(ROI_PDI_timecourse_tl, "Time [s]")
% ylabel(ROI_PDI_timecourse_tl, "PDI magnitude [au]")

% Plot multiple stackedplots to visualize ROI PDI timecourses
num_cols_per_sp = 16;
num_sps = ceil(num_regions/num_cols_per_sp); % # of stackedplot to use since they only allow 25 columns max

for spi = 1:num_sps
    % - NOTE: stackedplot only allows for 25 columns max - % 
    % Normal
    figure
    temp_ind_spi = (spi - 1)*num_cols_per_sp + 1:spi*num_cols_per_sp;
    % % ROI_PDI_timecourse_sp = stackedplot(t, PDI_ROI_timecourses_mat, 'DisplayLabels', region_acronyms);
    % ROI_PDI_timecourse_sp = stackedplot(t, [PDI_ROI_timecourses_mat(:, temp_ind_spi), GVTD], 'DisplayLabels', [roi.acronyms(temp_ind_spi); {'GVTD'}]);
    ROI_PDI_timecourse_sp = stackedplot(t, [PDI_ROI_GMS_timecourses_mat(:, temp_ind_spi), GVTD], 'DisplayLabels', [roi.acronyms(temp_ind_spi); {'GVTD'}]);
    title("ROI average PDI timecourses")
    xlabel("Time [s]")
    % fontsize(14, 'points')

end

for spi = 1:num_sps*2
    % - NOTE: stackedplot only allows for 25 columns max - % 
    
    % Hemisphere-separated
    figure
    temp_ind_spi = (spi - 1)*num_cols_per_sp + 1:spi*num_cols_per_sp;

    % ROI_PDI_timecourse_sp = stackedplot(t, [PDI_ROI_hemis_timecourses_mat(:, temp_ind_spi), GVTD], 'DisplayLabels', [roi.acronyms_hemis_interleaved(temp_ind_spi); {'GVTD'}]);
    ROI_PDI_timecourse_sp = stackedplot(t, [PDI_ROI_hemis_GMS_timecourses_mat(:, temp_ind_spi), GVTD], 'DisplayLabels', [roi.acronyms_hemis_interleaved(temp_ind_spi); {'GVTD'}]);
    title("ROI average (hemisphere-separated) PDI timecourses")
    xlabel("Time [s]")
    % fontsize(14, 'points')
end

% % - NOTE: stackedplot only allows for 25 columns max - % 
% % Normal
% figure
% % ROI_PDI_timecourse_sp = stackedplot(t, PDI_ROI_timecourses_mat, 'DisplayLabels', region_acronyms);
% ROI_PDI_timecourse_sp = stackedplot(t, [PDI_ROI_timecourses_mat, GVTD], 'DisplayLabels', [roi.acronyms; {'GVTD'}]);
% title("ROI average PDI timecourses")
% xlabel("Time [s]")
% % fontsize(14, 'points')
% 
% % Hemisphere-separated
% figure
% ROI_PDI_hemis_timecourse_sp = stackedplot(t, PDI_ROI_hemis_timecourses_mat, 'DisplayLabels', roi.acronyms_hemis_interleaved); 
% % ROI_PDI_timecourse_sp = stackedplot(t, [PDI_ROI_hemis_timecourses_mat, GVTD], 'DisplayLabels', [region_acronyms; {'GVTD'}]);
% title("ROI average (hemisphere-separated) PDI timecourses")
% xlabel("Time [s]")
% % fontsize(14, 'points')

%% FC correlation matrix (full timecourse)
% corr_PDI = corrcoef(PDI_ROI_timecourses_mat);
% corr_PDI_hemis = corrcoef(PDI_ROI_hemis_timecourses_mat);

corr_PDI_GMS = corrcoef(PDI_ROI_GMS_timecourses_mat);
corr_PDI_hemis_GMS = corrcoef(PDI_ROI_hemis_GMS_timecourses_mat);

%% Plot the FC correlation matrices (full timecourse)

plotCM(corr_PDI, roi)
plotCM(corr_PDI_hemis, roi, true)

plotCM(corr_PDI_GMS, roi)
plotCM(corr_PDI_hemis_GMS, roi, true)

%% -- Calculate changes in FC (seed correlation matrices) over time with sliding windows -- %

corr_ws = 30; % Correlation sliding window size (I need to convert this to be specific with actual time)
% For now, 30 represents around 30 sf * ~0.5 seconds per sf --> ~ 15 seconds

corr_sw_PDI = zeros(num_regions, num_regions, num_sf); % Sliding window PDI seed correlation matrices 
% THIS VERSION WORKS BASED ON THE SUPERFRAME INDICES, NOT TIME DIRECTLY
tic
% for wi = 1:num_sf - corr_ws
for wi = 1:num_sf
% for wi = 1:2
% for wi = 700

    sfis = wi:( wi + corr_ws - 1 ); % Superframe indices in the window

    corr_data_matrix_temp = zeros(corr_ws, num_regions); % Temporary data matrix for calculating the correlation matrix. Each column is the PDI timecourse for a region, within the time window defined by 'sfis'

    % Go through each region and temporarily store the timecourse in the
    % corresponding column of the data matrix
    for ri = 1:num_regions
        
        if wi - num_sf - corr_ws < 0 % Mirror/pad the data at the end of the sliding window range
            temp_PDI_ROI_timecourse_ri = padarray(PDI_ROI_timecourses{ri}', corr_ws, 'symmetric', 'post');
            % temp_PDI_ROI_timecourse_ri = paddata(PDI_ROI_timecourses{ri}, corr_ws, Dimension = 2, Pattern = 'reflect', Side = 'trailing');
            corr_data_matrix_temp(:, ri) = ( temp_PDI_ROI_timecourse_ri(sfis) );
        else % Otherwise, no need for padding
            corr_data_matrix_temp(:, ri) = ( PDI_ROI_timecourses{ri}(sfis) )';
        end
        
    end

    corr_sw_PDI(:, :, wi) = corrcoef(corr_data_matrix_temp); % Calculate the seed correlation matrix at that window

end
toc
clearvars wi ri sfis temp_PDI_ROI_timecourse_ri corr_data_matrix_temp

%% Plot the seed correlation matrices

% Plot the matrix for one window
figure; imagesc(squeeze(corr_sw_PDI(:, :, 100))); colormap jet; axis square; colorbar; clim([-1, 1])
% xticks(1:num_regions); xticklabels(region_names); yticks(1:num_regions); yticklabels(region_names) % Set the tick labels to be the ROI names
xticks(1:num_regions); xticklabels(region_acronyms); yticks(1:num_regions); yticklabels(region_acronyms) % Set the tick labels to be the ROI acronyms

%% Spaghetti plot of the correlation between each pair of regions, over time
corr_sw_legend = {}; % Legend for each pair
figure; hold on
for m = 1:num_regions
    for n = m + 1:num_regions
        plot(t, squeeze(corr_sw_PDI(m, n, :)))
        corr_sw_legend(end + 1) = {region_acronyms{m} + "-" + region_acronyms{n}};
    end
end
hold off
% xlabel('sf index')
xlabel('Time [s]')
ylabel('Correlation coefficient')
title("Correlation between ROIs, with sliding window size = " + num2str(corr_ws) + " superframes")
legend(corr_sw_legend)

%% Spaghetti plot of the correlation between each pair of regions, over time
% **** with the GVTD plotted below ****
corr_sw_legend = {}; % Legend for each pair
figure
tiledlayout('vertical')
nexttile; hold on
for m = 1:num_regions
    for n = m + 1:num_regions
        plot(t, squeeze(corr_sw_PDI(m, n, :)))
        corr_sw_legend(end + 1) = {region_acronyms{m} + "-" + region_acronyms{n}};
    end
end
hold off
% xlabel('sf index')
xlabel('Time [s]')
ylabel('Correlation coefficient')
title("Correlation between ROIs, with sliding window size = " + num2str(corr_ws) + " superframes")
legend(corr_sw_legend)

nexttile
plot(t, GVTD); xlabel("Time [s]"); ylabel("GVTD [au]")

%% Spaghetti plot of the correlation between each pair of regions, over time
% **** Testing: plot the % change in correlation coefficient compared to
%               its mean or median
corr_sw_legend = {}; % Legend for each pair
figure; hold on
for m = 1:num_regions
    for n = m + 1:num_regions
        % plot(t, (squeeze(corr_sw_PDI(m, n, :)) - mean(squeeze(corr_sw_PDI(m, n, :))) ./ squeeze(corr_sw_PDI(m, n, :))) .* 100)
        plot(t, (squeeze(corr_sw_PDI(m, n, :)) - median(squeeze(corr_sw_PDI(m, n, :))) ./ squeeze(corr_sw_PDI(m, n, :))) .* 100)
        corr_sw_legend(end + 1) = {region_acronyms{m} + "-" + region_acronyms{n}};
    end
end
hold off
ylim([-1, 1])
% xlabel('sf index')
xlabel('Time [s]')
ylabel('% change')
% title("Percent change in correlation coefficient (vs. temporal mean) between ROIs, with sliding window size = " + num2str(corr_ws) + " superframes")
title("Percent change in correlation coefficient (vs. temporal median) between ROIs, with sliding window size = " + num2str(corr_ws) + " superframes")
legend(corr_sw_legend)

%% Spaghetti plot of the correlation between each pair of regions, over time
% **** Testing: plot the relative value of the correlation coefficient compared to
%               its mean or median
corr_sw_legend = {}; % Legend for each pair
figure;
tiledlayout('vertical')
nexttile
hold on
for m = 1:num_regions
    for n = m + 1:num_regions
        % plot(t, squeeze(corr_sw_PDI(m, n, :)) ./ mean(squeeze(corr_sw_PDI(m, n, :))) .* 100)
        plot(t, squeeze(corr_sw_PDI(m, n, :)) ./ median(squeeze(corr_sw_PDI(m, n, :))) .* 100)
        corr_sw_legend(end + 1) = {region_acronyms{m} + "-" + region_acronyms{n}};
    end
end
hold off
% xlabel('sf index')
xlabel('Time [s]')
ylabel('Relative correlation coefficient')
% title("Relative correlation coefficient (vs. temporal mean) between ROIs, with sliding window size = " + num2str(corr_ws) + " superframes")
title("Relative correlation coefficient (vs. temporal median) between ROIs, with sliding window size = " + num2str(corr_ws) + " superframes")
legend(corr_sw_legend)

nexttile
plot(t, GVTD); xlabel("Time [s]"); ylabel("GVTD [au]")

%% Spaghetti plot of the correlation between each pair of regions, over time
% **** Testing: plot the normalized relative change in correlation coefficient compared to its mean
corr_sw_legend = {}; % Legend for each pair
figure; hold on
for m = 1:num_regions
    for n = m + 1:num_regions
        % temp = squeeze(corr_sw_PDI(m, n, :)) ./ mean(squeeze(corr_sw_PDI(m, n, :)));
        temp = (squeeze(corr_sw_PDI(m, n, :)) - mean(squeeze(corr_sw_PDI(m, n, :))) ./ squeeze(corr_sw_PDI(m, n, :)));
        temp = temp ./ max(abs(temp));
        plot(t, temp)
        corr_sw_legend(end + 1) = {region_acronyms{m} + "-" + region_acronyms{n}};
    end
end
clearvars temp
hold off
% xlabel('sf index')
xlabel('Time [s]')
ylabel('% change')
title("Normalized relative change in correlation coefficient (vs. temporal mean) between ROIs, with sliding window size = " + num2str(corr_ws) + " superframes")
legend(corr_sw_legend)

%% Spaghetti plot of the diff of the correlation between each pair of regions, over time
% **** Testing
corr_sw_legend = {}; % Legend for each pair
figure;
% tiledlayout('vertical')
% nexttile
hold on
for m = 1:num_regions
    for n = m + 1:num_regions
        % plot(t, squeeze(corr_sw_PDI(m, n, :)) ./ mean(squeeze(corr_sw_PDI(m, n, :))) .* 100)
        plot(t(1:end - 1), diff(squeeze(corr_sw_PDI(m, n, :))) )
        corr_sw_legend(end + 1) = {region_acronyms{m} + "-" + region_acronyms{n}};
    end
end
hold off
% xlabel('sf index')
xlabel('Time [s]')
ylabel('Temporal difference of correlation coefficient')
title("Temporal difference of correlation coefficient between ROIs, with sliding window size = " + num2str(corr_ws) + " superframes")
legend(corr_sw_legend)

%% Plot all the ROI average PDI timecourses
figure; hold on; xlabel('Time [s]'); ylabel('PDI ROI average')
for ind = 1:length(PDI_ROI_timecourses)
    % plot(sfTimeTags, filtfilt(ones(1,4),4,PDI_ROI_timecourses{ind}))
    plot(sfTimeTags, PDI_ROI_timecourses{ind})
end
legend(region_acronyms)




%% Video of the registered PDI across superframes (MIPs)
addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Processing')
% figure; imagesc(squeeze(max(PDIallSF_reg(:, :, :, 1), [], 1))'); colormap hot; axis equal
reg_volume_size_mm = size(AA_template_50um) .* 50e-6 .* 1e3;
generateTiffStack_acrossframes(PDIallSF_reg .^ 1, reg_volume_size_mm, 'hot', 1:size(PDIallSF_reg, 1), t)

%% Load pupil tracking video
[pupilData] = readMP4; % Select the pupil video and read it

% Get and plot the time points of the pupil data relative to the start of the ultrasound acquisition
time_diff = startTimetag - pupilData.startTimetag; % Time difference between the ultrasound acquisition and the pupil data start
time_diff.Format = 'hh:mm:ss.SSS'; % Show millisecond precision
pupilData.timestamps_relative_to_US_start = pupilData.timestamps - seconds(time_diff);

figure;
plot(sfTimeTags, 1, 'x')
hold on
plot(pupilData.timestamps_relative_to_US_start, 1, 'o')
hold off
xlabel('Time [s]')

%% Define ROIs for the eye
addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Processing\Pupil tracking')
eyeROIs = defineEyeROIs(pupilData.frames(:, :, 1));


