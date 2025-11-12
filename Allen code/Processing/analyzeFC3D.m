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
rigid_tform_50um = fUSmap_50um_rigid_reg.tform;

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

% region_masks_10um = cell(num_regions, 1);
region_masks_50um = cell(num_regions, 1);
for rn = 1:num_regions % region number
    region_mask_10um_temp = ( AA >= region_inds{rn, 1} & AA <= region_inds{rn, 2} );
    % region_masks_10um{rn} = region_mask_10um_temp;

    % Resize the 10um mask to 50um voxel size
    region_masks_50um{rn} = imresize3(region_mask_10um_temp, 10/50, 'Method', 'cubic');
end

clearvars region_mask_10um_temp rn

%% Load the timing data
[timingFilePathFN, timingFilePath] = uigetfile(['..\Timing data\TD.mat'], 'Select the timing data');
timingFilePath = [timingFilePath, timingFilePathFN];
load(timingFilePath)

% %% Upsample/interpolate over time to match the timing data...
% pupil_fr = 10; % Pupil data (behavioral camera) frame rate

%% Apply the registration/warping to the ultrasound data

Rout = imref3d(size(AA_template_50um)); % Reference for the output of the transformation
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

%% Correlation without resampling PDI in time
% figure; plot(sfTimeTags)
% figure; plot(diff(sfTimeTags))
num_sf = size(PDIallSF_reg, 4); % # of superframes

corr_ws = 30; % Correlation sliding window size (I need to convert this to be specific with actual time)
% For now, 30 represents around 30 sf * ~0.5 seconds per sf --> ~ 15 seconds

% PDI_ROI_timecourses = cell(num_regions, num_sf); % Store average ROI PDI timecourses in a cell array (each cell is an average timecourse)
PDI_ROI_timecourses = cell(num_regions, 1); % Store average ROI PDI timecourses in a cell array (each cell is an average timecourse)

% -- Calculate PDI [ROI average] timecourses -- %
tic
for ti = 1:num_sf % "time" index -- go through each superframe
% for ti = 11
    PDIallSF_reg_ti_temp = squeeze(PDIallSF_reg(:, :, :, ti)); % Registered volume at "time" index ti

    for ri = 1:num_regions % region/ROI index -- loop through each region

        ROI_mask_temp = region_masks_50um{ri}; % ROI #ri mask
        PDI_ri_masked_temp = PDIallSF_reg_ti_temp(ROI_mask_temp); % Vectorized voxels of the registered PDI at "time" index ti
        PDI_ROI_timecourses{ri}(ti) = mean(PDI_ri_masked_temp);

    end
end
clearvars ti ri PDIallSF_reg_ti_temp ROI_mask_temp PDI_ri_masked_temp
toc

% figure; plot(PDI_ROI_timecourses{1})

%% -- Calculate changes in FC (seed correlation matrices) over time with sliding windows -- %
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

%% Plot the seed correlation matrices

% Plot the matrix for one window
figure; imagesc(squeeze(corr_sw_PDI(:, :, 100))); colormap hot; axis square; colorbar
% xticks(1:num_regions); xticklabels(region_names); yticks(1:num_regions); yticklabels(region_names) % Set the tick labels to be the ROI names
xticks(1:num_regions); xticklabels(region_acronyms); yticks(1:num_regions); yticklabels(region_acronyms) % Set the tick labels to be the ROI acronyms

% Spaghetti plot of the correlation between each pair of regions, over time
corr_sw_legend = {}; % Legend for each pair
figure; hold on
for m = 1:num_regions
    for n = m + 1:num_regions
        plot(sfTimeTags, squeeze(corr_sw_PDI(m, n, :)))
        corr_sw_legend(end + 1) = {region_acronyms{m} + "-" + region_acronyms{n}};
    end
end
hold off
% xlabel('sf index')
xlabel('Time [s]')
ylabel('Correlation coefficient')
title("Correlation between ROIs, with sliding window size = " + num2str(corr_ws) + " superframes")
legend(corr_sw_legend)

%% Plot all the ROI average PDI timecourses
figure; hold on; xlabel('Time [s]'); ylabel('PDI ROI average')
for ind = 1:length(PDI_ROI_timecourses)
    % plot(sfTimeTags, filtfilt(ones(1,4),4,PDI_ROI_timecourses{ind}))
    plot(sfTimeTags, PDI_ROI_timecourses{ind})
end
legend(region_acronyms)





%% Load pupil tracking video
[pupilData] = readMP4;

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




%% Separate each trial
ah = 3; % Approximate a cutoff value for analog high

ind_above_ah = find(TD.airPuffOutput > ah); % Get indices of the air puff output above analog high
ind_shift_below_ah = find(TD.airPuffOutput(ind_above_ah - 1) < ah); % See which indices above analog high have an analog low when shifted by -1 (rising edge)
ind_rising_edge = ind_above_ah(ind_shift_below_ah); % Store the original indices for the rising edges
% hold on
% plot(ind_rising_edge, ones(size(ind_rising_edge)) .* 5, 'o')
% hold off

stim_starts_gap = (P.Mcr_fcp.apis.seq_length_s - P.Mcr_fcp.apis.stim_length_s) * P.daqrate; % How long we expect the stim gap to be between the end of one stim to the start of the next
stim_prestart_baseline = (P.Mcr_fcp.apis.delay_time_ms / 1e3) * P.daqrate; % The duration between the baseline period and the corresponding stim start
stim_starts = ind_rising_edge([true; diff(ind_rising_edge) > stim_starts_gap]); % Add a 1/true at the beginning index for the first stim

% Plot the air puff signal and the calculated start points of each stim period
figure; plot(TD.airPuffOutput)
hold on
plot(stim_starts, ones(size(stim_starts)) .* 5, 'o')
hold off

clearvars ind_above_ah ind_shift_below_ah ind_rising_edge
% figure; plot(TD.sfTimeTagsDAQStart_adj) % plot the time tags for each superframe, adjusted to match the DAQ sampling rate

trial_windows = cell(size(stim_starts)); % Cell array of size (# trials, 1). Each cell contains the time points (according to the DAQ rate) that correspond to that trial.
trial_sf = cell(size(trial_windows));    % Cell array of size (# trials, 1). Each cell contains the superframe indices that started within that trial.

sfStarts = (TD.sfTimeTagsDAQStart_adj - TD.sfWidth_adj); % Adjust the superframe time tags so each index is at the start of the superframe acquisition

% Go through each trial within the run and assign the trial timepoints and the corresponding superframe indices
for trial = 1:length(trial_windows)
    trial_windows{trial} = stim_starts(trial) - stim_prestart_baseline : stim_starts(trial) + stim_starts_gap;

    trial_sf{trial} = find(sfStarts >= trial_windows{trial}(1) & sfStarts <= trial_windows{trial}(end));
end
clearvars trial

%% Remove outliers
% Use the "median" method of the filloutliers function
ro_fillmethod = "linear"; %
ro_findmethod = "percentiles";
quartile_threshold = [0, 99];
ro_dim = 4;

PDIallSF_ro = filloutliers(PDIallSF, ro_fillmethod, ro_findmethod, quartile_threshold, ro_dim);
CBViallSF_ro = filloutliers(CBViallSF, ro_fillmethod, ro_findmethod, quartile_threshold, ro_dim);
CBFsiallSF_ro = filloutliers(CBFsiallSF, ro_fillmethod, ro_findmethod, quartile_threshold, ro_dim);

%% Resample the trials for the hemodynamic parameters
interp_factor = 100;
% interp_factor = 1000;
[trial_CBVi_usi] = resampleTrials(CBViallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);
[trial_CBFsi_usi] = resampleTrials(CBFsiallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);
[trial_PDI_usi] = resampleTrials(PDIallSF, trial_sf, trial_windows, sfStarts, P, interp_factor);

% [trial_CBVi_usi] = resampleTrials(CBViallSF_ro, trial_sf, trial_windows, sfStarts, P, interp_factor);
% [trial_CBFsi_usi] = resampleTrials(CBFsiallSF_ro, trial_sf, trial_windows, sfStarts, P, interp_factor);
% [trial_PDI_usi] = resampleTrials(PDIallSF_ro, trial_sf, trial_windows, sfStarts, P, interp_factor);


%% Select usable trials
trials_to_remove_dlg = inputdlg('Enter space-separated trial numbers to remove:',...
             'Sample', [1 50]);
trials_to_remove = str2num(trials_to_remove_dlg{1});

% trial_CBVi_usi(trials_to_remove) = [];
% trial_CBFsi_usi(trials_to_remove) = [];
% trial_PDI_usi(trials_to_remove) = [];

trials_to_keep = setdiff(1:P.numTrials, trials_to_remove);

%% Calculate the relative hemodynamic changes for each trial

[trial_CBVi_usi_baseline, trial_rCBV_usi] = fUS_calc_rHP(trial_CBVi_usi(trials_to_keep), P, interp_factor);
[trial_CBFsi_usi_baseline, trial_rCBFspeed_usi] = fUS_calc_rHP(trial_CBFsi_usi(trials_to_keep), P, interp_factor);
[trial_PDI_usi_baseline, trial_rPDI_usi] = fUS_calc_rHP(trial_PDI_usi(trials_to_keep), P, interp_factor);

% [trial_CBVi_usi_baseline_alltrials, trial_rCBV_usi_alltrials] = fUS_calc_rHP(trial_CBVi_usi, P, interp_factor);
% [trial_CBFsi_usi_baseline_alltrials, trial_rCBFspeed_usi_alltrials] = fUS_calc_rHP(trial_CBFsi_usi, P, interp_factor);
% [trial_PDI_usi_baseline_alltrials, trial_rPDI_usi_alltrials] = fUS_calc_rHP(trial_PDI_usi, P, interp_factor);

%% Inspect the trials
% fUS_plotTrials(trial_rPDI_usi, [48, 68, 12])
% fUS_plotTrials(trial_rCBV_usi, [48, 68, 12])
% fUS_plotTrials(trial_rCBFspeed_usi, [48, 68, 12])

%% Trial average the relative hemodynamic changes

rCBV_TA = fUS_trialAverage(trial_rCBV_usi);
rCBFspeed_TA = fUS_trialAverage(trial_rCBFspeed_usi);
rPDI_TA = fUS_trialAverage(trial_rPDI_usi);

%% Correlation on the trial averaged rCBV

% Resample the stim pattern/predicted HRF
trial_stim_pattern = zeros(P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor, 1);
trial_stim_pattern(P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor : ...
    P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor + ...
    P.Mcr_fcp.apis.stim_length_s * P.daqrate / interp_factor) = 1;
figure; plot((1:length(trial_stim_pattern)) .* interp_factor ./ P.daqrate, trial_stim_pattern); title('Trial stim pattern'); xlabel('Time [s]')

% zt = 2;
zt = 2.58;
[r_rCBV, z_rCBV, am_rCBV] = activationMap3D(rCBV_TA, trial_stim_pattern, zt);

% volumeViewer(r_rCBV)
% volumeViewer(z_rCBV)
% volumeViewer(am_rCBV)
figure; imagesc(squeeze(max(r_rCBV(:, :, :), [], 1))'); colorbar; colormap jet; title('Correlation map coronal MIP'); clim([0, 1]) %clim([-1, 1])]
figure; imagesc(squeeze(max(z_rCBV(:, :, :), [], 1))'); colorbar; colormap jet; title('z-score map coronal MIP');
% figure; imagesc(squeeze(mean(z_rCBV(:, :, :), 1))'); colormap jet; clim([0, 1]) % clim([-1, 1])
% figure; imagesc(am_rCBV); colormap jet; title("Activation Map (rCBV) with z threshold = " + num2str(zt))
figure; imagesc(squeeze(max(am_rCBV(:, :, :), [], 1))'); colorbar; colormap jet; title("Activation Map (rCBV) coronal MIP with z threshold = " + num2str(zt))
figure; imagesc(squeeze(max(am_rCBV(:, :, :), [], 3))'); colorbar; colormap jet; title("Activation Map (rCBV) axial MIP with z threshold = " + num2str(zt))

% generateTiffStack_multi({r_rCBV}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({z_rCBV}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({am_rCBV}, [8.8, 8.8, 8], 'jet', 5)

%% Correlation on the trial averaged rCBFspeed

% Resample the stim pattern/predicted HRF
trial_stim_pattern = zeros(P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor, 1);
trial_stim_pattern(P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor : ...
    P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor + ...
    P.Mcr_fcp.apis.stim_length_s * P.daqrate / interp_factor) = 1;
figure; plot((1:length(trial_stim_pattern)) .* interp_factor ./ P.daqrate, trial_stim_pattern); title('Trial stim pattern'); xlabel('Time [s]')

% zt = 2;
zt = 2.58;

% [r_rCBFspeed, z_rCBFspeed, am_rCBFspeed] = activationMap3D_boxfilt(rCBFspeed_TA, trial_stim_pattern, zt);
[r_rCBFspeed, z_rCBFspeed, am_rCBFspeed] = activationMap3D(rCBFspeed_TA, trial_stim_pattern, zt);

% volumeViewer(r_rCBFspeed)
% volumeViewer(z_rCBFspeed)
% volumeViewer(am_rCBFspeed)
figure; imagesc(squeeze(max(r_rCBFspeed(:, :, :), [], 1))'); colorbar; colormap jet; title('Correlation map coronal MIP'); clim([0, 1]) %clim([-1, 1])]
figure; imagesc(squeeze(max(z_rCBFspeed(:, :, :), [], 1))'); colorbar; colormap jet; title('z-score map coronal MIP');
% figure; imagesc(squeeze(mean(z_rCBFspeed(:, :, :), 1))'); colormap jet; clim([0, 1]) % clim([-1, 1])
% figure; imagesc(am_rCBFspeed); colormap jet; title("Activation Map (rCBV) with z threshold = " + num2str(zt))
figure; imagesc(squeeze(max(am_rCBFspeed(:, :, :), [], 1))'); colorbar; colormap jet; title("Activation Map (rCBFspeed) coronal MIP with z threshold = " + num2str(zt))
figure; imagesc(squeeze(max(am_rCBFspeed(:, :, :), [], 3))'); colorbar; colormap jet; title("Activation Map (rCBFspeed) axial MIP with z threshold = " + num2str(zt))

% generateTiffStack_multi({r_rCBFspeed}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({z_rCBFspeed}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({am_rCBFspeed}, [8.8, 8.8, 8], 'jet', 5)

%% Correlation on the trial averaged rPDI

% Resample the stim pattern/predicted HRF
trial_stim_pattern = zeros(P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor, 1);
trial_stim_pattern(P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor : ...
    P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor + ...
    P.Mcr_fcp.apis.stim_length_s * P.daqrate / interp_factor) = 1;
figure; plot((1:length(trial_stim_pattern)) .* interp_factor ./ P.daqrate, trial_stim_pattern); title('Trial stim pattern'); xlabel('Time [s]')

% zt = 2;
zt = 2.58;

[r_rPDI, z_rPDI, am_rPDI] = activationMap3D(rPDI_TA, trial_stim_pattern, zt);

% volumeViewer(r_rPDI)
% volumeViewer(z_rPDI)
% volumeViewer(am_rPDI)
% figure; imagesc(squeeze(mean(r_rPDI(:, :, :), 1))'); colorbar; colormap jet; title('Correlation map coronal Mean IP'); clim([0, 1]) %clim([-1, 1])]

figure; imagesc(squeeze(max(r_rPDI(:, :, :), [], 1))'); colorbar; colormap jet; title('Correlation map coronal MIP'); clim([0, 1]) %clim([-1, 1])]
figure; imagesc(squeeze(max(z_rPDI(:, :, :), [], 1))'); colorbar; colormap jet; title('z-score map coronal MIP');
% figure; imagesc(squeeze(mean(z_rPDI(:, :, :), 1))'); colormap jet; clim([0, 1]) % clim([-1, 1])
% figure; imagesc(am_rPDI); colormap jet; title("Activation Map (rCBV) with z threshold = " + num2str(zt))
figure; imagesc(squeeze(max(am_rPDI(:, :, :), [], 1))'); colorbar; colormap jet; title("Activation Map (rPDI) coronal MIP with z threshold = " + num2str(zt))
figure; imagesc(squeeze(max(am_rPDI(:, :, :), [], 3) .^ 1)'); colorbar; colormap jet; title("Activation Map (rPDI) axial MIP with z threshold = " + num2str(zt))

% generateTiffStack_multi({r_rPDI}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({z_rPDI}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({am_rPDI}, [8.8, 8.8, 8], 'jet', 5)
% generateTiffStack_multi({squeeze(mean(PDIallSF, 4)) .^ 0.5, am_rPDI}, [8.8, 8.8, 8], 'jet', 5)

%% Plot ROIs defined by the half-max of the activation maps
numPtsUSI = P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor; % # of time points per trial for the upsampling

fraction = 0.75;
roi_indices_rPDI = roi_prop_max(am_rPDI, fraction);
roi_rPDI_TA = calc_ROI_avg(rPDI_TA, roi_indices_rPDI);
roi_indices_rCBV = roi_prop_max(am_rCBV, fraction);
roi_rCBV_TA = calc_ROI_avg(rCBV_TA, roi_indices_rCBV);
roi_indices_rCBFspeed = roi_prop_max(am_rCBFspeed, fraction);
roi_rCBFspeed_TA = calc_ROI_avg(rCBFspeed_TA, roi_indices_rCBFspeed);

figure; plot((1:length(roi_rPDI_TA)) .* interp_factor ./ P.daqrate, roi_rPDI_TA); xlabel('Time [s]'); ylabel('rPDI'); title("rPDI ROI timecourse")
figure; plot((1:length(roi_rCBV_TA)) .* interp_factor ./ P.daqrate, roi_rCBV_TA); xlabel('Time [s]'); ylabel('rCBV'); title("rCBV ROI timecourse")
figure; plot((1:length(roi_rCBFspeed_TA)) .* interp_factor ./ P.daqrate, roi_rCBFspeed_TA); xlabel('Time [s]'); ylabel('rCBFspeed'); title("rCBFspeed ROI timecourse")

% % Look at the max point
% [m, ind] = max(am_rPDI, [], 'all')
% [i, j, k] = ind2sub(size(am_rPDI), ind)
% 
% ts = 2; % test size
% testsection = rPDI_TA(i - ts : i + ts, j - ts : j + ts, k - ts : k + ts, :);
% % figure; plot(squeeze(rPDI_TA(i, j, k, :)))
% figure; plot(squeeze(mean(mean(mean(testsection, 1), 2), 3)))
% 
% [m, ind] = max(r_rPDI, [], 'all')
% [i, j, k] = ind2sub(size(am_rPDI), ind)
% figure; plot(squeeze(rPDI_TA(i, j, k, :)))
% 
% am_rPDI_t = am_rPDI; % thresholded
% am_rPDI_t(am_rPDI_t < 1.3) = 0;
% am_rPDI_t_roi_mask = am_rPDI_t > 1.3;

%% Plot median-averaged ROIs defined by the half-max of the activation maps
numPtsUSI = P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor; % # of time points per trial for the upsampling

fraction = 0.75;
median_filter_windowsize = 51;

roi_indices_rPDI = roi_prop_max(am_rPDI, fraction);
roi_rPDI_TA = calc_ROI_avg(movmedian(rPDI_TA, median_filter_windowsize, 4), roi_indices_rPDI);
roi_indices_rCBV = roi_prop_max(am_rCBV, fraction);
roi_rCBV_TA = calc_ROI_avg(movmedian(rCBV_TA, median_filter_windowsize, 4), roi_indices_rCBV);
roi_indices_rCBFspeed = roi_prop_max(am_rCBFspeed, fraction);
roi_rCBFspeed_TA = calc_ROI_avg(movmedian(rCBFspeed_TA, median_filter_windowsize, 4), roi_indices_rCBFspeed);

figure; plot((1:length(roi_rPDI_TA)) .* interp_factor ./ P.daqrate, roi_rPDI_TA); xlabel('Time [s]'); ylabel('rPDI'); title("rPDI ROI timecourse with median filter; window size = " + num2str(median_filter_windowsize))
figure; plot((1:length(roi_rCBV_TA)) .* interp_factor ./ P.daqrate, roi_rCBV_TA); xlabel('Time [s]'); ylabel('rCBV'); title("rCBV ROI timecourse with median filter; window size = " + num2str(median_filter_windowsize))
figure; plot((1:length(roi_rCBFspeed_TA)) .* interp_factor ./ P.daqrate, roi_rCBFspeed_TA); xlabel('Time [s]'); ylabel('rCBFspeed'); title("rCBFspeed ROI timecourse with median filter; window size = " + num2str(median_filter_windowsize))

