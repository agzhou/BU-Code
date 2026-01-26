%% Description:
%   - 3D fUS and functional connectivity (FC) processing
%   - For use with continuously acquired data
%   - Timing data should be processed and saved with plotfUStiming_FC.m first

clearvars
%% load params and stuff
IQpath = uigetdir('D:\Allen\Data\', 'Select the IQ data path');
IQpath = [IQpath, '\'];

% Load parameters
% if ~exist('P', 'var')
%     load([IQpath, '..\params.mat'])
% end
% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [IQpath, '..\params.mat']);
    load([params_pathname, params_filename])
end

% Load Verasonics reconstruction parameters: datapath\PData.mat
if ~exist('PData', 'var')
    load([IQpath, 'PData.mat'])
end

IQfilenameStructure = ['IQ-', num2str(round(P.maxAngle)), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

savepath = uigetdir(IQpath + "..\", 'Select the save path');
savepath = [savepath, '\'];

addpath([cd, '\..\']) % Add the main "Processing" path

% Load the timing data
[timingFilePathFN, timingFilePath] = uigetfile([IQpath, '..\Timing data\TD.mat'], 'Select the timing data');
timingFilePath = [timingFilePath, timingFilePathFN];
load(timingFilePath)
% load(timingFilePath, 'acqStart', 'airPuffOutput', 'daqStartTimetag', 'sfTimeTags', 'sfTimeTagsDAQStart', 'sfTimeTagsDAQStart_adj', 'sfWidth', 'sfWidth_adj', 'timeStamp')

%% Define some parameters
max_freq_expected = 0.4; % Expected maximum frequency of the relevant signals[Hz]

parameterPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound', 'Block size [frames]', 'Block overlap %'};
parameterDefaults = {'1', '', '20', '', num2str(TD.sfWidth / max_freq_expected * P.numFramesPerBuffer), '50'};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});
numFiles = endFile - startFile + 1;
sv_threshold_lower = str2double(parameterUserInput{3});
sv_threshold_upper = str2double(parameterUserInput{4});
bs = str2double(parameterUserInput{5});        % Block size [frames]; the default is calculated to equal the length of 1 period of the highest frequency component expected (max_freq_expected)
if mod(bs, 2) ~= 0 % Check if the block size is even
    error('Block size must be even')
end
bo = str2double(parameterUserInput{6}) ./ 100; % Block overlap in decimal format

clearvars parameterPrompt parameterDefaults parameterUserInput

%% Look at the IQ MIP to get an idea for where to crop

% Load the first IQ data file
tic
load([IQpath, IQfilenameStructure, num2str(1)])

IQ = single(squeeze(IData + 1i .* QData));
clearvars IData QData

figure; imagesc(squeeze(max(abs(IQ(:, :, :, 2)), [], 1))')
% volumeViewer(squeeze(abs(IQ(:, :, :, 2))))

%% Prompt for the z crop
zCropPrompt = {'z start (voxel number)', 'z end (voxel number)'};
zCropDefaults = {'1', num2str(size(IQ, 3))};
zCropUserInput = inputdlg(zCropPrompt, 'Input Parameters', 1, zCropDefaults);

%     zstart = 45; zend = 135;
zstart = str2double(zCropUserInput{1});
zend = str2double(zCropUserInput{2});

clearvars zCropPrompt zCropDefaults zCropUserInput

%% Define how to slice each block
% nsfpb = bs / P.numFramesPerBuffer; % # of superframes needed per block
% nfpsfib = P.numFramesPerBuffer .* ones(ceil(nsfpb), 1); nfpsfib(end) = P.numFramesPerBuffer .* (nsfpb - floor(nsfpb)); % # of frames per superframe in each block
nfpbo = bs * bo; % # of frames per block overlap
if nfpbo ~= floor(nfpbo) % If nfpbo is not a natural number
    error("Block overlap must a natural number of frames")
end

numBlocks = floor( (numFiles * P.numFramesPerBuffer / bs) / (1-bo) ) - 1; % **** CHECK THIS ****

%% Set up the High Pass Filter
% fc = 50; % Cutoff frequency [Hz]
% fs = P.frameRate; % Sampling frequency [Hz]
% HPF_order = 3; % Butterworth filter order
% 
% [HPF_b, HPF_a] = butter(HPF_order, fc/(fs/2), 'high');

%% Save proc params
% numg1pts = 20; % Only calculate the first N points
% save([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'numg1pts', 'zstart', 'zend');
save([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'zstart', 'zend', 'bs', 'bo', 'max_freq_expected');

% Add band pass filter params later............

%% Main loop: go through each block
for bn = 1:numBlocks
% for bn = 237:numBlocks
% for bn = [1:4, 6:numBlocks]
    tic
    % Define which frame numbers (relative to the experiment start) should be used
    if bn == 1
        frames_bn = 1:bs;
    else
        frames_bn = ( (bn - 1)*bs + 1 : (bn)*bs ) - nfpbo*(bn - 1);
    end

    % Define which superframes (and which portions of each) to load and use
    % for each block
    sf_start_bn = ceil(frames_bn(1) ./ P.numFramesPerBuffer); % The superframe to start on for block bn (out of the whole experiment)
    if bn == 1 & bs > P.numFramesPerBuffer
        fracOfStartSFToUse = 1; % Special case for the first block
    else
        fracOfStartSFToUse = ceil( (frames_bn(1) - 1)./ P.numFramesPerBuffer ) - (frames_bn(1) - 1)./ P.numFramesPerBuffer; % Fraction of the first superframe to use in block bn (starting from the end of the superframe)
    end
    numFramesOfStartSFToUse = fracOfStartSFToUse * P.numFramesPerBuffer; % # of frames in the first superframe to use in block bn (starting from the end of the superframe)
    numFullSFToUseAfterStartSF = floor( (bs - numFramesOfStartSFToUse)/P.numFramesPerBuffer ); % # of full superframes to use after the first superframe
    if floor(numFullSFToUseAfterStartSF) == (bs - numFramesOfStartSFToUse)/P.numFramesPerBuffer % If there is no need for a partial end superframe
        numFramesPerSFToUse = [numFramesOfStartSFToUse, P.numFramesPerBuffer .* ones(1, numFullSFToUseAfterStartSF)]';
    else
        numFramesOfEndSFToUse = bs - numFramesOfStartSFToUse - numFullSFToUseAfterStartSF * P.numFramesPerBuffer;
        numFramesPerSFToUse = [numFramesOfStartSFToUse, P.numFramesPerBuffer .* ones(1, numFullSFToUseAfterStartSF), numFramesOfEndSFToUse]';
    end

    % Sometimes there might be a zero at the beginning if the multiples are clean
    if numFramesPerSFToUse(1) == 0
        numFramesPerSFToUse = numFramesPerSFToUse(2:end);
    end

    IQ = [];
    for sfi = sf_start_bn:sf_start_bn + length(numFramesPerSFToUse) - 1 % Go through and load each superframe, with slicing
        % Load the IQ data
        load([IQpath, IQfilenameStructure, num2str(sfi)])
        
        IQ_sfi = single(squeeze(IData + 1i .* QData));
        clearvars IData QData
        
        if sfi == sf_start_bn % Special case if it's the starting superframe, where the ending chunk needs to be added
            IQ = cat(4, IQ, IQ_sfi(:, :, :, P.numFramesPerBuffer - numFramesPerSFToUse(sfi - sf_start_bn + 1) + 1 : end));
        else
            IQ = cat(4, IQ, IQ_sfi(:, :, :, 1:numFramesPerSFToUse(sfi - sf_start_bn + 1)));
        end

    end

    % figure; imagesc(squeeze(max(abs(IQ(:, :, :, 2)), [], 1))')
    
    % Crop the IQ first 
    IQm = IQ(:, :, zstart:zend, :);

%     figure; imagesc(squeeze(max(abs(IQm(:, :, :, 2)), [], 1))')

    %%%%%%%%%%%%%% IF USING THE PREDEFINED MASK %%%%%%%%%%%%
%     IQm = IQ;
%     IQm(coronal_mask_rep) = 0; % Apply the brain mask to the IQ: set the non-brain voxels equal to 0

    % Apply the HPF
%     dim = length(size(IQm)); % Operate on the time dimension
%     IQm = filter(HPF_b, HPF_a, IQm, [], dim);

    % SVD decluttering
%     [PP, EVs, V_sort] = getSVs2D(IQm);
    [xp, yp, zp, nf] = size(IQm);
    PP = reshape(IQm, [xp*yp*zp, nf]);
    % tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
    [U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
    SVs = diag(S);
%     disp('Full SVD done')
    % toc
    % disp('SVs decomposed')

    % -- Some adaptive thresholding stuff -- %
    % Plot one SVD subspace as an image
%     subspace = 20;
%     subspace_img = reshape(U(:, subspace) * SVs(subspace) * V(:, subspace)', [xp, yp, zp, nf]);
%     figure; imagesc(squeeze(max(abs(subspace_img(:, :, :, 2)), [], 1))')
% %     volumeViewer(abs(subspace_img(:, :, :, 2)))
% 
%     SSM = plotSSM(U, false);
% %     SSM = plotSSM(U, true);
%     [~, a_opt, b_opt] = fitSSM(SSM, false); % Get the optimal singular value thresholds
% %     [~, a_opt, b_opt] = fitSSM(SSM, true); % Get the optimal singular value thresholds
%     

    [IQf, noise] = applySVs2D(IQm, PP, SVs, V, sv_threshold_lower, sv_threshold_upper);
%     [IQf, noise] = applySVs2D(IQm, PP, SVs, V, a_opt, b_opt);
    % disp('SVD filtered images put together')

%     volumeViewer(abs(IQf(:, :, :, 1)))
%     figure; imagesc(squeeze(abs(max(IQf(:, :, :, 1), [], 1)))'); colorbar
%     generateTiffStack_acrossframes(abs(IQf), [8.8, 8.8, 8], 'hot', 1:80)
    % clearvars IQ

    % Use the IQf with separated negative and positive frequency components
%     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
   
%     [PDI] = calcPowerDoppler(IQf_separated);
    PDI = sum(abs(IQf) .^ 2, 4) ./ size(IQf, 4);
%     [CDI] = calcColorDoppler(IQf_FT_separated, P);

%     figure; imagesc(squeeze(max(PDI, [], 1))' .^ 0.5); colormap hot; colorbar
%     figure; imagesc(squeeze(max(PDI ./ noise, [], 1))' .^ 0.5); colormap hot; colorbar
%     volumeViewer(PDI)
%     volumeViewer(PDI ./ noise)

%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'noise', '-v7.3', '-nocompression');
    save([savepath, 'fUSdata-', num2str(bn), '.mat'], 'PDI', 'noise', 'SVs', 'numFramesPerSFToUse', '-v7.3')

    disp("fUS result for block " + num2str(bn) + " saved" )
%     disp("g1 result for file " + num2str(filenum) + " saved" )

    toc
    
end

save([savepath, 'blocking_info.mat'], 'bn', 'bo', 'bs', 'startFile', 'endFile', 'nf', 'nfpbo', 'numBlocks', 'numFiles')


%% Store all the PDI across the experiment into one matrix
load([savepath, 'fUSdata-', num2str(1), '.mat'], 'PDI', 'noise')
PDIallBlocks = zeros([size(PDI), numBlocks]); % Matrix with the CBVi for every superframe
PDIallBlocks(:, :, :, 1) = PDI ./ noise;

for bn = 1:numBlocks
%     load([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI')
    load([savepath, 'fUSdata-', num2str(bn), '.mat'], 'PDI')

    PDIallBlocks(:, :, :, bn) = PDI ./ noise;
end

%% Save the PDIallBlocks
save([savepath, 'PDIallBlocks.mat'], "PDIallBlocks")

%% Visualize the PDI and CDI across the experiment
% generateTiffStack_acrossframes(PDIallSF{3} .^ 0.7, [8.8, 8.8, 8], 'hot', 1:80)
generateTiffStack_acrossframes(PDIallBlocks .^ 0.5, [8.8, 8.8, 8], 'hot', 1:80)

%% Check different MIPs across superframes
% yr = 20:40;
yr = 1:80;
generateTiffStack_acrossframes(PDIallBlocks .^ 0.5, [8.8, 8.8, 8], 'hot', yr)
% generateTiffStack_acrossframes(CBViallSF .^ 1, [8.8, 8.8, 8], 'hot', yr)

%% Prepare template(s) for atlas registration
% Create templates for each hemodynamic parameter, averaging across superframes
PDI_allSF_avg = mean(PDIallBlocks, 4);

voxel_size = PData.PDelta .* P.wl; % Voxel size (y, x, z) in meters
fUS_volume_dimensions_m = [P.Trans.numelements/2 * P.Trans.spacingMm / 1e3, P.Trans.numelements/2 * P.Trans.spacingMm / 1e3, (P.endDepthMM - P.startDepthMM)/1e3]; % Volume size in meters
fUS_volume_dimensions_voxels = PData.Size; % Volume size in voxels (from the recon PData)

% Adjust the sizes based on the pre-SVD/clutter filtering cropping
fUS_cropped_volume_dimensions_voxels = size(PDI_allSF_avg);
fUS_cropped_volume_dimensions_m = fUS_cropped_volume_dimensions_voxels ./ fUS_volume_dimensions_voxels .* fUS_volume_dimensions_m;

% User input for target voxel size (post-interpolation)
targetVoxelSizePrompt = {'y Target Voxel Size [um]', 'x Target Voxel Size [um]', 'z Target Voxel Size [um]'};
% targetVoxelSizeDefaults = {'10', '10', '10'};
targetVoxelSizeDefaults = {'50', '50', '50'};
targetVoxelSizeUserInput = inputdlg(targetVoxelSizePrompt, 'Input Target Voxel Size', 1, targetVoxelSizeDefaults);

% Store target voxel size inputs and convert to meters
target_voxel_size(1) = str2double(targetVoxelSizeUserInput{1}) ./ 1e6;
target_voxel_size(2) = str2double(targetVoxelSizeUserInput{2}) ./ 1e6;
target_voxel_size(3) = str2double(targetVoxelSizeUserInput{3}) ./ 1e6;

prereg_interp_factor = voxel_size ./ target_voxel_size;

% Resample hemodynamic parameter template maps to the desired voxel size
PDI_allSF_avg_rs = imresize3(PDI_allSF_avg, 'Scale', prereg_interp_factor, 'Method', 'cubic');

% Store pre-registration parameters
prereg_params.orig_voxel_size = voxel_size;
prereg_params.fUS_volume_dimensions_m = fUS_volume_dimensions_m;
prereg_params.fUS_volume_dimensions_voxels = fUS_volume_dimensions_voxels;
prereg_params.fUS_cropped_volume_dimensions_voxels = fUS_cropped_volume_dimensions_voxels;
prereg_params.fUS_cropped_volume_dimensions_m = fUS_cropped_volume_dimensions_m;
prereg_params.target_voxel_size = target_voxel_size;
prereg_params.prereg_interp_factor = prereg_interp_factor;
% prereg_params. = 

save([savepath, 'prereg_PDI_params_50um.mat'], "PDI_allSF_avg_rs", "prereg_params")

%% Prepare template(s) for atlas registration (using the preloaded prereg_params struct)
% Create templates for each hemodynamic parameter, averaging across superframes
PDI_allSF_avg = mean(PDIallBlocks, 4);

% voxel_size = PData.PDelta .* P.wl; % Voxel size (y, x, z) in meters
% fUS_volume_dimensions_m = [P.Trans.numelements/2 * P.Trans.spacingMm / 1e3, P.Trans.numelements/2 * P.Trans.spacingMm / 1e3, (P.endDepthMM - P.startDepthMM)/1e3]; % Volume size in meters
% fUS_volume_dimensions_voxels = PData.Size; % Volume size in voxels (from the recon PData)

% Adjust the sizes based on the pre-SVD/clutter filtering cropping
fUS_cropped_volume_dimensions_voxels = size(PDI_allSF_avg);
fUS_cropped_volume_dimensions_m = fUS_cropped_volume_dimensions_voxels ./ prereg_params.fUS_volume_dimensions_voxels .* prereg_params.fUS_volume_dimensions_m;

% User input for target voxel size (post-interpolation)
targetVoxelSizePrompt = {'y Target Voxel Size [um]', 'x Target Voxel Size [um]', 'z Target Voxel Size [um]'};
% targetVoxelSizeDefaults = {'10', '10', '10'};
targetVoxelSizeDefaults = {'50', '50', '50'};
targetVoxelSizeUserInput = inputdlg(targetVoxelSizePrompt, 'Input Target Voxel Size', 1, targetVoxelSizeDefaults);

% Store target voxel size inputs and convert to meters
target_voxel_size(1) = str2double(targetVoxelSizeUserInput{1}) ./ 1e6;
target_voxel_size(2) = str2double(targetVoxelSizeUserInput{2}) ./ 1e6;
target_voxel_size(3) = str2double(targetVoxelSizeUserInput{3}) ./ 1e6;

prereg_interp_factor = prereg_params.orig_voxel_size ./ target_voxel_size;

% Resample hemodynamic parameter template maps to the desired voxel size
PDI_allSF_avg_rs = imresize3(PDI_allSF_avg, 'Scale', prereg_interp_factor, 'Method', 'cubic');

% Store pre-registration parameters
% prereg_params.orig_voxel_size = voxel_size;
% prereg_params.fUS_volume_dimensions_m = fUS_volume_dimensions_m;
% prereg_params.fUS_volume_dimensions_voxels = fUS_volume_dimensions_voxels;
% prereg_params.fUS_cropped_volume_dimensions_voxels = fUS_cropped_volume_dimensions_voxels;
% prereg_params.fUS_cropped_volume_dimensions_m = fUS_cropped_volume_dimensions_m;
prereg_params.target_voxel_size = target_voxel_size;
prereg_params.prereg_interp_factor = prereg_interp_factor;
% prereg_params. = 























%% Helper functions

function [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF)

    g1A_T = {};
    
    g1A_T{1} = abs(g1(:, :, :, 2)) > g1_tau1_cutoff; % First treatment: tau1 is above some cutoff (make sure there is some actual blood signal there)
%     g1A_T{2} = abs(g1(:, :, :, tau1_index_CBF)) > abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > |g1(tau2)| (noise might have the g1 randomly increase with tau, but it should not happen with a voxel where there is a real blood signal)
%     g1A_T{3} = abs(g1(:, :, :, tau1_index_CBF)) > 2 .* abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
%     % g1A_T{4} = abs(g1(:, :, :, tau1_index_CBF)) - 1 .* abs(g1(:, :, :, tau2_index_CBF)) > tau_difference_cutoff; % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    
    g1A_mask = true(size(g1A_T{1})); % Mask of voxels to keep for the g1 treatments
    for i = 1:length(g1A_T)
        g1A_mask = and(g1A_mask, g1A_T{i});
    end

end

