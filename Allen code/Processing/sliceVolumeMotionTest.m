% Description: look at different slices of a volume over time and see if
% the motion is localized

% Motivation: the linear array showed fewer motion artifacts than the RCA.
% Hypothesis: there is localized motion (throat/tongue, blinking) that the
% linear array was positioned away from

clearvars
%% Define parameters
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

%% Define some parameters

parameterPrompt = {'Start file number', 'End file number', 'Number of files to stitch together', 'Width of slices [mm]'};
parameterDefaults = {'1', '', '', '1.1'};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});
numFiles = endFile - startFile + 1;
nfpb = str2double(parameterUserInput{3}); % # files per "block"
sliceWidthMM = str2double(parameterUserInput{4});
sliceWidth = sliceWidthMM / 1e3;

clearvars parameterPrompt parameterDefaults parameterUserInput

%% Calculate the number of total "blocks" to use
numBlocks = floor(numFiles / nfpb);
save([savepath, 'blocking_info.mat'], 'startFile', 'endFile', 'numFiles', 'nfpb', 'numBlocks')

%% Calculate the spatial slicing (coronal slices)
voxel_size = PData.PDelta .* P.wl; % Voxel size (y, x, z) in meters (from beamforming)
nyp = PData.Size(1); % # of y pixels
volSize = [P.numElements.*P.Trans.spacingMm, P.numElements.*P.Trans.spacingMm, P.endDepthMM - P.startDepthMM] ./ 1e3; % Volume extents [m]

numSlices = floor(volSize(1) / sliceWidth);
csm = cell(numSlices, 1); % Initialize coronal slice masks variable
numVoxelsPerSlice = floor(sliceWidth / voxel_size(1));

for si = 1:numSlices
    csm{si} = (si - 1) * numVoxelsPerSlice + 1 : (si) * numVoxelsPerSlice;
end

save([savepath, 'slicing_info.mat'], 'csm', 'numSlices', 'numVoxelsPerSlice', 'volSize', 'voxel_size', 'sliceWidth', 'sliceWidthMM')

%% Main loop
for bi = 1:numBlocks
% for bi = [81:numBlocks]
% for bi = 1
    IQ = [];
    filenumsToUse = (bi - 1) * nfpb + 1 : bi * nfpb;

    % Load the IQ data
    tic
    for fn = filenumsToUse
        IQtemp = load([IQpath, IQfilenameStructure, num2str(fn)], 'IQ').('IQ');
        IQ = cat(4, IQ, IQtemp);
    end
    clearvars IQtemp

    % figure; imagesc(squeeze(max(abs(IQ(:, :, :, 2)), [], 1))')

    ixc = zeros(size(IQ, 4), numSlices);
    SVs = zeros(size(IQ, 4), numSlices);
    % Go through each slice for that IQ ensemble
    for si = 1:numSlices
        IQs = IQ(csm{si}, :, :, :);

        % Calculate the cross correlation of raw IQ (masked) to look at motion
        ixc(:, si) = calcIXC_simple(IQs);

        % SVD decluttering
        [xp, yp, zp, nf] = size(IQs);
        PP = reshape(IQs, [xp*yp*zp, nf]);
    %     [U, S, V] = svd(PP); % Already sorted in decreasing order
        [U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
        SVs(:, si) = diag(S);

    end
    % save([savepath, 'ixc-sliced-block', num2str(bi)], 'ixc')
    save([savepath, 'metrics-sliced-block', num2str(bi)], 'ixc', 'SVs')
     % figure; plot(abs(ixc)); xlabel('Frame'); ylabel('|Cross correlation of images|')
     % figure; plot(abs(ixc) - abs(ixc(2, :))); xlabel('Frame'); ylabel('|Cross correlation of images|')



    % % Choose 2 frames to evaluate
    % ref_fn = 1;     % Reference frame #
    % moving_fn = 251; % Moving frame #
    % ref_vol = squeeze(IQ(:, :, :, ref_fn));
    % moving_vol = squeeze(IQ(:, :, :, moving_fn));
    % % % Try using the phase
    % % ref_vol = squeeze(angle(IQ(:, :, :, ref_fn)));
    % % moving_vol = squeeze(angle(IQ(:, :, :, moving_fn)));
    % 
    % % Upsample the images
    % us_factor = 1; % Upsampling factor
    % ref_vol_us = imresize3(ref_vol, us_factor, 'Method', 'cubic');
    % moving_vol_us = imresize3(moving_vol, us_factor, 'Method', 'cubic');
    % 
    % % Test: look at MIPs of the upsampled IQ volumes
    % figure; imagesc(squeeze(max(abs(ref_vol_us), [], 1))'); colormap gray
    % figure; imagesc(squeeze(max(abs(moving_vol_us), [], 1))'); colormap gray
    % 
    % % Try a few shifts (THESE MUST BE INTEGERS)
    % shift.inc = [1, 1, 1]; % y, x, z shift increments in units of [upsampled voxels]
    % shift.max = [5, 5, 5]; % Maximum y, x, z |shift| in units of [upsampled voxels]
    % % shift.max = [0, 0, 5]; % Maximum y, x, z |shift| in units of [upsampled voxels]
    % % shift.max = [2, 1, 1]; % Maximum y, x, z |shift| in units of [upsampled voxels]
    % shift.yspan = -shift.max(1):shift.inc(1):shift.max(1);
    % shift.xspan = -shift.max(2):shift.inc(2):shift.max(2);
    % shift.zspan = -shift.max(3):shift.inc(3):shift.max(3);
    % 
    % [shift.ygrid, shift.xgrid, shift.zgrid] = meshgrid(shift.yspan,  shift.xspan,  shift.zspan);
    % % Squeeze in case the shift in any dimension is disabled, and vectorize
    % shift.ygrid = squeeze(shift.ygrid); shift.ygrid = shift.ygrid(:);
    % shift.xgrid = squeeze(shift.xgrid); shift.xgrid = shift.xgrid(:);
    % shift.zgrid = squeeze(shift.zgrid); shift.zgrid = shift.zgrid(:);
    % shift.shifts = [shift.ygrid, shift.xgrid, shift.zgrid];
    % 
    % shift.numShifts = length(shift.ygrid); % Total # of shifts to try
    % 
    % % shift.ixc = zeros(P.numFramesPerBuffer, shift.numShifts); % Initialize the post-shift ixc matrix. Each column is the ixc timecourse for that shift.
    % shift.ixc = zeros(1, shift.numShifts); % Initialize the post-shift ixc matrix. Each column is the ixc timecourse for that shift.
    % vs_us = size(ref_vol_us); % Upsampled volume's size
    % for sn = 1:shift.numShifts % shift number
    % % for sn = 1
    %     disp(sn)
    %     shift_sn = [shift.ygrid(sn), shift.xgrid(sn), shift.zgrid(sn)];
    %     moving_vol_us_sn = imtranslate(moving_vol_us, shift_sn, 'OutputView','same'); % Shifted (upsampled) moving volume at shift number #sn
    %     % figure; imagesc(squeeze(max(abs(moving_vol_us_sn), [], 1))'); colormap gray
    % 
    %     shift.ixc(:, sn) = calcIXC_shift(ref_vol_us, moving_vol_us_sn, true);
    % end   
    % 
    % shift.abs_ixc = abs(shift.ixc);
    % [shift.opt_shift_ixc, shift.opt_shift_ind] = max(shift.abs_ixc);
    % shift.opt_shift = shift.shifts(shift.opt_shift_ind, :);
    % 
    % figure; plot(shift.abs_ixc)
    % figure; plot(shift.shifts(:, 3), shift.abs_ixc); xlabel('z shift [voxels]'); ylabel('Cross correlation') % z shift only
    % 
    % %% Compare phase of the moving and fixed frames
    % generateTiffStack_acrossframes(abs(IQ) .^ 0.3, [8.8, 8.8, 8], 'gray', 1:size(IQ, 1))
    % 
    % figure; imagesc(squeeze(angle(ref_vol(size(ref_vol, 1)/2, :, :)))'); colormap jet; xlabel('x'); ylabel('z'); colorbar; title('Phase [rad]')
    % figure; imagesc(squeeze(angle(moving_vol(size(ref_vol, 1)/2, :, :)))'); colormap jet; xlabel('x'); ylabel('z'); colorbar; title('Phase [rad]')
    % figure; imagesc(squeeze(angle(ref_vol_us(size(ref_vol, 1)/2, :, :)))'); colormap jet; xlabel('x'); ylabel('z'); colorbar; title('Phase [rad]')
    % figure; imagesc(squeeze(angle(moving_vol_us(size(ref_vol, 1)/2, :, :)))'); colormap jet; xlabel('x'); ylabel('z'); colorbar; title('Phase [rad]')
    % 
    % %% TEST: Translate the frame, downsample back, and then do SVD
    % 

    

    % -- Some adaptive thresholding stuff -- %
    % Plot one SVD subspace as an image
%     subspace = 20;
%     subspace_img = reshape(U(:, subspace) * SVs(subspace) * V(:, subspace)', [xp, yp, zp, nf]);
%     figure; imagesc(squeeze(max(abs(subspace_img(:, :, :, 2)), [], 1))')
% %     volumeViewer(abs(subspace_img(:, :, :, 2)))
% 
    % SSM = plotSSM(U, false);
% %     SSM = plotSSM(U, true);
%     [~, a_opt, b_opt] = fitSSM(SSM, false); % Get the optimal singular value thresholds
% %     [~, a_opt, b_opt] = fitSSM(SSM, true); % Get the optimal singular value thresholds
%     

%     [IQf, noise] = applySVs2D(IQm, PP, SVs, V, sv_threshold_lower, sv_threshold_upper);
% %     [IQf, noise] = applySVs2D(IQm, PP, SVs, V, a_opt, b_opt);
%     disp('SVD filtered images put together')
% 
% %     volumeViewer(abs(IQf(:, :, :, 1)))
% %     figure; imagesc(squeeze(abs(max(IQf(:, :, :, 1), [], 1)))'); colorbar
% %     generateTiffStack_acrossframes(abs(IQf), [8.8, 8.8, 8], 'hot', 1:80)
%     % clearvars IQ
% 
%     % Use the IQf with separated negative and positive frequency components
% %     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
% 
% %     [PDI] = calcPowerDoppler(IQf_separated);
%     PDI = sum(abs(IQf) .^ 2, 4) ./ size(IQf, 4);
% %     [CDI] = calcColorDoppler(IQf_FT_separated, P);
% 
% %     figure; imagesc(squeeze(max(PDI, [], 1))' .^ 0.5); colormap hot; colorbar
% %     figure; imagesc(squeeze(max(PDI ./ noise, [], 1))' .^ 0.5); colormap hot; colorbar
% %     volumeViewer(PDI)
% %     volumeViewer(PDI ./ noise)

    % save([savepath, 'metrics-', num2str(bi), '.mat'], 'ixc', 'SVs', 'SSM', '-v7.3')

%     disp("fUS result for file " + num2str(filenum) + " saved" )
    disp("info for file " + num2str(bi) + " saved" )
%     disp("g1 result for file " + num2str(filenum) + " saved" )

    toc
    
end

%% Store each block's metrics (THIS WORKS FOR NON-OVERLAPPING BLOCKS ONLY)
ixc_allblocks = zeros(P.numFramesPerBuffer * numBlocks, numSlices); % Make an ixc matrix. Each column is the ixc of one slice over the whole experiment's timecourse.
for bi = 1:numBlocks
% for bi = [81:numBlocks]
% for bi = 1


    metrics_bi = load([savepath, 'metrics-sliced-block', num2str(bi)]);
    ixc_allblocks((bi - 1)*P.numFramesPerBuffer + 1:bi*P.numFramesPerBuffer, :) = metrics_bi.ixc;
end

%% Plot the metrics
dt = 1/P.frameRate; % Time step calculated from the block size, overlap, and frame rate
t = 0:dt:(numBlocks*P.numFramesPerBuffer - 1)*dt; % Time stamps of each block

ixcmin = min(abs(ixc_allblocks), [], 'all');
ixcmax = max(abs(ixc_allblocks), [], 'all');
figure
sp = stackedplot(t, abs(ixc_allblocks));
for ai = 1:size(sp.AxesProperties, 1) % Go through the index for each axes object
    sp.AxesProperties(ai).YLimits = [ixcmin, ixcmax];
end
title("Volume cross correlation over different coronal slices (width = " + num2str(sliceWidthMM) + "mm)")
xlabel('Time [s]')

%% Make a video of the mouse in some slices over one superframe
testSliceInd = 4;
testSlice = squeeze(abs(IQ(csm{testSliceInd}, :, :, :)));
generateTiffStack_acrossframes(testSlice .^ 0.3, [sliceWidth, volSize(2), volSize(3)], 'gray', 1:size(testSlice, 1))