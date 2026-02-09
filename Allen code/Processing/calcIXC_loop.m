%% Description:
%       IXC processing for multiple superframes (volumetric data)

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

%% Define some parameters

parameterPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound'};
parameterDefaults = {'1', '', '20', ''};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});
numFiles = endFile - startFile + 1;
sv_threshold_lower = str2double(parameterUserInput{3});
sv_threshold_upper = str2double(parameterUserInput{4});

clearvars parameterPrompt parameterDefaults parameterUserInput


%% Main loop
for filenum = startFile:endFile
% for filenum = [2:endFile]
% for filenum = 18:endFile
% for filenum = 1

    % Load the IQ data
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    
    IQ = single(squeeze(IData + 1i .* QData));
%     clearvars IData QData

    % figure; imagesc(squeeze(max(abs(IQ(:, :, :, 2)), [], 1))')

    % Calculate the cross correlation of raw IQ (masked) to look at motion
    ixc = calcIXC_simple(IQ);
    ut_ms = (1:size(IQ, 4)) ./ P.frameRate .* 1e3; % micro time [ms]
%     % figure; plot(ut_ms, abs(ixc)); xlabel('Micro time [ms]'); ylabel('|Cross correlation of images|')
%     figure; plot(abs(ixc)); xlabel('Frame'); ylabel('|Cross correlation of images|')
%     title("Superframe " + num2str(filenum))


%     % Choose 2 frames to evaluate
%     ref_fn = 1;     % Reference frame #
%     moving_fn = 21; % Moving frame #
%     ref_vol = squeeze(IQ(:, :, :, ref_fn));
%     moving_vol = squeeze(IQ(:, :, :, moving_fn));
% 
%     % Upsample the images
%     us_factor = 4; % Upsampling factor
%     ref_vol_us = imresize3(ref_vol, us_factor, 'Method', 'cubic');
%     moving_vol_us = imresize3(moving_vol, us_factor, 'Method', 'cubic');
%     
%     % Test: look at MIPs of the upsampled IQ volumes
%     figure; imagesc(squeeze(max(abs(ref_vol_us), [], 1))'); colormap gray
%     figure; imagesc(squeeze(max(abs(moving_vol_us), [], 1))'); colormap gray
% 
%     % Try a few shifts (THESE MUST BE INTEGERS)
%     shift.inc = [1, 1, 1]; % y, x, z shift increments in units of [upsampled voxels]
%     % shift.max = [0, 0, 5]; % Maximum y, x, z |shift| in units of [upsampled voxels]
%     % shift.max = [5, 0, 2]; % Maximum y, x, z |shift| in units of [upsampled voxels]
%     shift.max = [2, 1, 1]; % Maximum y, x, z |shift| in units of [upsampled voxels]
%     shift.yspan = -shift.max(1):shift.inc(1):shift.max(1);
%     shift.xspan = -shift.max(2):shift.inc(2):shift.max(2);
%     shift.zspan = -shift.max(3):shift.inc(3):shift.max(3);
% 
%     [shift.ygrid, shift.xgrid, shift.zgrid] = meshgrid(shift.yspan,  shift.xspan,  shift.zspan);
%     % Squeeze in case the shift in any dimension is disabled, and vectorize
%     shift.ygrid = squeeze(shift.ygrid); shift.ygrid = shift.ygrid(:);
%     shift.xgrid = squeeze(shift.xgrid); shift.xgrid = shift.xgrid(:);
%     shift.zgrid = squeeze(shift.zgrid); shift.zgrid = shift.zgrid(:);
%     shift.shifts = [shift.ygrid, shift.xgrid, shift.zgrid];
% 
%     shift.numShifts = length(shift.ygrid); % Total # of shifts to try
% 
%     % shift.ixc = zeros(P.numFramesPerBuffer, shift.numShifts); % Initialize the post-shift ixc matrix. Each column is the ixc timecourse for that shift.
%     shift.ixc = zeros(1, shift.numShifts); % Initialize the post-shift ixc matrix. Each column is the ixc timecourse for that shift.
%     vs_us = size(ref_vol_us); % Upsampled volume's size
%     for sn = 1:shift.numShifts % shift number
%     % for sn = 1
%         disp(sn)
%         shift_sn = [shift.ygrid(sn), shift.xgrid(sn), shift.zgrid(sn)];
%         moving_vol_us_sn = imtranslate(moving_vol_us, shift_sn, 'OutputView','same'); % Shifted (upsampled) moving volume at shift number #sn
%         % figure; imagesc(squeeze(max(abs(moving_vol_us_sn), [], 1))'); colormap gray
%         
%         shift.ixc(:, sn) = calcIXC_shift(ref_vol_us, moving_vol_us_sn, true);
%     end   
%     shift.abs_ixc = abs(shift.ixc);








    % SVD decluttering
    [xp, yp, zp, nf] = size(IQ);
    PP = reshape(IQ, [xp*yp*zp, nf]);
    tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
    [U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
    SVs = diag(S);
%     disp('Full SVD done')

    % -- Some adaptive thresholding stuff -- %
    % Plot one SVD subspace as an image
%     subspace = 20;
%     subspace_img = reshape(U(:, subspace) * SVs(subspace) * V(:, subspace)', [xp, yp, zp, nf]);
%     figure; imagesc(squeeze(max(abs(subspace_img(:, :, :, 2)), [], 1))')
% %     volumeViewer(abs(subspace_img(:, :, :, 2)))
% 
    SSM = plotSSM(U, false);
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

%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'noise', '-v7.3', '-nocompression');
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'noise', 'SVs', 'SSM', '-v7.3')
    save([savepath, 'metrics-', num2str(filenum), '.mat'], 'ixc', 'SVs', 'SSM', '-v7.3')

%     disp("fUS result for file " + num2str(filenum) + " saved" )
    disp("info for file " + num2str(filenum) + " saved" )
%     disp("g1 result for file " + num2str(filenum) + " saved" )

    toc
    
end

%% Make a plot of all the ixcs
ixc_allfiles = [];
SVs_allfiles = [];
for filenum = startFile:endFile
    load([savepath, 'metrics-', num2str(filenum), '.mat'])
    ixc_allfiles = cat(2, ixc_allfiles, ixc);
    SVs_allfiles = cat(2, SVs_allfiles, SVs);
end

% figure; plot(abs(ixc_allfiles)); xlabel("Frame"); ylabel('|Cross correlation of images|')
figure; plot(ut_ms, abs(ixc_allfiles)); xlabel("Micro time [ms]"); ylabel('|Cross correlation of images|')
figure; semilogy(SVs_allfiles); xlabel("Singular value number"); ylabel("Singular value magnitude")

%% Calculate metrics for how often or largely the XC drops for each "superframe"
ixc_threshold = 0.9; % Threshold for XC dropping to be "significant"
ixc_under_threshold_mat = ixc_allfiles < ixc_threshold;
ixc_under_threshold = sum(ixc_under_threshold_mat, 1); % How many time points for each superframe is the XC under the threshold
ixc_min = min(ixc_allfiles, [], 1); % The minimum ixc value for each superframe

figure; plot(ixc_under_threshold); xlabel("Superframe #"); ylabel("Number of frames under threshold = " + num2str(ixc_threshold))

%% Save the ixc and singular value data
save([savepath, 'ixc_SV_data.mat'], "ixc_allfiles", "SVs_allfiles", "ixc_threshold", "ixc_under_threshold", "ixc_min")