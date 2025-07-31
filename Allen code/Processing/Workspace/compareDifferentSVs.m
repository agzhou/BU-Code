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

IQfilenameStructure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

savepath = uigetdir('D:\Allen\Data\', 'Select the save path');
savepath = [savepath, '\'];

addpath([cd, '\Speckle tracking']) % add path for the g1 calculation functions

% Load the timing data
% [timingFilePathFN, timingFilePath] = uigetfile([IQpath, '..\Timing data\TD.mat'], 'Select the timing data');
% timingFilePath = [timingFilePath, timingFilePathFN];
% load(timingFilePath)

%% Define some parameters

parameterPrompt = {'Start file number', 'End file number', 'Tau 1 index for CBFspeed', 'Tau 2 index for CBFspeed', 'Tau 1 index for CBV'};
parameterDefaults = {'1', '', '2', '10', '2'};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});
numFiles = endFile - startFile + 1;
tau1_index_CBF = str2double(parameterUserInput{3});
tau2_index_CBF = str2double(parameterUserInput{4});
tau1_index_CBV = str2double(parameterUserInput{5});

clearvars parameterPrompt parameterDefaults parameterUserInput

taustep = 1/P.frameRate;
% tau = taustep:taustep:(P.numFramesPerBuffer * taustep);
tau = 0:taustep:((P.numFramesPerBuffer - 1) * taustep);
tau_ms = tau .* 1000; % Assuming even time spacing between frames

% tau1_index_CBF = 2;
% tau2_index_CBF = 6;
% tau1_index_CBV = 2;

%%
load('E:\Allen BME-BOAS-27 Data Backup\AZ01 fUS\07-21-2025 awake RC15gV manual right whisker stim\coronal_mask_rep_07_24_2025.mat')

%% Main loop with the Adaptive SVD Thresholding
% for filenum = startFile:endFile
% for filenum = 2:endFile
% for filenum = [285:-1:189]
for filenum = 8
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    IQ = single(squeeze(IData + 1i .* QData));
    clearvars IData QData

    IQm = IQ;
    IQm(coronal_mask_rep) = 0; % Apply the brain mask to the IQ: set the non-brain voxels equal to 0
    
    % Determine the optimal SV thresholds with the spatial similarity matrix
    [xp, yp, zp, nf] = size(IQm);
    PP = reshape(IQm, [xp*yp*zp, nf]);
    tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
    [U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
    disp('Full SVD done')
    toc

    SSM = plotSSM(U, false);
    [~, a_opt, b_opt] = fitSSM(SSM, false); % Get the optimal singular value thresholds
    SVs = diag(S);

%     [PP, EVs, V_sort] = getSVs2D(IQ);
%     disp('SVs decomposed')
    [IQf_opt] = applySVs2D(IQm, PP, SVs, V, a_opt, b_opt);
    a_subopt = a_opt*2;
    b_subopt = b_opt;
    [IQf_subopt] = applySVs2D(IQm, PP, SVs, V, a_subopt, b_opt);

%     volumeViewer(abs(IQf(:, :, :, 1)))
%     figure; imagesc(squeeze(abs(max(IQf(:, :, :, 1), [], 1)))')
    % clearvars IQ

    % Use the IQf with separated negative and positive frequency components
%     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
    
    numg1pts = 20; % Only calculate the first N points
%     g1_n = g1T(IQf_separated{1}, numg1pts);
% %     [CBFsi_n, CBVi_n] = g1_to_CBi(g1_n, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%     g1_p = g1T(IQf_separated{2}, numg1pts);
%     [CBFsi_p, CBVi_p] = g1_to_CBi(g1_p, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
    g1_opt = g1T(IQf_opt, numg1pts);
    g1_subopt = g1T(IQf_subopt, numg1pts);
%     [PDI] = calcPowerDoppler(IQf_separated);
    PDI_opt = sum(abs(IQf_opt) .^ 2, 4) ./ size(IQf_opt, 4);
    PDI_subopt = sum(abs(IQf_subopt) .^ 2, 4) ./ size(IQf_subopt, 4);
end


%% Look at the optimized one
g1_tau1_cutoff = 0.3;
[g1A_mask_opt] = createg1mask(g1_opt, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);
%     [g1A_mask] = createg1mask(g1Avg, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);

[CBFsi_opt, CBVi_opt] = g1_to_CBi(g1_opt, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)

%     CBFsi(~g1A_mask) = -Inf; % Remove noisy points from the CBFspeed index (in theory)
CBFsi_opt(~g1A_mask_opt) = 0; % Remove noisy points from the CBFspeed index (in theory)

%% Suboptimized one
[g1A_mask_subopt] = createg1mask(g1_subopt, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);
%     [g1A_mask] = createg1mask(g1Avg, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);

[CBFsi_subopt, CBVi_subopt] = g1_to_CBi(g1_subopt, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)

%     CBFsi(~g1A_mask) = -Inf; % Remove noisy points from the CBFspeed index (in theory)
CBFsi_subopt(~g1A_mask_subopt) = 0; % Remove noisy points from the CBFspeed index (in theory)

%% erm

SSM = plotSSM(U);

% Test to look at the individual "weighted images"
k_test = 18; % Which column vector to use
ss_wi = reshape(U(:, k_test) * V(:, k_test)', [xp, yp, zp, nf]); % Subspace weighted image
volumeViewer(abs(ss_wi(:, :, :, 1)))
%     figure; imagesc(abs(mean(test, 4)))

%% Compare the optimal and suboptimal thresholded images

figure; imagesc(squeeze(max(CBVi_opt, [], 1))')
figure; imagesc(squeeze(max(CBVi_subopt, [], 1))')

avg_CBVi_opt = mean(CBVi_opt, 'all', 'omitnan')
avg_CBVi_subopt = mean(CBVi_subopt, 'all', 'omitnan')


avg_PDI_opt = mean(PDI_opt, 'all', 'omitnan')
avg_PDI_subopt = mean(PDI_subopt, 'all', 'omitnan')

avg_PDI_opt/avg_PDI_subopt
(b_opt - a_opt + 1)/(b_subopt - a_subopt + 1)

% Average values of the filtered IQ intensities
avg_IQf_opt = squeeze(mean(IQf_opt, [1, 2, 3], 'omitnan'));
avg_IQf_subopt = squeeze(mean(IQf_subopt, [1, 2, 3], 'omitnan'));
figure; plot(1:size(IQf_opt, 4), abs(avg_IQf_opt), '-', 1:size(IQf_subopt, 4), abs(avg_IQf_subopt), '--');
title('Average intensity across the filtered IQ volume, for optimal vs. suboptimal SV thresholds')
legend("Lower threshold = " + num2str(a_opt), "Lower threshold = " + num2str(a_subopt))
xlabel('Frame number (within the superframe)')
ylabel('Intensity [au]')

% Difference of the average IQ intensities
figure; plot(abs(avg_IQf_opt) - abs(avg_IQf_subopt), '-o'); hold on; yline(0); hold off
title('Average intensity across the filtered IQ volume, optimal minus suboptimal')
xlabel('Frame number (within the superframe)')
ylabel('Intensity difference [au]')

%%
noblood_roi = {30:50, 20:30, 80:90};
g1_noblood_roi_opt = squeeze(mean(g1_opt(noblood_roi{1}, noblood_roi{2}, noblood_roi{3}, :), [1, 2, 3]));
g1_noblood_roi_subopt = squeeze(mean(g1_subopt(noblood_roi{1}, noblood_roi{2}, noblood_roi{3}, :), [1, 2, 3]));
figure; plot((1:size(g1_opt, 4)) - 1, abs(g1_noblood_roi_opt), '-o'); hold on; plot((1:size(g1_opt, 4)) - 1, abs(g1_noblood_roi_subopt), '-o'); hold off; legend("Lower threshold = " + num2str(a_opt), "Lower threshold = " + num2str(a_subopt)); xlabel('Tau index'); ylabel('|g1|')

%% Try normalizing the filtered IQ
IQf_opt_n = IQf_opt ./ (b_opt - a_opt);
IQf_subopt_n = IQf_subopt ./ (b_subopt - a_subopt);

% Post-IQ normalization PDI
PDI_opt_n = sum(abs(IQf_opt_n) .^ 2, 4) ./ size(IQf_opt_n, 4);
PDI_subopt_n= sum(abs(IQf_subopt_n) .^ 2, 4) ./ size(IQf_subopt_n, 4);


avg_PDI_opt_n = mean(PDI_opt_n, 'all', 'omitnan')
avg_PDI_subopt_n = mean(PDI_subopt_n, 'all', 'omitnan')

avg_PDI_opt_n/avg_PDI_subopt_n

% Post-IQ normalization tl-fUS
g1_opt_n = g1T(IQf_opt_n, numg1pts);
g1_subopt_n = g1T(IQf_subopt_n, numg1pts);

[g1A_mask_opt_n] = createg1mask(g1_opt_n, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);
[CBFsi_opt_n, CBVi_opt_n] = g1_to_CBi(g1_opt_n, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
CBFsi_opt_n(~g1A_mask_opt_n) = 0; % Remove noisy points from the CBFspeed index (in theory)

[g1A_mask_subopt_n] = createg1mask(g1_subopt_n, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);
[CBFsi_subopt_n, CBVi_subopt_n] = g1_to_CBi(g1_subopt_n, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
CBFsi_subopt_n(~g1A_mask_subopt_n) = 0; % Remove noisy points from the CBFspeed index (in theory)

% Compare the post-IQ-normalization CBVi images
figure; imagesc(squeeze(max(CBVi_opt_n, [], 1))')
figure; imagesc(squeeze(max(CBVi_subopt_n, [], 1))')

figure; imagesc(squeeze(CBVi_opt_n(40, :, :))')
figure; imagesc(squeeze(CBVi_subopt_n(40, :, :))')

(b_opt - a_opt + 1)/(b_subopt - a_subopt + 1)
%% Helper functions

function [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF)

    g1A_T = {};
    
    g1A_T{1} = abs(g1(:, :, :, 2)) > g1_tau1_cutoff; % First treatment: tau1 is above some cutoff (make sure there is some actual blood signal there)
    g1A_T{2} = abs(g1(:, :, :, tau1_index_CBF)) > abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > |g1(tau2)| (noise might have the g1 randomly increase with tau, but it should not happen with a voxel where there is a real blood signal)
    g1A_T{3} = abs(g1(:, :, :, tau1_index_CBF)) > 2 .* abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    % g1A_T{4} = abs(g1(:, :, :, tau1_index_CBF)) - 1 .* abs(g1(:, :, :, tau2_index_CBF)) > tau_difference_cutoff; % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    
    g1A_mask = true(size(g1A_T{1})); % Mask of voxels to keep for the g1 treatments
    for i = 1:length(g1A_T)
        g1A_mask = and(g1A_mask, g1A_T{i});
    end

end