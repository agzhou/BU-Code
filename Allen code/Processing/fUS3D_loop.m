%% Description:
%       3D (tl-)fUS processing
%       Timing data should be processed with plotfUStiming.m first

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
[timingFilePathFN, timingFilePath] = uigetfile([IQpath, '..\Timing data\TD.mat'], 'Select the timing data');
timingFilePath = [timingFilePath, timingFilePathFN];
load(timingFilePath)
% load(timingFilePath, 'acqStart', 'airPuffOutput', 'daqStartTimetag', 'sfTimeTags', 'sfTimeTagsDAQStart', 'sfTimeTagsDAQStart_adj', 'sfWidth', 'sfWidth_adj', 'timeStamp')
%% Define some parameters

parameterPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound', 'Tau 1 index for CBFspeed', 'Tau 2 index for CBFspeed', 'Tau 1 index for CBV'};
parameterDefaults = {'1', '', '20', '', '2', '11', '2'};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});
numFiles = endFile - startFile + 1;
sv_threshold_lower = str2double(parameterUserInput{3});
sv_threshold_upper = str2double(parameterUserInput{4});
tau1_index_CBF = str2double(parameterUserInput{5});
tau2_index_CBF = str2double(parameterUserInput{6});
tau1_index_CBV = str2double(parameterUserInput{7});

clearvars parameterPrompt parameterDefaults parameterUserInput

taustep = 1/P.frameRate;
% tau = taustep:taustep:(P.numFramesPerBuffer * taustep);
tau = 0:taustep:((P.numFramesPerBuffer - 1) * taustep);
tau_ms = tau .* 1000; % Assuming even time spacing between frames

% tau1_index_CBF = 2;
% tau2_index_CBF = 6;
% tau1_index_CBV = 2;

%% Set up the High Pass Filter
% fc = 50; % Cutoff frequency [Hz]
% fs = P.frameRate; % Sampling frequency [Hz]
% HPF_order = 3; % Butterworth filter order
% 
% [HPF_b, HPF_a] = butter(HPF_order, fc/(fs/2), 'high');

%% Define the mask manually for now

% load('E:\Allen BME-BOAS-27 Data Backup\AZ01 fUS\07-21-2025 awake RC15gV manual right whisker stim\coronal_mask_rep_07_24_2025.mat')
% load('I:\Ultrasound Data from 04-11-2025 to 05-08-2025\05-06-2025 AZ03 fUS pre-stroke\run 1 all frames stacked\coronal_mask_rep_07_31_2025.mat')
% load('J:\Ultrasound data from 7-21-2025\08-06-2025 AZ01 RCA fUS\coronal_mask_rep.mat')

%% Define other cropping
%     zstart = 40;
% %     zstart = 50;
%     zend = size(IQ, 3);
%     zstart = 15;
%     zstart = 45;
    zstart = 52;
%     zend = 105;
    zend = 127;

%% Save proc params
numg1pts = 20; % Only calculate the first N points
save([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'numg1pts', 'zstart', 'zend');

%% Main loop
for filenum = startFile:endFile
% for filenum = [2:endFile]
% for filenum = [endFile - 1:-1:startFile]
% for filenum = [116:endFile]
% for filenum = 1

    % Load the IQ data
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    
    IQ = single(squeeze(IData + 1i .* QData));
    clearvars IData QData

    % figure; imagesc(squeeze(max(abs(IQ(:, :, :, 2)), [], 1))')
    
    % Crop the IQ first 

    IQm = IQ(:, :, zstart:zend, :);

%     figure; imagesc(squeeze(max(abs(IQm(:, :, :, 2)), [], 1))')

    %%%%%%%%%%%%%% IF USING THE PREDEFINED MASK %%%%%%%%%%%%
%     IQm = IQ;
%     IQm(coronal_mask_rep) = 0; % Apply the brain mask to the IQ: set the non-brain voxels equal to 0

    % Apply the HPF
%     dim = length(size(IQm)); % Operate on the time dimension
%     IQm_HPF = filter(HPF_b, HPF_a, IQm, [], dim);
    IQm_HPF = IQm;

    % SVD decluttering
%     [PP, EVs, V_sort] = getSVs2D(IQm);
    [xp, yp, zp, nf] = size(IQm_HPF);
    PP = reshape(IQm_HPF, [xp*yp*zp, nf]);
    tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
    [U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
    SVs = diag(S);
%     disp('Full SVD done')
    toc
    disp('SVs decomposed')

    [IQf_HPF, noise] = applySVs2D(IQm_HPF, PP, SVs, V, sv_threshold_lower, sv_threshold_upper);
%     [IQf, noise] = applySVs2D(IQm, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')

%     volumeViewer(abs(IQf(:, :, :, 1)))
%     figure; imagesc(squeeze(abs(max(IQf(:, :, :, 1), [], 1)))')
    % clearvars IQ

    % Use the IQf with separated negative and positive frequency components
%     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
    
%     g1_n = g1T(IQf_separated{1}, numg1pts);
% %     [CBFsi_n, CBVi_n] = g1_to_CBi(g1_n, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%     g1_p = g1T(IQf_separated{2}, numg1pts);
%     [CBFsi_p, CBVi_p] = g1_to_CBi(g1_p, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
    g1 = g1T(IQf_HPF, numg1pts);
%     g1 = g1T(IQf);
%     [CBFsi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
% %     savefast([savepath, 'fUSdata-', num2str(filenum), '.mat'], g1, CBFi, CBVi);

%     [PDI] = calcPowerDoppler(IQf_separated);
    PDI = sum(abs(IQf_HPF) .^ 2, 4) ./ size(IQf_HPF, 4);
%     [CDI] = calcColorDoppler(IQf_FT_separated, P);

%     figure; imagesc(squeeze(max(PDI, [], 1))' .^ 0.5); colormap hot
%     figure; imagesc(squeeze(max(PDI ./ noise, [], 1))' .^ 0.5); colormap hot
%     volumeViewer(PDI)

%     save([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI', '-v7.3', '-nocompression');
%     disp("PDI and CDI for file " + num2str(filenum) + " saved" )
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'CBFsi', 'CBVi', 'PDI', 'CDI', '-v7.3', '-nocompression');
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'CBFsi', 'CBVi', 'PDI', 'CDI', 'g1_n', 'g1_p', 'CBFsi_n', 'CBVi_n', 'CBFsi_p', 'CBVi_p',  '-v7.3', '-nocompression');
%     save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'g1_n', 'g1_p', 'PDI', 'CDI', '-v7.3', '-nocompression');
    save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1', 'PDI', 'noise', '-v7.3', '-nocompression');
%     save([savepath, 'g1-', num2str(filenum), '.mat'], 'g1', 'g1_n', 'g1_p', '-v7.3', '-nocompression');
%     save([savepath, 'g1-', num2str(filenum), '.mat'], 'g1', '-v7.3', '-nocompression');

    disp("fUS result for file " + num2str(filenum) + " saved" )
%     disp("g1 result for file " + num2str(filenum) + " saved" )

    toc
    
end
% savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'tau1_index_CBF', 'tau2_index_CBF', 'tau1_index_CBV');
% savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'numg1pts');
% save([savepath, 'fUS_proc_params.mat'], 'HPF_a', 'HPF_b', 'tau', 'tau_ms', 'numg1pts', 'zstart', 'zend');
% save([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'numg1pts', 'zstart', 'zend');
% save([savepath, 'fUS_proc_params.mat'], 'HPF_a', 'HPF_b', 'tau', 'tau_ms', 'numg1pts', 'coronal_mask_rep');
% save([savepath, 'fUS_proc_params.mat'], 'tau', 'tau_ms', 'numg1pts', 'coronal_mask_rep');
% savefast([savepath, 'PDI_CDI_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');

%% Testing
testIQf = applySVs2D(IQ, PP, EVs, V_sort, 35, sv_threshold_upper);
figure; imagesc(squeeze(max(abs(testIQf(30:50, :, :, 1)), [], 1) .^ 0.5)'); colormap hot
% volumeViewer(abs(testIQf(:, :, :, 1)))
testPDI = mean(abs(IQf) .^ 2, 4);
figure; imagesc(squeeze(max(testPDI(30:50, :, :, 1), [], 1) .^ 0.5)'); colormap hot
figure; imagesc(squeeze(max(testPDI(:, :, :, 1), [], 3) .^ 0.5)'); colormap hot

%% Convert g1 into CBV, CBFspeed, etc.

g1_tau1_cutoff = 0.2;
% g1_tau1_cutoff = 0.1;

% g1_tau1_cutoff = 0.0;
% tau_difference_cutoff = 0.2;

for filenum = startFile:endFile
% for filenum = [endFile]
% for filenum = 1
%     load([savepath, 'g1-', num2str(filenum)], 'g1') % Load the saved g1 mat files
    load([savepath, 'fUSdata-', num2str(filenum)], 'g1') % Load the saved g1 mat files

    [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);
%     [g1A_mask] = createg1mask(g1Avg, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);

    [CBFsi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)

%     CBFsi(~g1A_mask) = -Inf; % Remove noisy points from the CBFspeed index (in theory)
    CBFsi(~g1A_mask) = 0; % Remove noisy points from the CBFspeed index (in theory)

    save([savepath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
%     save([savepath, 'tlfUSdatatest-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
    disp("tl-fUS result for file " + num2str(filenum) + " saved" )

end
save([savepath, 'tlfUS_proc_params.mat'], 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff', 'g1A_mask');
% save([savepath, 'tlfUStest_proc_params.mat'], 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff');
figure; imagesc(squeeze(max(CBVi(:, :, :), [], 1) .^ 0.3)'); colormap hot
figure; imagesc(squeeze(max(CBVi(:, :, :), [], 3) .^ 0.5)'); colormap hot
vcmap = colormap_ULM;
figure; imagesc(squeeze(mean(CBFsi(:, :, :), 1))'); colormap(vcmap)

% generateTiffStack_multi({CBVi .^ 0.7}, [8.8, 8.8, 8], 'hot', 5)
%% Get and save PDI, CDI only
for filenum = startFile:endFile
% for filenum = 1
% for filenum = [285:-1:189]
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
%     load(['E:\Allen BME-BOAS-27 Data Backup\AZ03 Stroke RC15gV\fUS\05-06-2025 pre-stroke\IQ Data - Verasonics recon\', IQfilenameStructure, num2str(filenum)])
    
    IQ = squeeze(IData + 1i .* QData);
    clearvars IData QData
    
    % SVD decluttering
    [xp, yp, zp, nf] = size(IQ);
    
    [PP, EVs, V_sort] = getSVs2D(IQ);
    disp('SVs decomposed')
    [IQf] = applySVs2D(IQ, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')

    clearvars IQ

    % Use the IQf with separated negative and positive frequency components
%     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
    
    [PDI] = calcPowerDoppler(IQf);
%     [CDI] = calcColorDoppler(IQf_FT_separated, P);

%     save([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI', '-v7.3', '-nocompression');
    save([savepath, 'PDI-', num2str(filenum), '.mat'], 'PDI', '-v7.3', '-nocompression');
    disp("PDI and CDI for file " + num2str(filenum) + " saved" )

    toc
    
end
savefast([savepath, 'PDI_CDI_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');


%% CBFspeed index with the Derivative Method
% vtest = sqrt( (abs(g1(:, :, :, tau1_index_CBF)) - abs(g1(:, :, :, tau2_index_CBF))) ./ (abs(g1(:, :, :, tau1_index_CBF)) .* ( tau(tau2_index_CBF)^2 - tau(tau1_index_CBF)^2 )) );
CBFspeed_DM = (abs(g1(:, :, :, tau1_index_CBF)) - abs(g1(:, :, :, tau2_index_CBF))) ./ (abs(g1(:, :, :, tau1_index_CBF)) .* ( tau(tau2_index_CBF)^2 - tau(tau1_index_CBF)^2 ));
CBFspeed_DM(CBFspeed_DM < 0) = 0;
CBFspeed_DM = sqrt(CBFspeed_DM);

CBFspeed_DM_masked = CBFspeed_DM;
CBFspeed_DM_masked(~g1A_mask) = 0;

% vtest_NLMF = NLMF(vtest ./ max(vtest, [], 'all')); % Try using nonlinear means filtering to get rid of those high-noise points
% vtest_gaussLPF = imgaussfilt3(vtest);

%% Compare the derivative CBFspeed index method to the old log one
CBFsi_all_masked = CBFsi_all;
CBFsi_all_masked(~g1A_mask) = 0;

figure; imagesc(squeeze(max(CBFsi_all_masked(30:50, :, :), [], 1))'); colormap(vcmap)
figure; imagesc(squeeze(mean(CBFsi_all_masked(30:50, :, :), 1))'); colormap(vcmap)

% generateTiffStack_MeanIPs_multi({CBFsi_all_masked .^ 1}, [8.8, 8.8, 8], vcmap, 10)
% generateTiffStack_multi({CBFsi_all_masked .^ 1}, [8.8, 8.8, 8], vcmap, 10)

generateTiffStack_MeanIPs_multi({CBFsi_all_masked .^ 1, CBFspeed_DM_masked}, [8.8, 8.8, 8], vcmap, 20)


%% Trying some alternate CBFspeed index calculations
% CBFsi_all = squeeze(abs(g1(:, :, :, tau1_index_CBF)) - abs(g1(:, :, :, tau2_index_CBF)));
% CBFsi_p = squeeze(abs(g1_p(:, :, :, tau1_index_CBF)) - abs(g1_p(:, :, :, tau2_index_CBF)));
% CBFsi_n = squeeze(abs(g1_n(:, :, :, tau1_index_CBF)) - abs(g1_n(:, :, :, tau2_index_CBF)));
% volumeViewer(CBFsi_all)
% figure; imagesc(squeeze(max(CBFsi_all(30:50, :, :), [], 1))')

vcmap = colormap_ULM; % velocity colormap
% generateTiffStack_multi({CBFsi_all}, [8.8, 8.8, 8], vcmap, 5)
% generateTiffStack_multi({CBFsi_p}, [8.8, 8.8, 8], vcmap, 5)
% generateTiffStack_multi({CBFsi_n}, [8.8, 8.8, 8], vcmap, 5)

% CBFsi_all_test = squeeze(abs(g1(:, :, :, tau1_index_CBF)) ./ abs(g1(:, :, :, tau2_index_CBF))); %./ (tau2_index_CBF - tau1_index_CBF);
% CBFsi_all_test = squeeze(log(abs(g1(:, :, :, tau1_index_CBF)) ./ abs(g1(:, :, :, tau2_index_CBF)))); %./ (tau2_index_CBF - tau1_index_CBF);
% CBFsi_all_test = squeeze(abs(g1(:, :, :, tau1_index_CBF)) ./ abs(g1(:, :, :, tau2_index_CBF))); %./ (tau2_index_CBF - tau1_index_CBF);
CBFsi_all_test = squeeze(log(abs(g1(:, :, :, tau1_index_CBF))) - log(abs(g1(:, :, :, tau2_index_CBF)))); %./ (tau2_index_CBF - tau1_index_CBF);
figure; imagesc(squeeze(max(CBFsi_all_test(30:50, :, :), [], 1))'); colormap(vcmap)

figure; imagesc(squeeze(mean(CBFsi_all_test(30:50, :, :), 1))'); colormap(vcmap)

volumeViewer(CBFsi_all_test)

%% Store all the updated CBVi and CBFsi across the experiment into one matrix
load([savepath, 'tlfUSdata-', num2str(1), '.mat'], 'CBFsi', 'CBVi')
CBViallSF = zeros([size(CBVi), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
CBViallSF(:, :, :, 1) = CBVi;

CBFsiallSF = zeros([size(CBFsi), endFile - startFile + 1]); % Matrix with the CBFsi for every superframe
CBFsiallSF(:, :, :, 1) = CBFsi;

for filenum = startFile + 1:endFile
    load([savepath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi')
    CBViallSF(:, :, :, filenum) = CBVi;
    CBFsiallSF(:, :, :, filenum) = CBFsi;
end

%% Store all the PDI across the experiment into one matrix (with separated frequencies)
% % load([savepath, 'PDI_CDI-', num2str(1), '.mat'], 'PDI', 'CDI')
% load([savepath, 'fUSdata-', num2str(1), '.mat'], 'PDI', 'CDI')
% % PDIallSF = cell([length(PDI), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
% PDIallSF = cell([size(PDI)]); % Matrix with the CBVi for every superframe
% % PDIallSF(:,  1) = PDI;
% CDIallSF = cell([size(CDI)]); % Matrix with the CBVi for every superframe
% % CDIallSF(:,  1) = CDI;
% 
% % for filenum = startFile + 1:endFile
% for filenum = startFile:endFile
% %     load([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI')
%     load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'CDI')
% %     PDI = load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI')
% %     CDI = load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'CDI')
% 
%     for i = 1:3
% %         PDIallSF1 = cat(4, PDIallSF1, PDI{1});
% %         CDIallSF1 = cat(4, CDIallSF1, CDI{1});
% %         PDIallSF2 = cat(4, PDIallSF2, PDI{2});
% %         CDIallSF2 = cat(4, CDIallSF2, CDI{2});
% %         PDIallSF3 = cat(4, PDIallSF3, PDI{3});
% %         CDIallSF3 = cat(4, CDIallSF3, CDI{3});
%         PDIallSF{i} = cat(4, PDIallSF{i}, PDI{i});
%         CDIallSF{i} = cat(4, CDIallSF{i}, CDI{i});
%     end
% end

%% Store all the PDI across the experiment into one matrix
% load([savepath, 'PDI_CDI-', num2str(1), '.mat'], 'PDI', 'CDI')
% load([savepath, 'fUSdata-', num2str(1), '.mat'], 'PDI', 'CDI')
load([savepath, 'fUSdata-', num2str(1), '.mat'], 'PDI')
PDIallSF = zeros([size(PDI), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
PDIallSF(:, :, :, 1) = PDI;

% for filenum = startFile + 1:endFile
for filenum = startFile:endFile
%     load([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI')
    load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI')

    PDIallSF(:, :, :, filenum) = PDI;
end

%% Visualize the PDI and CDI across the experiment
% mwr = 30:50; % MIP window range
% mdim = 1; % MIP dimension
% newsize = size(CBViallSF);
% newsize(mdim) = 1; % Set the size of the new variable to 1
% CBViMIPStack = zeros(newsize);
PDIMIPStack = squeeze(max(PDIallSF{1}(30:50, :, :, :), [], 1));
CDIMIPStack = squeeze(max(CDIallSF{1}(30:50, :, :, :), [], 1));

%
generateTiffStack_acrossframes(PDIallSF{3} .^ 0.7, [8.8, 8.8, 8], 'hot', 1:80)
generateTiffStack_acrossframes(PDIallSF .^ 0.5, [8.8, 8.8, 8], 'hot', 1:80)

%% Visualize the CBVi across the experiment
% mwr = 30:50; % MIP window range
% mdim = 1; % MIP dimension
% newsize = size(CBViallSF);
% newsize(mdim) = 1; % Set the size of the new variable to 1
% CBViMIPStack = zeros(newsize);
CBViMIPStack = squeeze(max(CBViallSF(30:50, :, :, :), [], 1));
%% Check different MIPs across superframes
yr = 20:40;
generateTiffStack_acrossframes(CBViallSF .^ 0.7, [8.8, 8.8, 8], 'hot', yr)
% generateTiffStack_acrossframes(CBViallSF .^ 1, [8.8, 8.8, 8], 'hot', yr)


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

%% Prepare template(s) for atlas registration
% Create templates for each hemodynamic parameter, averaging across superframes
CBVi_allSF_avg = mean(CBViallSF, 4);
CBFsi_allSF_avg = mean(CBFsiallSF, 4);
PDI_allSF_avg = mean(PDIallSF, 4);

voxel_size = PData.PDelta .* P.wl; % Voxel size (y, x, z) in meters
fUS_volume_dimensions_m = [P.Trans.numelements/2 * P.Trans.spacingMm / 1e3, P.Trans.numelements/2 * P.Trans.spacingMm / 1e3, (P.endDepthMM - P.startDepthMM)/1e3]; % Volume size in meters
fUS_volume_dimensions_voxels = PData.Size; % Volume size in voxels (from the recon PData)

% Adjust the sizes based on the pre-SVD/clutter filtering cropping
fUS_cropped_volume_dimensions_voxels = size(CBVi_allSF_avg);
fUS_cropped_volume_dimensions_m = fUS_cropped_volume_dimensions_voxels ./ fUS_volume_dimensions_voxels .* fUS_volume_dimensions_m;

targetVoxelSizePrompt = {'y Target Voxel Size [um]', 'x Target Voxel Size [um]', 'z Target Voxel Size [um]'};
targetVoxelSizeDefaults = {'10', '10', '10'};
targetVoxelSizeUserInput = inputdlg(targetVoxelSizePrompt, 'Input Target Voxel Size', 1, targetVoxelSizeDefaults);

% Store target voxel size inputs and convert to meters
target_voxel_size(1) = str2double(targetVoxelSizeUserInput{1}) ./ 1e6;
target_voxel_size(2) = str2double(targetVoxelSizeUserInput{2}) ./ 1e6;
target_voxel_size(3) = str2double(targetVoxelSizeUserInput{3}) ./ 1e6;

prereg_interp_factor = voxel_size ./ target_voxel_size;

% Resample hemodynamic parameter template maps to the desired voxel size
CBVi_allSF_avg_rs = imresize3(CBVi_allSF_avg, 'Scale', prereg_interp_factor, 'Method', 'cubic');
CBFsi_allSF_avg_rs = imresize3(CBFsi_allSF_avg, 'Scale', prereg_interp_factor, 'Method', 'cubic');
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

%% Look at a ROI (rPDI thresholded)
numPtsUSI = P.Mcr_fcp.apis.seq_length_s * P.daqrate / interp_factor; % # of time points per trial for the upsampling
% Calculate the timecourse from the average within that ROI
roi_rPDI_TA = zeros(size(rPDI_TA, 3), 1);
% repmat(roi_mask, [1, 1, stim_pattern.trial_duration])
for ti = 1:numPtsUSI
% for ti = 1
     temp_rPDI_TA = rPDI_TA(:, :, :, ti);
%      temp_roi_rPDI_avg = mean(temp_rPDI_TA(roi_mask));
     temp_roi_rPDI_avg = mean(temp_rPDI_TA(am_rPDI_t_roi_mask));
     roi_rPDI_TA(ti) = temp_roi_rPDI_avg;
end

% Plot the average timecourse in the ROI
figure; plot((1:length(roi_rPDI_TA)) .* interp_factor ./ P.daqrate, roi_rPDI_TA); xlabel('Time [s]'); ylabel('rPDI'); title("rPDI ROI timecourse")
mmws = 30; % Movmean window size (in units of the trial interpolation rate)
figure; plot((1:length(roi_rPDI_TA)) .* interp_factor ./ P.daqrate, smoothdata(roi_rPDI_TA, 'movmean', mmws)); xlabel('Time [s]'); ylabel('rPDI'); title("rPDI ROI timecourse, moving mean over " + num2str(mmws/length(roi_rPDI_TA) * P.Mcr_fcp.apis.seq_length_s) + "s")

%% Plot activation at each slice
for slice = 1:10
    my_inds = (slice-1)*5:slice*5;
    my_inds = my_inds+1;
    figure; imagesc(squeeze(mean(test(my_inds, :, :), 1))')
    title(num2str(mean(my_inds)))
end

kernel = ones(3, 3, 3);
kernel(2, 2, 2) = 3;%sum(kernel, 'all');

test_r_CBVi_relative_change_conv = convn(test_r_CBVi_relative_change, kernel, 'same');

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

