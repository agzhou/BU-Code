%% Description:
%       3D fUS processing: loop through IQ data files and calculate PDI
%       with XDoppler

clearvars

%% 1. Load acquisition parameters, timing info, beamforming parameters, etc.
IQpath = uigetdir('G:\', 'Select the IQ data path');
IQpath = [IQpath, '\'];

codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
AllenProcessingCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing");
addpath(AllenProcessingCodePath)

% Load parameters
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
% Create vectors defining x and z coordinates in actual length units
% x_mm = (1:PData.Size(2)) .* PData.PDelta(1) .* P.wl .* 1e3; % x [mm]
% z_mm = (1:PData.Size(1)) .* PData.PDelta(3) .* P.wl .* 1e3; % z [mm]

IQfilenameStructure = ['IQ-RCCR-', num2str(round(P.maxAngle)), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

savepath = uigetdir([IQpath, '\..'], 'Select the save path');
savepath = [savepath, '\'];

% addpath([cd, '\Speckle tracking']) % add path for the g1 calculation functions

% Load or define the timing data
% TDchoices = {'Yes', 'No - Manually Define'};
% TDopts.Default = TDchoices{1}; TDopts.Interpreter = 'none';
% TDanswer = questdlg('Use automatically-acquired stim timing data?', 'Timing data source', TDchoices{1}, TDchoices{2}, TDopts);
% switch TDanswer
%     case TDchoices{1}
%         [timingFilePathFN, timingFilePath] = uigetfile([IQpath, '..\Timing data\TD.mat'], 'Select the timing data');
%         timingFilePath = [timingFilePath, timingFilePathFN];
%         load(timingFilePath)
%     case TDchoices{2}
%         manualTimingPrompt = {'Baseline duration [s]', 'Stim duration [s]', 'Rest duration [s]', 'Offset duration [s]', 'Number of trials'};
%         manualTimingDefaults = {'5', '5', '20', '0', '10'};
%         manualTimingUserInput = inputdlg(manualTimingPrompt, 'Input Parameters', 1, manualTimingDefaults);
% 
%         % Store the user inputs for stim timing parameters into the corresponding variables
%         baseline_duration = str2double(manualTimingUserInput{1});
%         stim_duration = str2double(manualTimingUserInput{2});
%         rest_duration = str2double(manualTimingUserInput{3});
%         offset_duration = str2double(manualTimingUserInput{4});
%         num_trials = str2double(manualTimingUserInput{5});
%         stim_sample_rate = 1000; % [Hz]
% 
%         % trial_stim_pattern = 
% 
% end

% load(timingFilePath, 'acqStart', 'airPuffOutput', 'daqStartTimetag', 'sfTimeTags', 'sfTimeTagsDAQStart', 'sfTimeTagsDAQStart_adj', 'sfWidth', 'sfWidth_adj', 'timeStamp')

%% 2. Define some processing parameters

procParamsPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound', 'z start', 'z end'};
procParamsDefaults = {'1', '', '30', num2str(P.numFramesPerBuffer), '20', '130'};
procParamsUserInput = inputdlg(procParamsPrompt, 'Input Parameters', 1, procParamsDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(procParamsUserInput{1});
endFile = str2double(procParamsUserInput{2});
numFiles = endFile - startFile + 1;
sv_threshold_lower = str2double(procParamsUserInput{3});
sv_threshold_upper = str2double(procParamsUserInput{4});
zstart = str2double(procParamsUserInput{5});
zend = str2double(procParamsUserInput{6});

clearvars procParamsPrompt procParamsDefaults procParamsUserInput

% %% Define time lags
% % startTau = t1i; % Index for the first tau point (tau1) for subsequent analysis.
% 
% nTau = ceil(20e-3 *P.frameRate); % # of time lags to consider; empirically set by assuming all g1 for voxels containing actual flow decay within 10 ms
% % nTau = ceil(100e-3 *P.frameRate); % # of time lags to consider; empirically set by assuming all g1 for voxels containing actual flow decay within 10 ms
% tau = (0:nTau - 1)' ./ P.frameRate; % Time lag vector [s]

%% Set up the High Pass Filter (parameters from the 2020 vUS paper)
HPF.fc = 25; % Cutoff frequency [Hz]
% 25 Hz corresponds to 1 mm/s

HPF.fs = P.frameRate; % Sampling frequency [Hz]
HPF.order = 4; % Butterworth filter order

[HPF.b, HPF.a] = butter(HPF.order, HPF.fc/(HPF.fs/2), 'high');

%% Loop through files: Clutter Filtering and Power Doppler (XDoppler) calculation
for filenum = startFile:endFile
% for filenum = endFile:-1:startFile
% for filenum = 2:endFile
% for filenum = startFile
% for filenum = 1
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])

    % Mask the region to actually process
    [xpo, ypo, zpo, nfo] = size(IQ_CR); % Original sizes
    % figure; imagesc(squeeze(max(abs(IQ_CR(:, :, :, 1)), [], 1)).')
    
    xrange = 1:xpo;
    yrange = 1:ypo;
    % zrange = 1:100;
    % zrange = 1:zpo;
    zrange = zstart:zend;
    voxelRange = {xrange, yrange, zrange};
    xp = length(xrange); yp = length(yrange); zp = length(zrange); nf = size(IQ_CR, 4);
    
    % Mask volumes
    IQ_CR = IQ_CR(xrange, yrange, zrange, :);
    IQ_RC = IQ_RC(xrange, yrange, zrange, :);
    
    % [brainMask, maskInfo] = createBrainVoxelMask(sum(abs(IQ), 4), 'MarginVoxels', 0);
    % IQ = IQ .* brainMask;   % zero out non-brain voxels before the SVD sees them
    
    % SVD decluttering - apply on both RC and CR
    IQf_RC = SVD_CF(IQ_RC, sv_threshold_lower, sv_threshold_upper);
    IQf_CR = SVD_CF(IQ_CR, sv_threshold_lower, sv_threshold_upper);
    
    % IQf_OPW = IQf_RC + IQf_CR;
    % PDI_OPW = mean(abs(IQf_OPW).^2, 4);
    
    % XD = IQf_RC .* conj(IQf_CR);
    PDI = abs(mean(IQf_RC .* conj(IQf_CR), 4)); % Not exactly the estimator from the paper, but close
    % PDI_XD = 2.*real(mean(XD, 4)); % This version is sensitive to phase shifts from flow, between RC and CR acquisitions
    % PDI = abs(mean(IQf_RC .* conj(IQf_CR) + IQf_CR .* conj(IQf_RC), 4)); % This version is sensitive to phase shifts from flow, between RC and CR acquisitions

    % figure; imagesc(squeeze(max(PDI, [], 1)).' .^ 0.5); colormap hot

    % save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'CDI', 'g1', 'noise', '-v7.3', '-nocompression');
    % save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'CDI', 'g1', '-v7.3', '-nocompression');
    save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', '-v7.3', '-nocompression');

    disp("fUS data for file " + num2str(filenum) + " saved" )
%     disp("g1 result for file " + num2str(filenum) + " saved" )

    toc

    clearvars PDI IQ_RC IQ_CR IQf_RC IQf_CR
    
end
% savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'tau1_index_CBF', 'tau2_index_CBF', 'tau1_index_CBV');
save([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'voxelRange', 'HPF', 'xp', 'yp', "zp", "nf");
% savefast([savepath, 'PDI_CDI_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');

%% Loop through files: Clutter Filtering and Power Doppler (OPW) calculation
mkdir([savepath, 'OPW'])
% for filenum = startFile:endFile
% for filenum = endFile:-1:startFile
% for filenum = 2:endFile
% for filenum = startFile
for filenum = 1
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])

    % OPW: coherently sum CR and RC
    IQ = IQ_CR + IQ_RC;

    % Mask the region to actually process
    [xpo, ypo, zpo, nfo] = size(IQ); % Original sizes
    % figure; imagesc(squeeze(max(abs(IQ(:, :, :, 1)), [], 1)).')
    
    xrange = 1:xpo;
    yrange = 1:ypo;
    % zrange = 1:100;
    % zrange = 1:zpo;
    zrange = zstart:zend;
    voxelRange = {xrange, yrange, zrange};
    xp = length(xrange); yp = length(yrange); zp = length(zrange); nf = size(IQ, 4);
    
    % Mask volumes
    IQ = IQ(xrange, yrange, zrange, :);
    
    % [brainMask, maskInfo] = createBrainVoxelMask(sum(abs(IQ), 4), 'MarginVoxels', 0);
    % IQ = IQ .* brainMask;   % zero out non-brain voxels before the SVD sees them
    
    % SVD decluttering - apply on both RC and CR
    IQf = SVD_CF(IQ, sv_threshold_lower, sv_threshold_upper);
    
    % IQf_OPW = IQf_RC + IQf_CR;
    % PDI_OPW = mean(abs(IQf_OPW).^2, 4);
    
    PDI = squeeze(mean(abs(IQf).^2, 4));
    % figure; imagesc(squeeze(max(PDI, [], 1)).' .^ 0.5); colormap hot

    % save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'CDI', 'g1', 'noise', '-v7.3', '-nocompression');
    % save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'CDI', 'g1', '-v7.3', '-nocompression');
    save([savepath, 'OPW\fUSdata-', num2str(filenum), '.mat'], 'PDI', '-v7.3', '-nocompression');

    disp("fUS data for file " + num2str(filenum) + " saved" )
%     disp("g1 result for file " + num2str(filenum) + " saved" )

    toc

    % clearvars PDI IQ_RC IQ_CR IQ IQf
    
end
% savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'tau1_index_CBF', 'tau2_index_CBF', 'tau1_index_CBV');
save([savepath, 'OPW\fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'voxelRange', 'HPF', 'xp', 'yp', "zp", "nf");
% savefast([savepath, 'PDI_CDI_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');
