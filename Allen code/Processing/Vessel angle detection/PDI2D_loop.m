%% Description:
%       2D (tl-)fUS processing
%       Timing data should be processed with plotfUStiming.m first

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
x_mm = (1:PData.Size(2)) .* PData.PDelta(1) .* P.wl .* 1e3; % x [mm]
z_mm = (1:PData.Size(1)) .* PData.PDelta(3) .* P.wl .* 1e3; % z [mm]

IQfilenameStructure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

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

procParamsPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound'};
procParamsDefaults = {'1', '', '20', num2str(P.numFramesPerBuffer)};
procParamsUserInput = inputdlg(procParamsPrompt, 'Input Parameters', 1, procParamsDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(procParamsUserInput{1});
endFile = str2double(procParamsUserInput{2});
numFiles = endFile - startFile + 1;
sv_threshold_lower = str2double(procParamsUserInput{3});
sv_threshold_upper = str2double(procParamsUserInput{4});

clearvars procParamsPrompt procParamsDefaults procParamsUserInput

%% 3. Loop through files: Clutter Filtering and Power Doppler calculation
for filenum = startFile:endFile
% for filenum = 2:endFile
% for filenum = startFile
% for filenum = 1
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    
    % IQ = squeeze(IData + 1i .* QData);
    % clearvars IData QData
    
    % SVD decluttering
%     [xp, yp, zp, nf] = size(IQ);
    
    [PP, EVs, V_sort] = getSVs1D(IQ);
    disp('SVs decomposed')
    [IQf, noise] = applySVs1D(IQ, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')

%     figure; imagesc(squeeze(abs(IQf(:, :, 1))) .^ 0.5)

    % clearvars IQ

    % Use the IQf with separated negative and positive frequency components
    [IQf_separated, IQf_FT_separated, nFTpts] = separatePosNegFreqs(IQf); % Outputs are cell arrays in the order of: negative, positive, all frequencies
    [PDI] = calcPowerDoppler(IQf_separated, noise);
    [CDI] = calcColorDoppler(IQf_FT_separated, P);

    % PDI = sum(abs(IQf) .^ 2, 3) ./ size(IQf, 3);
    % PDI = sum(abs(IQf) .^ 2, 3) ./ size(IQf, 3) ./ noise;
    % figure; imagesc(x_mm, z_mm, squeeze(PDI .^ 0.5)); colormap hot; colorbar; title('Power Doppler'); xlabel('x [mm]'); ylabel('z [mm]')
    % figure; imagesc(x_mm, z_mm, squeeze(abs(IQ(:, :, 1)))); colorbar; title('IQ'); xlabel('x [mm]'); ylabel('z [mm]')

%     save([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI', '-v7.3', '-nocompression');
%     disp("PDI and CDI for file " + num2str(filenum) + " saved" )
    % save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', '-v7.3', '-nocompression');
    save([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'CDI', '-v7.3', '-nocompression');

    disp("fUS result for file " + num2str(filenum) + " saved" )
%     disp("g1 result for file " + num2str(filenum) + " saved" )

    toc
    
end
% savefast([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'tau', 'tau_ms', 'tau1_index_CBF', 'tau2_index_CBF', 'tau1_index_CBV');
save([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');
% savefast([savepath, 'PDI_CDI_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');

%% 4. Store all the PDI across the experiment into one matrix
% load([savepath, 'PDI_CDI-', num2str(1), '.mat'], 'PDI', 'CDI')
% load([savepath, 'fUSdata-', num2str(1), '.mat'], 'PDI', 'CDI')
load([savepath, 'fUSdata-', num2str(2), '.mat'], 'PDI')
% PDIallSF = cell([length(PDI), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
% PDIallSF = cell([size(PDI)]); 
PDIallSF = zeros([size(PDI), endFile - startFile + 1]); % Matrix with the CBVi for every superframe
% PDIallSF(:,  1) = PDI;
% PDIallSF(:, :, 1) = PDI;
% CDIallSF = cell([size(CDI)]); % Matrix with the CBVi for every superframe
% CDIallSF(:,  1) = CDI;

% for filenum = startFile + 1:endFile
for filenum = startFile:endFile
%     load([savepath, 'PDI_CDI-', num2str(filenum), '.mat'], 'PDI', 'CDI')
%     load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI', 'CDI')
    load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI')
%     PDI = load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'PDI')
%     CDI = load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'CDI')

%     for i = 1:3
%         PDIallSF{i} = cat(3, PDIallSF{i}, PDI{i});
%         CDIallSF{i} = cat(3, CDIallSF{i}, CDI{i});
%     end

    if iscell(PDI)
        PDIallSF(:, :, filenum) = PDI{3};
    else
        PDIallSF(:, :, filenum) = PDI;
    end
end

