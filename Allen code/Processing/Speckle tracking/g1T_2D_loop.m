%% Only use if needed: convert Jianbo/Bingxue's acquisition parameters to something I can use
% First: manually load an IQ file, like: load('E:\PROJ_tlfUS\IQdata\0806_2021_BL3_vUS_run1(good)\IQ-10-5-5000-1000-1-BL3-1.mat')
if ~exist('P_old', 'var')
    P_old = P; clearvars P
end
[P] = oldP2P(P_old);

%% Add the Speckle tracking folder to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
AllenProcessingCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing\");
addpath(genpath(AllenProcessingCodePath))

%% Load data and things
if ~exist('IQpath', 'var')
    IQpath = uigetdir('D:\Allen\Data\', 'Select the IQ data path');
    IQpath = [IQpath, '\'];
end

if ~exist('savepath', 'var')
    savepath = uigetdir([IQpath, '..\'], 'Select the save path');
    savepath = [savepath, '\'];
end

% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [IQpath, '..\params.mat']);
    load([params_pathname, params_filename])
end

fDim = 3; % Dimension of the data corresponding to frequency (or time)
zDim = 1; % Dimension of the data corresponding to z (axial direction)
xDim = 2; % Dimension of the data corresponding to x (lateral direction)

% Prompt for parameter user input
parameterPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound'};
parameterDefaults = {'1', '', '20', num2str(P.numFramesPerBuffer)};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

startFile = str2double(parameterUserInput{1}); % File to start reconstructing from
endFile = str2double(parameterUserInput{2});   % File to stop reconstructing on
sv_threshold_lower = str2double(parameterUserInput{3});
sv_threshold_upper = str2double(parameterUserInput{4});

%% Set up the High Pass Filter (parameters from the 2020 vUS paper)
HPF.fc = 25; % Cutoff frequency [Hz]
% 25 Hz corresponds to 1 mm/s

HPF.fs = P.frameRate; % Sampling frequency [Hz]
HPF.order = 4; % Butterworth filter order

[HPF.b, HPF.a] = butter(HPF.order, HPF.fc/(HPF.fs/2), 'high');

%% Go through superframes and process
IQfilenameStructure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];
for fi = startFile:endFile
    load([IQpath, IQfilenameStructure, num2str(fi)], 'IQ'); % Load IQ data

    % ========= 1. Preprocessing ========= %%
    % IQ = squeeze(complex(IData, QData));
    % clearvars IData QData
    
    % sv_threshold_lower = 20; sv_threshold_upper = size(IQ, 3);
    
    % Mask the region to actually process
    [zpo, xpo, nfo] = size(IQ); % Original sizes
    % figure; imagesc(squeeze(abs(IQ(:, :, 1))))
    % zrange = 1:100;
    zrange = 1:zpo;
    % zrange = 20:140;
    xrange = 1:xpo;
    % xrange = 40:60;
    IQm = IQ(zrange, xrange, :);
    
    % 1.1 SVD clutter filter
    %     [CM, EVs, V_sort] = getSVs2D(IQ);
    [zp, xp, nf] = size(IQm);
    
    % CM = reshape(IQm, [zp*xp, nf]); % Covariance matrix
    % tic
    % %     [U, S, V] = svd(PP); % Already sorted in decreasing order
    % [U, S, V] = svd(CM, 'econ'); % Already sorted in decreasing order
    % SVs = diag(S);
    % %     disp('Full SVD done')
    % toc
    % disp('SVs decomposed')
    % 
    % [IQf, noise] = applySVs1D(IQm, CM, SVs, V, sv_threshold_lower, sv_threshold_upper);

    [CM, EVs, V_sort] = getSVs1D(IQ);
    disp('SVs decomposed')
    [IQf, noise] = applySVs1D(IQ, CM, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')
    
    % 1.2 High pass filter (apply to the post-SVD clutter filtered data)
    HPF.dim = length(size(IQf)); % Operate on the time dimension
    IQf_HPF = filter(HPF.b, HPF.a, IQf, [], HPF.dim);
    
    % Testing
    % figure; imagesc(squeeze(abs(IQf(:, :, 1))))
    % tempPDI = sum(abs(IQf).^2, 3);
    % figure; imagesc(tempPDI .^ 0.5)
    % tp = [128, 39]; % Test point
    % tp = [87, 26]; % Test point
    % figure; plot(squeeze(abs(IQf_HPF(tp(1), tp(2), :))))
    % figure; plot(squeeze(real(IQf(tp(1), tp(2), :))))
    % figure; plot(squeeze(real(IQf_HPF(tp(1), tp(2), :))))
    
    % ========= 2. Directional flow filtering ========= %%
    
    % 2.1 Separate positive and negative frequencies
    [IQf_separated_noHPF, IQf_FT_separated_noHPF, nFTpts_noHPF] = separatePosNegFreqs(IQf); % Outputs are cell arrays in the order of: negative, positive, all frequencies

    [IQf_separated, IQf_FT_separated, nFTpts] = separatePosNegFreqs(IQf_HPF); % Outputs are cell arrays in the order of: negative, positive, all frequencies
    % frameDim = length(size(IQf)); % Get the dimension corresponding to time/frames
    
    ctp = 1:length(IQf_FT_separated); % Indices of which frequency Components To Process (typically [1, 2, 3]: negative, positive, all)
    ctp_labels = {"Down flows", "Up flows", "All flows"};
    
    % Testing: plot the separated and full Fourier spectrums and reconstructed IQ signals
    faxis = linspace(-P.frameRate/2, P.frameRate/2, nFTpts)';
    % figure; plot(faxis, squeeze(abs(IQf_FT_separated{1}(tp(1), tp(2), :))))
    % figure; plot(faxis, squeeze(abs(IQf_FT_separated{2}(tp(1), tp(2), :))))
    % figure; plot(faxis, squeeze(abs(IQf_FT_separated{3}(tp(1), tp(2), :))))
    % figure; plot(1:P.numFramesPerBuffer, squeeze(abs(IQf_separated{1}(tp(1), tp(2), tp(3), :))))
    % figure; plot(1:P.numFramesPerBuffer, squeeze(abs(IQf_separated{2}(tp(1), tp(2), tp(3), :))))
    % figure; plot(1:P.numFramesPerBuffer, squeeze(abs(IQf_separated{3}(tp(1), tp(2), tp(3), :))))
    
    % 2.2 Mask out all frequencies outside of +/- 1100 Hz, and other bands
    % considered 'system noise' -- see Jianbo's sysNoiseRemove.m
    % freqMask = abs(faxis) > 1100; % [Hz]
    % % *********** add the system noise stuff later ***********
    % 
    % % Create IQf_FT_separated_masked: the Fourier-transformed filtered IQ data,
    % % with some frequencies masked out
    % IQf_FT_separated_masked = cell(size(IQf_FT_separated));
    % IQf_separated_masked = cell(size(IQf_separated));
    % for j = ctp
    %     IQf_FT_separated_masked{j} = IQf_FT_separated{j};
    %     IQf_FT_separated_masked{j}(:, :, freqMask) = 0;
    % 
    %     IQf_separated_masked{j} = ifft(ifftshift(IQf_FT_separated_masked{j}, fDim), nFTpts, fDim);
    % end
    % % figure; plot(faxis, squeeze(abs(IQf_FT_separated{3}(tp(1), tp(2), :))), '-', 'LineWidth', 2)
    % % hold on
    % % plot(faxis, squeeze(abs(IQf_FT_separated_masked{3}(tp(1), tp(2), :))), '--', 'LineWidth', 1)
    % % hold off
    
    % PDI and CDI
    [PDI] = calcPowerDoppler(IQf_separated, noise);
    [CDI] = calcColorDoppler(IQf_FT_separated, P);

    save([savepath, 'PDI_CDI-', num2str(fi), '.mat'], 'PDI', 'CDI', '-v7.3', '-nocompression');
       
    % 3. Calculate g1
    % startTau = 1; % Index for the first tau point (tau1) for subsequent analysis. Changed this from 2 to 1 on 7/8/26 because I changed the g1T.m function to output g1 starting from tau = tau1 instead of tau = 0.
    startTau = 2; % Index for the first tau point (tau1) for subsequent analysis.
    
    nTau = ceil(20e-3 *P.frameRate); % # of time lags to consider; empirically set by assuming all g1 for voxels containing actual flow decay within 10 ms
    % nTau = ceil(100e-3 *P.frameRate); % # of time lags to consider; empirically set by assuming all g1 for voxels containing actual flow decay within 10 ms
    tau = (0:nTau - 1)' ./ P.frameRate; % Time lag vector [s]
    
    % g1neg = g1T(IQf_separated{1}, nTau + startTau - 1); % Add the startTau-1 because the values start at startTau, but we still want nTau points total
    % g1pos = g1T(IQf_separated{2}, nTau + startTau - 1); % Add the startTau-1 because the values start at startTau, but we still want nTau points total
    
    % Store g1 for each frequency component in a cell array
    g1 = cell(size(IQf_separated));
    g1_noHPF = cell(size(IQf_separated_noHPF));
    
    for j = ctp
        g1{j} = g1T(IQf_separated{j}, nTau); % Use the base filtered IQ
        g1_noHPF{j} = g1T(IQf_separated_noHPF{j}, nTau); % Use the base filtered IQ
        % g1{j} = g1T(IQf_separated_masked{j}, nTau); % Use the filtered IQ with system noise removed
    end

    % save([savepath, 'g1-', num2str(fi), '.mat'], 'g1', 'IQf_separated', 'IQf_FT_separated')
    % save([savepath, 'g1-', num2str(fi), '.mat'], 'g1')
    save([savepath, 'g1-', num2str(fi), '.mat'], 'g1', 'g1_noHPF')
    % save([savepath, 'g1-', num2str(fi), '.mat'], 'g1', 'IQf_separated_masked', 'IQf_FT_separated_masked')
end
save([savepath, 'g1_proc_params.mat'], 'nFTpts', 'nTau', 'tau', 'startTau', 'freqMask', 'faxis')
save([savepath, 'fUS_proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper');

%% Testing
% Load a few superframes of the g1 and take the average
g1_sum = cell(3, 1);
si = 1; % Start index
ei = 85; % End index
for fi = si:ei
    load([savepath, 'g1-', num2str(fi), '.mat'], 'g1')

    for j = 1:3
        if isempty(g1_sum{j})
            g1_sum{j} = g1{j};
        else
            g1_sum{j} = g1_sum{j} + g1{j};
        end
    end
end

%%
g1_avg = cell(size(g1_sum));
nfie = ei - si + 1; % # of files in ensemble
for j = 1:3
    g1_avg{j} = g1_sum{j} ./ nfie;
end

pixelTimeseriesGUI(g1_avg{3}, squeeze(abs(g1_avg{3}(:, :, 2))), 'ComplexMode', 'abs')