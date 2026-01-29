%% Description:
%       3D fUS and functional connectivity (FC) processing
%       Timing data should be processed and saved with plotfUStiming_FC.m first

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
% for filenum = 3

    % Load the IQ data
    tic
    load([IQpath, IQfilenameStructure, num2str(filenum)])
    
    IQ = single(squeeze(IData + 1i .* QData));
%     clearvars IData QData

    % figure; imagesc(squeeze(max(abs(IQ(:, :, :, 2)), [], 1))')

    % Calculate the cross correlation of raw IQ (masked) to look at motion
    ixc = calcIXC(IQ);
    figure; plot((1:size(IQ, 4)) ./ P.numFramesPerBuffer .* 1e3, abs(ixc)); xlabel('Micro time [ms]'); ylabel('|Cross correlation of images|')
    % figure; plot(abs(ixc)); xlabel('Frame'); ylabel('|Cross correlation of images|')
    title("Superframe " + num2str(filenum))

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

