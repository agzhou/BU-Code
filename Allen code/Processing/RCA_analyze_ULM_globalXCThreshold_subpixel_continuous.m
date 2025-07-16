%% Description
% ULM data analysis:
% Take reconstructed data from individual buffers/batches/superframes, each with
% some number of frames/subframes, each containing N acquisitions/angles.

% Re-sample with a rolling method to get more effective frames.

% SVD to separate the bubble signals from the tissue signal and other
% clutter

% ...

% **** Make sure that the current directory is set to the directory where this
% script is located!!! ****

%% Use parallel processing for speed

% https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox

pp = gcp('nocreate');
if isempty(pp)
    % There is no parallel pool
    parpool LocalProfile1

end

%% Load parameters and make folder for saving the processed data

% Get data path of the reconstructed IQ data
datapath = uigetdir('G:\Allen\Data\', 'Select the IQ data path');
datapath = [datapath, '\'];

% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [datapath, '..\params.mat']);
    load([params_pathname, params_filename])
end

% Load Verasonics reconstruction parameters: datapath\PData.mat
if ~exist('PData', 'var')
    load([datapath, 'PData.mat'])
end

% Prompt for parameter user input
parameterPrompt = {'Start file number', 'End file number', 'SVD lower bound', 'SVD upper bound', 'Image refinement factor - x', 'Image refinement factor - y', 'Image refinement factor - z', 'XC Threshold Factor', 'x pixel spacing [um]', 'y pixel spacing [um]', 'z pixel spacing [um]', 'Ensemble size'};
parameterDefaults = {'1', '', '20', '150', '2', '2', '2', '0.2', num2str(PData.PDelta(1) * P.wl * 1e6), num2str(PData.PDelta(2) * P.wl * 1e6), num2str(PData.PDelta(3) * P.wl * 1e6), '5000'};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% define # of files manually for now
% str2double(parameterUserInput{});
startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});
numFiles = endFile - startFile + 1;
sv_threshold_lower = str2double(parameterUserInput{3});
sv_threshold_upper = str2double(parameterUserInput{4});
imgRefinementFactor = [str2double(parameterUserInput{5}), str2double(parameterUserInput{6}), str2double(parameterUserInput{7})];
if any(floor(imgRefinementFactor) ~= imgRefinementFactor) || any(imgRefinementFactor < 1)
    error('Image refinement factors must be whole numbers')
end
XCThreshold = str2double(parameterUserInput{8});
xpix_spacing = str2double(parameterUserInput{9});
ypix_spacing = str2double(parameterUserInput{10});
zpix_spacing = str2double(parameterUserInput{11});
es = str2double(parameterUserInput{12}); % Ensemble size (for SVD)

% clearvars parameterPrompt parameterDefaults parameterUserInput

savepath = uigetdir([datapath, '..\'], 'Select the save path');
savepath = [savepath, '\'];

filename_structure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

addpath([cd, '\normxcorr3.m'])
addpath([cd, '\Parthasarathy_Radial_Particle_Localization'])

%% Other parameters for processing the data

% Region of interest
% xrange = int16(1:80);
% yrange = int16(1:80);
% zrange = int16(1:142);

% framerange = 1:200;
% framerange = 1:size(IQf, 3);
% range = {xrange, yrange, zrange, framerange};
% range = {xrange, yrange, zrange};

% Load and refine simulated PSF
if ~exist('PSF', 'var')
%     load('G:\Allen\Data\RC15gV PSF sim\PSF.mat', 'PSF')
    % figure; imagesc(squeeze(abs(PSF(40, :, :)))')
    datapath_split = split(string(datapath), filesep);
    PSF_path = fullfile(join(datapath_split(1:find(contains(datapath_split, 'Data'))), '\') + "\RC15gV PSF sim\PSF.mat");
    load(PSF_path, 'PSF')
end

% PSFs = PSF(190:210, 118:138, :); % PSF section
PSFs = PSF(30:50, 30:50, 92:110); % PSF section
refPSF = imresize3(PSFs, [size(PSFs, 1) * imgRefinementFactor(1), size(PSFs, 2) * imgRefinementFactor(2), size(PSFs, 3) * imgRefinementFactor(3)]);
% volumeViewer(abs(refPSF))

%% Process the data

% Calculate the total # of frames in the experiment
totalNumFrames = numFiles * P.numFramesPerBuffer;

% # of ensembles to process, if we don't overlap
nes = floor(totalNumFrames / es);
%%%%%%%%%%%% Add a check for the ensemble size being an integer multiple of the # of
% frames per sf %%%%%%%%%%%%%%%%%%%%%%%%%%
nsfpes = es/P.numFramesPerBuffer; % # of superframes per ensemble

%% Go through each ensemble and process the data
for en = 1:nes % en = ensemble number
% for en = 1
    tic
    % tic
    IQen = []; % IQ stacked over the ensemble
    
    for sfi = 1:nsfpes % Go through all the superframes for one ensemble (sfi = superframe index)
        load([datapath, filename_structure, num2str((en - 1) * nsfpes + sfi), '.mat'])  % load each reconstructed buffer/batch/superframe
        
        % IQ = squeeze(IData + 1i .* QData);   % Combine I and Q, which are saved separately. It's easier to save the big reconstructed data with savefast, which doesn't support complex values. The data is already a coherent sum.
        % clear IData QData

        IQen = cat(4, IQen, squeeze(IData + 1i .* QData)); % Add the next superframe's IQ to the ensemble IQ
        disp("File " + num2str((en - 1) * nsfpes + sfi) + " loaded")
    end
    clearvars IQ
    
    % toc

%     if filenum == 1
        [xp, yp, zp, nf] = size(IQen);
        xrange = int16(1:xp);
        yrange = int16(1:yp);
        zrange = int16(1:zp);
        range = {xrange, yrange, zrange};
        range{4} = int16(1:nf); % set frame range after rolling on the first file
%     end
    zpixfactor = zpix_spacing ./ xpix_spacing; % relative z pixel spacing vs x or y, for the radial centers algorithm

    % SVD proc part 1
    tic
    [PP, EVs, V_sort] = getSVs2D(IQen);
    disp('SVs decomposed')
    toc

    plot_FFT_SVs_function(V_sort, P)
    figure; plot(abs(log10(EVs))); title('Singular value magnitude'); xlabel('Singular value number'); ylabel('log10(Singular value magnitude)')

    % SVD proc part 2
%     tic
    [IQenf] = applySVs2D(IQen, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')

    % generateTiffStack_acrossframes(abs(IQenf) .^ 0.5, [8.8, 8.8, 8], 'gray', 1:80)
    % figure; imagesc(squeeze(max(abs(IQenf(30:50, :, :, 10))))')
    % volumeViewer(abs(IQenf(:, :, :, 10)))
    % generateIQVideo(IQenf, [8.8, 8.8, 8], 10)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % LOOKING AT THE FULL SVD
%     IQen_single = single(IQen);
% 
%     % Determine the optimal SV thresholds with the spatial similarity matrix
%     [xp, yp, zp, nf] = size(IQen_single);
%     PP = reshape(IQen_single, [xp*yp*zp, nf]);
%     tic
% %     [U, S, V] = svd(PP); % Already sorted in decreasing order
%     [U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
%     disp('Full SVD done')
%     toc
% 
%     %
%     SSM = zeros(nf, nf); % Initialize the spatial similarity matrix
% 
%     SSM_const = 1/(xp * yp * zp); % constant in front of the summation term
% %
%     tic
%     for n = 1:nf
% %     for n = 1:10
%         abs_u_n = abs(U(:, n)); % The nth column vector from U
%         mean_abs_u_n = sum(abs_u_n) / length(abs_u_n);
%         stddev_abs_u_n = std(abs_u_n);
% %         for m = 1:nf
%         for m = 1:n % leverage the symmetry of the SSM
%             abs_u_m = abs(U(:, m)); % The mth column vector from U
%             mean_abs_u_m = sum(abs_u_m) / length(abs_u_m);
%             SSM(n, m) = sum( ((abs_u_n - mean_abs_u_n) .* (abs_u_m - mean_abs_u_m)) ...
%                         ./ stddev_abs_u_n ...
%                         ./ std(abs_u_m) );
%         end
%     end
%     SSM = SSM .* SSM_const; % Normalize
%     SSM = SSM + SSM'; % Apply the symmetry to fill out the "missing" values. 
%     % NOTE: I think this adds an extra 1 to all diagonal values
%     toc
%     figure; imagesc(SSM); axis square % Show the SSM
% 
% 
%     % Test to look at the individual "weighted images"
%     k_test = 20; % Which column vector to use
%     test = reshape(U(:, k_test) * V(:, k_test)', [xp, yp, zp, nf]);
%     volumeViewer(abs(mean(test, 4)))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %% XC test
    % rfnX = imgRefinementFactor(1); % refinement pixel increase factor
    % rfnY = imgRefinementFactor(2);
    % rfnZ = imgRefinementFactor(3);
    % 
    % if ~all(imgRefinementFactor == 1)
    %     % refIQs = zeros(size(IQs, 1) * rfnX, size(IQs, 2) * rfnY, size(IQs, 3) * rfnZ, size(IQs, 4)); % refined IQ section
    % 
    %     % go through all frames and refine
    %     % parfor f = 1:size(IQs, 4)
    %     for f = 1
    %         I_temp = IQs(:, :, :, f);
    %         refIQs(:, :, :, f) = imresize3(I_temp, [size(I_temp, 1) * rfnX, size(I_temp, 2) * rfnY, size(I_temp, 3) * rfnZ], 'cubic');
    %     end
    % else
    %     refIQs = IQs;
    % end
    % disp('IQ upsampled')
    % %% cross correlation
    % 
    % % normxcorr3 from file exchange
    % % (https://www.mathworks.com/matlabcentral/fileexchange/73946-normxcorr3-fast-3d-ncc)
    % XC = normxcorr3(abs(refPSF), abs(refIQs(:, :, :, 1)), 'same'); % Cross correlate the filtered/refined images and the simulated PSF
    % for f = 2:size(refIQs, 4)
    %     XC(:, :, :, f) = normxcorr3(abs(refPSF), abs(refIQs(:, :, :, f)), 'same');
    % end
    % disp('Upsampled IQ cross correlated with the PSF')

    % clearvars PP EVs V_sort

%     IQd = diff(IQ, 1, 4); % Frame subtraction

    % Could chunk it
    chunk_size = P.numFramesPerBuffer;
    range{4} = 1:chunk_size;
    for ci = 1:es/chunk_size % chunk index
    % for ci = 1
        % [centersRC, refIQs, XC] = localizeBubbles3D_globalXCThreshold_subpixel_noparfor(IQenf(:, :, :, (ci - 1) * chunk_size + 1 : ci * chunk_size), refPSF, range, imgRefinementFactor, XCThreshold, zpixfactor);
        [centersRC, refIQs, XC] = localizeBubbles3D_globalXCThreshold_subpixel(IQenf(:, :, :, (ci - 1) * chunk_size + 1 : ci * chunk_size), refPSF, range, imgRefinementFactor, XCThreshold, zpixfactor);
    
        savefast([savepath, 'centers-', num2str((en - 1) * nsfpes + ci)], 'centersRC')
    
        disp(strcat("Center finding done: chunk ", num2str((en - 1) * nsfpes + ci)))
    end
    toc
end
img_size = [xp, yp, zp] .* imgRefinementFactor;
save([savepath, 'proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'PSF', 'range', 'imgRefinementFactor', 'XCThreshold', 'xpix_spacing', 'ypix_spacing', 'zpix_spacing', 'img_size', 'es', 'chunk_size', 'nsfpes')

%% Helper functions

function [fileCount] = countFiles(fileName,filePath)

    fileInfo = strsplit(fileName, '-');

    % This keeps everything except the last numeric component
    prefix = strjoin(fileInfo(1:end-1), '-'); 
    
    % Construct the search pattern
    searchPattern = fullfile(filePath, prefix + "-*.mat"); % Wildcard for different numbers
    
    % Get a list of matching files
    fileList = dir(searchPattern);
    
    % Count the number of matching files
    fileCount = numel(fileList);


end

function [] = generateIQVideo(IQen, actualSize, frameRate)
    IQen = abs(IQen);
    % Set up the figure for the video
    savepath = uigetdir('D:\Allen\Data\', 'Select the save path for the video');
    savepath = [savepath, '\'];

    vo = VideoWriter([savepath, 'volumeVideo']);

    vo.Quality = 100;
    vo.FrameRate = frameRate;
    open(vo);

%     V = volshow(IQen, 'RenderingStyle', 'MaximumIntensityProjection');
%     viewer = V.Parent;
%     % V_old = V;
% %     V.Alphamap(1:100) = 0;          % Change transparency
%     viewer.BackgroundColor = [1, 1, 1];  % Make background white
%     viewer.BackgroundGradient = 'off'; % Turn off the background gradient
%     % V.ScaleFactors(3) = size(IQen, 1) / size(IQen, 3) * actualSize(1) / actualSize(3); % scale with # pixels and region size
% %     V.CameraPosition = V.CameraPosition ./ 2;
% %     V.CameraPosition = [2.1161 -3.7332 -0.1764];
%     viewer.CameraUpVector = [0 0 -1];
% %     V.CameraViewAngle = 15;

    plotPower = 1;
    
    tic
    for fi = 1:size(IQen, 4) % Go through each frame

        %%%%%%%%%%%%%%%%
        if mod(fi - 1, 10) == 0 % only get a video frame every N frames
            disp(fi)
            if fi == 1
                V = volshow(IQen(:, :, :, fi) .^ plotPower, 'RenderingStyle', 'MaximumIntensityProjection');
                viewer = V.Parent;
            else
                volshow(IQen(:, :, :, fi) .^ plotPower, 'RenderingStyle', 'MaximumIntensityProjection', Parent=viewer);
            end
%             V.Alphamap(1:100) = 0;          % Change transparency
            viewer.BackgroundColor = [1, 1, 1];  % Make background white
            viewer.BackgroundGradient = 'off'; % Turn off the background gradient

            viewer.CameraUpVector = [0, 0, -1]; % Flip the z axis
            viewer.CameraZoom = 1.2; % Zoom
            % V.CameraViewAngle = 15;
            % V.ScaleFactors(3) = size(densityMapInterpolated, 1) / size(densityMapInterpolated, 3) * actualSize(1) / actualSize(3); % scale with # pixels and region size
            % V.ScaleFactors = V.ScaleFactors .* 1.5; % zoom

            cv = getframe(viewer.Parent);     % get the current volume
            rgb = frame2im(cv);      % convert the frame to rgb data
            writeVideo(vo, rgb);
        end
    end

    close(vo)
    toc

end
