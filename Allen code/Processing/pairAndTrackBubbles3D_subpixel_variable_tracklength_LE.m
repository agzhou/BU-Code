% Description: Part 2 of the ULM processing for a 2D array. This script
%   takes the result of the bubble localization from 'RCA_analyze_ULM_globalXCThreshold_subpixel.m' ('centers-filenum.mat') and
%   performs pairing and tracking, plus further image refinement for
%   plotting purposes.

%%%%% Long Ensemble version (loading) %%%%%%

%   centersRC is a a cell array of dimensions (# frames per file, 1) with each cell being an (N x 3)
%   matrix. N corresponds to the number of bubbles detected, and each column corresponds to the
%   y, x, and z coordinates where the bubbles were detected. The
%   coordinates will be fractional from the subpixel method.
%   There is one 'centers' for each file that the
%   localization was performed on.

%   This script stacks all the coordinates across every frame in the measurement,
%   then does pairing with the Hungarian method
%   (assignmunkres). 
% 
%   Tracks are created with a VARIABLE frame persistence condition (only keep a
%   track if the same bubble is paired across p frames). Velocities are calculated on this data
%   according to the change in bubble position and the time elapsed between frames. 
%   The tracks are then smoothed and estimated with a Kalman filter, then further refined
%   with acceleration and direction constraints (throw away a track if
%   any acceleration or any direction change exceeds some threshold). The
%   velocity maps are plotted by interpolating between points in a track,
%   setting each interpolated point's velocity equal to that of the
%   original value between the two points. Pixels with multiple tracks will
%   take the average velocity value across all the tracks that intersect
%   it.   

% Dependencies:
%   parpool (Parallel Computing Toolbox)
%   assignmunkres (Sensor Fusion and Tracking Toolbox)
%   trimmean (Statistics and Machine Learning Toolbox)
%   ULM_interp3D_linear.m
%   colormap_ULM
%   smoothn.m (From Damien Garcia: https://www.mathworks.com/matlabcentral/fileexchange/25634-smoothn/)
%   densityMap3D.m      
% Acknowledgement: using Jianbo Tang's ULM code, the Song group's ULM papers, and Jean-Yves Tinevez's simpletracker as references
clearvars

%% 0. Use parallel processing for speed
% https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox

pp = gcp('nocreate');
if isempty(pp)
    % There is no parallel pool
    parpool LocalProfile1

end

%% 1. Add dependencies and load parameters
% Get data path of the localized bubble centers
% datapath = uigetdir('F:\Allen\Data\', 'Select the data path');
datapath = uigetdir('D:\Allen\Data\', 'Select the processed data path');
datapath = [datapath, '\'];

% Load localization processing parameters: proc_params.mat
% load([datapath, 'proc_params.mat'])

% Choose and load the params.mat file (from the acquisition)
[params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [datapath, '..\params.mat']);
load([params_pathname, params_filename])

% Add Jianbo's functions to the Matlab path for the Colormaps_fUS function
try  % Attempt to auto-find the directory in the current path (as Allen's Github is structured)
    oldcodePath_split = split(string(cd), filesep);
    oldcodePath = fullfile(join(oldcodePath_split(1:find(contains(oldcodePath_split, "BU-Code"))), "\") + "\Previous lab code\A-US-ULM\SubFunctions");
    addpath(oldcodePath)
catch % If that doesn't work, manually choose the folder
    oldcodePath = uigetdir('C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Previous lab code\A-US-ULM\SubFunctions\', 'Select the old ULM SubFunctions path');
    oldcodePath = [oldcodePath, '\'];
    addpath(oldcodePath)
end

% % [add line to load the recon pixel spacing .mat file]
% % For now, hard coding and assuming equal x, y, and z spacing
% pix_spacing = P.wl/2;
% [PData_filename, PData_pathname, ~] = uigetfile('*.mat', 'Select the PData file', [datapath, '..\PData.mat']);
% load([PData_pathname, PData_filename])

addpath([cd, '\munkres'])

% Define some region sizes for adjusting the aspect ratio of figures
lateral_width = P.Trans.elementLength * P.wl;
axial_depth = (P.endDepthMM - P.startDepthMM) / 1e3;

%% Look for all center files in the data folder
datapath_dir = dir(datapath);
% center_filenames = datapath_dir.name;
center_filenames = {};
for cfi = 1:length(datapath_dir)% go through center file indices
    center_filename_cfi = datapath_dir(cfi).name;
    % Keep the filename if it's a centers_.mat file
    if contains(center_filename_cfi, 'centers') & contains(center_filename_cfi, '.mat')
        center_filenames{end + 1} = center_filename_cfi;
    end
end
clearvars cfi center_filename_cfi



%% 1.5 Prompt for parameter user input
parameterPrompt = {'Start file number', 'End file number', 'x pixel spacing [um]', 'y pixel spacing [um]', 'z pixel spacing [um]', 'Maximum expected flow speed [mm/s]', 'Persistence frames', 'Moving window size [frames]', 'Acceleration constraint factor', 'Trimmed mean percentage', 'Direction constraint'};
%%% NEED TO CHANGE THE PIX SPACING TO USE THE saved PData %%%
half_pi = pi/2;
parameterDefaults = {'', '', num2str(P.Trans.spacingMm * 1e3), num2str(P.Trans.spacingMm * 1e3), num2str(P.wl/2 * 1e6), '50', '5', '3', '2', '20', num2str(half_pi)};
% parameterDefaults = {'', '', num2str(xpix_spacing * 1e6), num2str(ypix_spacing * 1e6), num2str(zpix_spacing * 1e6), '50', '3', '3', '2', '20', num2str(half_pi)};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

% Store the user inputs for parameters into the corresponding variables
startFile = str2double(parameterUserInput{1});
endFile = str2double(parameterUserInput{2});
xpix_spacing = str2double(parameterUserInput{3}) / 1e6;
ypix_spacing = str2double(parameterUserInput{4}) / 1e6;
zpix_spacing = str2double(parameterUserInput{5}) / 1e6;
maxSpeedExpectedMMPerS = str2double(parameterUserInput{6});           % max expected flow speed [mm/s]
pers = str2double(parameterUserInput{7}); % # of frames a track needs to persist through to keep it
mmws = str2double(parameterUserInput{8}); % Moving mean window size [frames]
aThresholdFactor = str2double(parameterUserInput{9});         % Acceleration change threshold factor
vTrimmedMeanPercentage = str2double(parameterUserInput{10});   % Trimmed mean percentage for the acceleration change threshold [%]
angleChangeThreshold = str2double(parameterUserInput{11});    % Angle change threshold [radians]

%% 2. Store some parameters and create a new grid for the ULM
numFiles = endFile - startFile + 1;  % # of files (superframes/buffers) to process
% totalFrames = numFiles * P.numFramesPerBuffer; % Total frames to process

% Cell array with an entry for each frame. Each entry contains (# bubbles) of coordinate pairs (z, x) of the detected bubble centers
% centerCoords = cell(totalFrames, 1); 
centerCoords = {};

% Concatenate all the centers- files
for n = startFile:endFile   % Go through each center file (for each buffer)
% for n = startFile
    tic
    load([datapath, center_filenames{n}])
    for ind = 1:length(centersALL)

        centerCoords = [centerCoords; centersALL{ind}'];
    end
    disp(strcat("Center coordinates for file ", num2str(n), " stored."))
    toc
end
totalFrames = length(centerCoords);

old_img_size = img_size; % Save the old upsampled image size before we
clearvars img_size
% create a new grid

clear centersTemp indTemp xc yc zc n bfi tsl

% Define the grid/voxel size
gridPrompt = {'x voxel size [um]', 'y voxel size [um]', 'z voxel size [um]'};
gridDefault = {'10', '10', '10'};
gridInput = inputdlg(gridPrompt, 'Define the ULM grid size', 1, gridDefault);
xvs = str2double(gridInput{1}) / 1e6; % x voxel size
yvs = str2double(gridInput{2}) / 1e6; % y voxel size
zvs = str2double(gridInput{3}) / 1e6; % z voxel size

% ovs = [xpix_spacing, ypix_spacing, zpix_spacing] ./ imgRefinementFactor; % old voxel sizes

nxv = lateral_width / xvs; % # of x voxels
nyv = lateral_width / yvs; % # of x voxels
nzv = axial_depth / zvs; % # of x voxels
vs_factor = [nxv, nyv, nzv] ./ old_img_size; % Factor to multiply the old coordinates with to interpolate to the new grid

% Map the coordinates to the new grid
centerCoords_newgrid = cell(size(centerCoords));
for f = 1:length(centerCoords)
    centerCoords_newgrid{f} = round(centerCoords{f} .* repmat(vs_factor, size(centerCoords{f}, 1), 1));
end

%% 3. Calculate bubble count and correct for frames that have bubbles at every voxel
bubbleCount = zeros(length(centerCoords_newgrid), 1);   % numFiles/# buffers x # frames per buffer. Count of bubbles in each frame
centerCoords_corrected = centerCoords_newgrid;          % Correct the centerCoords because some frames have every pixel identified as        a bubble

parfor fi = 1:length(centerCoords_newgrid) % frame index - go through every frame
    bufTemp = centerCoords_newgrid{fi};

    % Correct for some error that makes every pixel a bubble
%     if size(bufTemp, 1) >= img_size(1) * img_size(2) * img_size(3)
%         centerCoords_corrected{fi} = NaN;
%         bubbleCount(fi) = 0;
%     else
        bubbleCount(fi) = size(bufTemp, 1);
%     end
end

totalCount = sum(bubbleCount, 'all');

% Plot the bubble count
% figure; plot(1:length(bubbleCount), bubbleCount, '.')
figure; plot(1:length(bubbleCount), bubbleCount)
title('Bubble count')
xlabel('Frame number')
ylabel('Bubble count')

clear fi bufTemp

%% 3.5 Plot the raw bubble density map
img_size = [nyv, nxv, nzv];
numPadVoxels = 40;
pad_dims = [numPadVoxels, numPadVoxels, numPadVoxels]; % # of voxels to pad with in each dimension
bubbleDensityMapRaw = padarray(zeros(img_size(1), img_size(2), img_size(3)), pad_dims, 0, ['both']);
img_size = size(bubbleDensityMapRaw);

% Add the padding to the centerCoords to avoid negative values
centerCoords_corrected = cell(size(centerCoords_newgrid));
for f = 1:length(centerCoords_newgrid)
    centerCoords_corrected{f} = centerCoords_newgrid{f} + pad_dims;
end

for cci = 1:length(centerCoords_corrected) % centerCoords index
    cc = centerCoords_corrected{cci};
%     cc = cc + pad_dims;
%     cc = round(cc); %%%%%%%%% FOR TESTING %%%%%%%%%
    for nbcci = 1:size(cc, 1) % # bubbles in centerCoords_corrected at index cci
        bubbleDensityMapRaw(cc(nbcci, 1), cc(nbcci, 2), cc(nbcci, 3)) = bubbleDensityMapRaw(cc(nbcci, 1), cc(nbcci, 2), cc(nbcci, 3)) + 1;
    end
end
% 
% % volumeViewer(bubbleDensityMapRaw .^ 0.5)
% % figure; imagesc(squeeze(sum(bubbleDensityMapRaw, 1))' .^ 0.5); colormap hot; title('Raw bubble density, sum across y'); colorbar
figure; imagesc(squeeze(max(bubbleDensityMapRaw, [], 1))' .^ 0.5); colormap hot; title('Raw bubble density, MIP across y'); colorbar
% % figure; imagesc(squeeze(sum(bubbleDensityMapRaw(70:90, :, :), 1))' .^ 0.25); colormap hot; title('Raw bubble density, MIP across y = 70:90 \^0.25'); colorbar


%% 4. Define max speed (distance per frame) threshold and initialize variables
startFramePrompt = {'Start frame'}; % User input for the frame to start processing at
startFrameDefault = {'1'};
startFrameUserInput = inputdlg(startFramePrompt, 'Choose start frame #', 1, startFrameDefault);
startFrame = str2double(startFrameUserInput{1});                      % Frame to start processing at
timePerFrame = 1 / P.frameRate;                                       % time elapsed per frame [s]
maxDistPerFrameM = (maxSpeedExpectedMMPerS / 1000) * timePerFrame;    % max distance traveled per frame [m], according to the max expected flow speed and frame rate
% pixelsPerM = 1 / pix_spacing * imgRefinementFactor(1);              % # of pixels per meter, which depends on the pixel spacing from reconstruction and the image refinement factor from the localization
% maxPixelDistPerFrame = maxDistPerFrameM * pixelsPerM;               % max distance traveled per frame in units of pixels
xpixelsPerM = 1 / (xvs);      % # of x pixels per meter
ypixelsPerM = 1 / (yvs);      % # of y pixels per meter
zpixelsPerM = 1 / (zvs);      % # of z pixels per meter
maxXPixelDistPerFrame = maxDistPerFrameM * xpixelsPerM;               % max x distance traveled per frame in units of pixels
maxYPixelDistPerFrame = maxDistPerFrameM * ypixelsPerM;               % max y distance traveled per frame in units of pixels
maxZPixelDistPerFrame = maxDistPerFrameM * zpixelsPerM;               % max z distance traveled per frame in units of pixels
maxPixelDistPerFrame = [maxXPixelDistPerFrame, maxYPixelDistPerFrame, maxZPixelDistPerFrame];

%% 5. Calculate pairing
% ccc_copy = [centerCoords_corrected(1:end-1), centerCoords_corrected(2:end)]; % Copy of the centerCoords_corrected with the next file in the same row so it doesn't need to pass the whole cell array to each worker in the parpool
ccc_copy_source = [centerCoords_corrected(1:end-1)]; % Copy of the centerCoords_corrected with the next file in the same row so it doesn't need to pass the whole cell array to each worker in the parpool
ccc_copy_target = [centerCoords_corrected(2:end)];
bubblePairs = cell(totalFrames - 1, 1);   % Initialize cell vector of paired bubble indices

ubS = cell(totalFrames - 1, 1);             % unassigned bubbles from the source frames
ubT = cell(totalFrames - 1, 1);             % unassigned bubbles from the target frames

tic
parfor f = startFrame:totalFrames - 1       % Go through frames
% parfor f = startFrame:startFrame+100
% for f = startFrame:startFrame+0
    sourceFrame = ccc_copy_source{f};     % Get the coordinates for the source frame (f)
    targetFrame = ccc_copy_target{f}; % Get the coordinates for the target frame (f + 1)

    nbS = size(sourceFrame, 1);             % number of bubbles in the source frame
    nbT = size(targetFrame, 1);             % number of bubbles in the target frame
    D = NaN(nbS, nbT);                      % Initialize the distance matrix comparing source and target points' distances

    if nbS > 1 && nbT > 1                   % Only go through the pairing if the source and target frames have at least one bubble
        % Go through each source frame point and get the distance to each
        % point in the target frame. Store in the distance matrix D.
        for spi = 1:nbS % source point index
            sourcePoint = sourceFrame(spi, :);      % x, y, and z coords of the source frame's point "spi"
            d = targetFrame - sourcePoint;          % vectorized difference between the x, y, and z coords of all the points in the target frame and point spi from the source frame
            d(:, 1) = d(:, 1) ./ xpixelsPerM;       % Convert the x distance differences into natural distance units [m]
            d(:, 2) = d(:, 2) ./ ypixelsPerM;       % Convert the y distance differences into natural distance units [m]
            d(:, 3) = d(:, 3) ./ zpixelsPerM;       % Convert the z distance differences into natural distance units [m]
            D(spi, :) = sqrt(sum((d .^ 2), 2));     % distance formula on the above
    
        end
        
%         D(D > maxDistPerFrameM) = Inf;      % Set the elements above the distance per frame threshold to Inf so they aren't considered for pairing
%         [assignment, unassignedrows, unassignedcolumns] = assignmunkres(D, 10000000000000); % Pair with the Munkres algorithm, which minimizes the total cost (total paired distance)
%         bubblePairs{f} = assignment; % assignment will always be sorted to make the second column in order
%         ubS{f} = unassignedrows;
%         ubT{f} = unassignedcolumns;

        % Use Yi Cao's munkres for efficiency, which might break if there
        % are Infs?
        [assignment, cost] = munkres(D);
        indTemp = find(assignment);
        [sourceInd, targetInd] = ind2sub(size(assignment), indTemp);
        bubblePairs{f} = [sourceInd, targetInd]; % assignment will always be sorted to make the second column in order

    end
end
toc
disp('Pairing done')

clear f nbS nbT D spi assignment unassignedrows unassignedcolumns ccc_copy_source ccc_copy_target

%% 5.5 Plot the paired bubble count
bubbleCountPaired = zeros(size(bubblePairs, 1), 1);   % numFiles/# buffers x # frames per buffer. Count of bubbles in each frame
parfor bci = 1:length(bubbleCountPaired)
    bubbleCountPaired(bci) = size(bubblePairs{bci}, 1);
end
figure
plot(1:length(bubbleCount), bubbleCount, 1:length(bubbleCountPaired), bubbleCountPaired)
title('Bubble count: raw vs. paired')
legend('Raw', 'Paired')

%% 6. Create tracks with the lower bound of variable persistence

tic
% Separate the pairs of coordinates so we can change their sizes independently
bubblePairsPers = cell(length(bubblePairs), 2); % bubble pairs with the pairs separated into another cell dimension
for n = startFrame:length(bubblePairs)
    if ~isempty(bubblePairs{n})
        bubblePairsPers{n, 1} = bubblePairs{n}(:, 1); % Store the info for the source frame
        bubblePairsPers{n, 2} = bubblePairs{n}(:, 2); % Store the info for the target frame
    end
end

% Store the bubble indices that contribute to each track of at least (pers) # of frames
tracks = cell(size(bubblePairsPers, 1) - pers, 1); % Initialize the tracks cell array

for n = 1:length(bubblePairsPers) - pers                        % Go through each frame (and the following pers # of frames)
    bubblePairsPersTemp = bubblePairsPers;                      % Create a temporary version of the bubble pairs so we can take the info and edit it for the current track (I could make it more efficient by only taking a section)
    for pfc = 1:pers - 1                                        % persistence frame count
        startIndex = bubblePairsPersTemp{n + pfc - 1, 2};       % indices of the paired "target" bubbles in frame n + 1 (pair n), which will be sorted in ascending order
        endIndex = bubblePairsPersTemp{n + pfc, 1};             % indices of the paired "source" bubbles in frame n + 2 (pair n + 1), not necessarily in ascending order because it's aligned with the "target" indices in frame n + 1
        % find the common values in the start and end vectors, and the corresponding indices for each
        [trackContinuesIndices, is, ie] = intersect(startIndex, endIndex, 'stable'); 
        
        % Update the paired and tracked list for this immediate set of pairs
        bubblePairsPersTemp{n + pfc - 1, 1} = bubblePairsPersTemp{n + pfc - 1, 1}(is);
        bubblePairsPersTemp{n + pfc - 1, 2} = bubblePairsPersTemp{n + pfc - 1, 2}(is);
        bubblePairsPersTemp{n + pfc, 1} = bubblePairsPersTemp{n + pfc, 1}(ie);
        bubblePairsPersTemp{n + pfc, 2} = bubblePairsPersTemp{n + pfc, 2}(ie);
        
    end
    % Go back and update the previous pairs too
    for recpfc = 1:pfc - 1   % recursive persistence frame count
        recStartIndex = bubblePairsPersTemp{n + pfc - recpfc - 1, 2};    % indices of the paired "target" bubbles in frame n + 1 (pair n), which will be sorted in ascending order
        recEndIndex = bubblePairsPersTemp{n + pfc - recpfc, 1};          % indices of the paired "source" bubbles in frame n + 2 (pair n + 1), not necessarily in ascending order because it's aligned with the "target" indices in frame n + 1
        % find the common values in the start and end vectors, and the corresponding indices for each
        [recTrackContinuesIndices, ris, rie] = intersect(recStartIndex, recEndIndex, 'stable');
    
        % Trim/update the previous pairs
        bubblePairsPersTemp{n + pfc - recpfc - 1, 1} = bubblePairsPersTemp{n + pfc - recpfc - 1, 1}(ris);
        bubblePairsPersTemp{n + pfc - recpfc - 1, 2} = bubblePairsPersTemp{n + pfc - recpfc - 1, 2}(ris);
    end
    tracks{n} = bubblePairsPersTemp(n : n + pfc, :); % Update the track
end
disp('Tracks with persistence created')
toc
clear bubblePairsPersTemp n pfc recpfc ris rie recTrackContinuesIndices recStartIndex recEndIndex trackContinuesIndices is ie startIndex endIndex

%% 7. Clean tracks

% Turn tracks into a proper link of coordinates and indices - remove the
%   redundant cross-frame stuff.
% Each element in tracksClean is a matrix with each row corresponding to a
%   bubble within the track. Each row represents: 
%   [bubble index within the respective frame, x coordinate, y coordinate, z coordinate, frame number]

tracksClean = cell(size(tracks));   % Initialize the cleaned tracks cell array
nbitAll = zeros(size(tracks));      % # bubbles in each track
for n = startFrame:size(tracks, 1)  % Go through each track
    tracksTemp = tracks{n};         % Get the track at index n
    if ~isempty(tracksTemp)         % Only process if there are bubbles in the current track
        stt = size(tracksTemp);     % Size of temp track [# persistence frames, 2]
        nbif = length(tracksTemp{1}); % # of paired bubbles in each frame of the track (is the same for each entry in the track)
        nbitAll(n) = nbif;          % Store the # bubbles in each frame of track n
        tracksClean{n} = zeros((stt(1) + 1) * nbif, 5); % There are stt(1) + 1 frames represented in each entry of tracksTemp
        for fn = 1:stt(1)           % Go through each frame in the track and store the data
            tracksClean{n}((fn) * nbif + 1 : (fn + 1) * nbif, 1) = tracksTemp{fn, 2}; % Store bubble indices
            tracksClean{n}((fn) * nbif + 1 : (fn + 1) * nbif, 2:4) = centerCoords_corrected{n + fn}(tracksTemp{fn, 2}, :); % Store bubble coordinates
            tracksClean{n}((fn) * nbif + 1 : (fn + 1) * nbif, 5) = repmat(n + fn, nbif, 1); % Store frame numbers
        end
        % Store the info for the first frame in the track
        tracksClean{n}((0) * nbif + 1 : (1) * nbif, 1) = tracksTemp{1, 1};
        tracksClean{n}((0) * nbif + 1 : (1) * nbif, 2:4) = centerCoords_corrected{n}(tracksTemp{1, 1}, :);
        tracksClean{n}((0) * nbif + 1 : (1) * nbif, 5) = repmat(n, nbif, 1);
    else
        tracksClean{n} = NaN;
    end
end
disp('Tracks cleaned')
clear n tracksTemp nbif stt

%% 8: Combining tracks
tracksCleanDynamic = tracksClean; % Make a copy of the tracksClean that we can delete from as we combine tracks
tracksIndividualCombined = {}; % Matrix storage of the bubble data. Each element in the cell is a [# bubbles per frame in the track, 9, # persistence frames] matrix. Each row corresponds to [bubble f x coord, bubble f y coord, bubble f z coord, bubble f+1 x coord, bubble f+1 y coord, bubble f+1 z coord, x velocity, y velocity, z velocity]
nbitAllDynamic = nbitAll;
tic
for tci = startFrame:length(tracksClean)     % track index - go through each collection of tracks
% for tci = 43232:length(tracksClean)
    disp("tci = " + num2str(tci))
    tracksTemp = tracksCleanDynamic{tci};           % Get track collection at tci
    nbitci = nbitAllDynamic(tci);                    % # of bubbles in the track collection starting in index tci

    if nbitci > 0

        for bi = 1:nbitci % bubble index: go through each individual bubble in that track collection
%         for bi = 2
            nbitci = nbitAllDynamic(tci); 
            track_bi_indices = (0:pers) .* nbitci + bi; % Get the indices within tracksTemp for every frame within the p+1-frame track for a single tracked bubble
%             track_bi = tracksTemp(track_bi_indices, :);
            track_bi_updated = tracksTemp(track_bi_indices, :);
        
%             % Start storing each individual track_bi, and we can add to it
%             % if it keeps going in subsequent frame collections
%             tracksIndividualCombined{end + 1} = track_bi;

            % Attempt to find overlapping tracks with track_bi in the next "frame"s of tracksClean 
%             track_bi_updated = track_bi;
            ind = 0;
            flag = true;
            while flag & ((tci + ind + 1) <= length(tracksClean))
                ind = ind + 1;
                tracksNextTemp = tracksCleanDynamic{tci + ind}; % Get each subsequent collection of frames, ind away from tci
                nbitcipi = nbitAllDynamic(tci + ind);                    % # of bubbles in the track collection starting in index tci + ind
%                 checkFrame_inInd = tracksNextTemp( (ind - 1) * (nbitcipi) + 1 : ind * nbitcipi, :); % Get all the bubbles in frame tci + ind to check for overlaps
                checkFrame_inInd = tracksNextTemp( 1 : nbitcipi, :); % Get all the bubbles in frame tci + ind to check for overlaps

                bubbleInfoAt_tci_plus_ind = track_bi_updated(ind + 1, :); % Get the bubble info in the original track_bi at frame tci + ind (template to try to overlay onto the next pers-"frameshift")
                
                [~, ~, ie] = intersect(bubbleInfoAt_tci_plus_ind, checkFrame_inInd, 'rows', 'stable'); % See if there are any intersecting rows (if the track continues into that next "frameshift"
                if ~isempty(ie) % If there is an intersection
                    % Get the new bubble position
                    newBubble_index_frame_tci_plus_ind_plus_pers = pers .* nbitcipi + ie; % Get the indices within tracksTemp for every frame within the p+1-frame track for a single tracked bubble

                    newBubbleInfo = tracksNextTemp(newBubble_index_frame_tci_plus_ind_plus_pers, :);
                    track_bi_updated(end + 1, :) = newBubbleInfo; % Store the new bubble info into the continuing track_bi (NOT SURE IF I NEED TO MAKE A NEW VARIABLE)
                
                    % Delete tracks from the subsequent collection of
                    % frames so we don't get repeated individual combined
                    % tracks
                    tempRemoveIndices = (0:pers) .* nbitcipi + ie; % Remove all the pers # of bubbles within that next collection
                    tracksCleanDynamic{tci + ind}(tempRemoveIndices, :) = [];
%                     nbitAllDynamic(tci + ind) = size(tracksCleanDynamic{tci + ind}, 1);
                    nbitAllDynamic(tci + ind) = nbitAllDynamic(tci + ind) - 1; % Adjust the # of bubbles array so the indexing stays correct

                else
                    flag = false; % No intersection --> stop trying to find further overlaps
                    % ---- Could add a gap filling parameter here ---- %

%                     % ==== CAN ALSO DELETE THE REST OF A TRACK IN THE
%                     % SUBSEQUENT "FRAME COLLECTIONS" IF THERE IS NO INITIAL
%                     % INTERSECTION (FOR SPEED) ==== %
%                     if ind == 1
%                         
%                     end
                end
                
            end

            % Store the updated track_bi as its own individual track
            tracksIndividualCombined{end + 1} = track_bi_updated;

        end

%         for fn = 1:pers                     % Go through all the frames in the tracks with origin frame ti
%             startPoints = tracksTemp((fn - 1) * nbiti + 1 : fn * nbiti, 2:4);
%             endPoints = tracksTemp((fn) * nbiti + 1 : (fn + 1) * nbiti, 2:4);
%             vfn = (endPoints - startPoints) ./ timePerFrame; % velocity = displacement/time
% %             bVelocityC{ti, fn} = [startPoints, endPoints, vfn];  % each row is [x start coord, y start coord, z start coord, x end coord, y end coord, z end coord, x velocity, y velocity, z velocity]
%             bVelocityM{ti}(:, :, fn) = [startPoints, endPoints, vfn];
%         end
    end
end
toc
tracksIndividualCombined = tracksIndividualCombined'; % make it a vertical cell array

%% 8.5. Test: plot all of the individual combined tracks
figure; hold on
% for test_ind = 1:length(tracksIndividualCombined)
for test_ind = 35000:36000
%     plot3(tracksIndividualCombined{ test_ind }(:, 2), tracksIndividualCombined{ test_ind }(:, 3), tracksIndividualCombined{ test_ind }(:, 4), '-o')
    plot3(tracksIndividualCombined{ test_ind }(:, 2), tracksIndividualCombined{ test_ind }(:, 3), tracksIndividualCombined{ test_ind }(:, 4), '.')
end
hold off
clearvars test_ind
%% 9. Create the velocity map

% tracksV: adds the x, y, z velocity within a track, in [voxels/second]
tracksV = cell(size(tracksIndividualCombined)); % Matrix storage of the bubble data. Each element in the cell is a [# bubbles per frame in the track, 9, # persistence frames] matrix. Each row corresponds to [bubble f x coord, bubble f y coord, bubble f z coord, bubble f+1 x coord, bubble f+1 y coord, bubble f+1 z coord, x velocity, y velocity, z velocity]
tic

% Note: could use cellfun to make this more efficient
for ti = 1:length(tracksV)     % track index - go through each track
% for ti = 1:10000
    trackTemp = tracksIndividualCombined{ti};           % Get track ti
    nfit = length(trackTemp);                    % # of frames in trackTemp

    % Reorganize the data matrix for each track. 
    % The columns are: [frame #, paired bubble index in that frame (end side), x, y, z]
    trackReorganized = [trackTemp(:, 5), trackTemp(:, 1:4)]; % Move the frame column to the left

    trackReorganized = [trackReorganized(1:end-1, :), trackReorganized(2:nfit, 3:5)]; % Add "end points" next to the "start points" for each velocity calculation.
    % Note: have to remove the last indices for the last point since the end point is included in the row above

    trackReorganized(:, 9:11) = ( trackReorganized(:, 6:8) - trackReorganized(:, 3:5) ) ./ timePerFrame; % velocity = displacement/time [voxels per second]

    tracksV{ti} = trackReorganized; % Store the reorganized track info
end
toc
disp('Velocity map created')
clear ti fn trackTemp startPoints endPoints vfn nbif

%% 9.5.1. 3D plots of the tracks
track3DPlot(tracksV, 1:length(tracksV))
% track3DPlot(tracksV, 35000:36000)

%% 9.5.2. Spaghetti plot of speed along each track
trackSpeedSpaghettiPlot(tracksV, 35000:35010)

%% 10. Smooth velocities across each track
% Optionally, remove combined tracks that violate a distance criterion after applying a moving mean (low pass)

tracksVS = tracksV; % Distance criterion applied (?)
tic
tracksToKeep_DC = true(size(tracksVS)); % Mask for which tracks to keep, for the distance criterion

% Note: could use cellfun to make this more efficient
for ti = 1:length(tracksVS)     % track index - go through each track
% for ti = 30000
    trackTemp = tracksV{ti};           % Get track ti
    vTemp = trackTemp(:, 9:11);
    vTempMM = movmean(vTemp, mmws, 5); % Moving mean in time on the velocities
%     figure; plot(sqrt(sum(vTemp .^ 2, 2)))
%     figure; plot(sqrt(sum(vTempMM .^ 2, 2)))

    % Distance criteria
%     if any(abs(vTempMM) > maxPixelDistPerFrame / timePerFrame, 'all')
%         tracksToKeep_DC(ti) = false;
%     end

    tracksVS{ti}(:, 9:11) = vTempMM; % Store the reorganized track info
end

tracksVS(~tracksToKeep_DC) = [];
toc
disp('Velocity map smoothed')
% disp('Distance (velocity) criteria applied to velocity map')
clear ti fn trackTemp startPoints endPoints vfn nbif

%% 10.5.2 Scale the smoothed velocity into [mm/s]
tracksVS_MMS = tracksVS;
for ti = 1:length(tracksVS_MMS)     % track index - go through each track
    trackTemp = tracksVS_MMS{ti};
    tracksVS_MMS{ti}(:, 9:11) = tracksVS_MMS{ti}(:, 9:11) ./ [xpixelsPerM, ypixelsPerM, zpixelsPerM] .* 1e3;
    % Note: I'm being sloppy with which dimension is x and y because they
    % will usually be the same in voxel size or # of voxels
end

disp('Velocity map converted to mm/s')
%% 10.5.1. 3D plots of the tracks after distance criterion applied
% track3DPlot(tracksVS, 1:length(tracksVS))

track3DPlot(tracksVS, 42000:42100)

% My finding on 10/20/25 of using the distance threshold on the raw
% smoothed velocities: removes way too many (probably too noisy)

%% 10.5.2. Spaghetti plot of speed along each track after smoothing
% trackSpeedSpaghettiPlot(tracksVS_MMS, 140000:140010)
trackSpeedSpaghettiPlot(tracksVS_MMS, 10000 + [0:10])

%% Make a new variable that keeps only tracks with some velocity characteristics
low_prctile = 25;
low_speed_threshold = 2; % mm/s

high_speed_threshold = 50;

tracksVS_MMS_thresholded_mask = false(size(tracksVS_MMS));
tic
for ti = 1:length(tracksVS_MMS)     % track index - go through each track
    trackTemp = tracksVS_MMS{ti};
    speedsTemp = sqrt(sum(trackTemp(:, 9:11) .^ 2, 2));
    
    speed_low_prctile = prctile(speedsTemp, low_prctile);
    if speed_low_prctile < low_speed_threshold & max(speedsTemp) <= high_speed_threshold
        tracksVS_MMS_thresholded_mask(ti) = true;
    end

end
toc
tracksVS_MMS_thresholded = tracksVS_MMS(tracksVS_MMS_thresholded_mask);
clearvars ti trackTemp speedsTemp speed_low_prctile

%% Spaghetti plot for the low velocity tracks
% trackSpeedSpaghettiPlot(tracksVS_MMS_thresholded, 1:length(tracksVS_MMS_thresholded))
trackSpeedSpaghettiPlot(tracksVS_MMS_thresholded, 100 + [1:10])

%% 11. Kalman filter as a function
vMMStoPixelDispPerFrame = timePerFrame / 1e3 * [xpixelsPerM, ypixelsPerM, zpixelsPerM];
tic
tracksVS_KF_MMS = applyKF(tracksVS_MMS, vMMStoPixelDispPerFrame, img_size, pers);
toc

% Cellfun attempt
% vMMStoPixelDispPerFrame = timePerFrame / 1e3 * [xpixelsPerM, ypixelsPerM, zpixelsPerM];
% 
% tic
% % tracksVS_KF_MMS = cellfun(@applyKF, tracksVS_MMS, vMMStoPixelDispPerFrame, img_size);
% applyKFwrapper = @(tracksVS_MMS) applyKF_cf(tracksVS_MMS, vMMStoPixelDispPerFrame, img_size, pers);
% tracksVS_KF_MMS = cellfun(applyKFwrapper, tracksVS_MMS(1:80), 'UniformOutput', false);
% toc
%% 11a. Kalman filter with velocity in the state
% (==== NEED TO FIX THIS, PLUS VELOCITY UNITS ====)

% tracksVS_KF_MMS = cell(size(tracksVS_MMS));
tracksVS_KF_MMS = tracksVS_MMS;

numDims = 3;
numStates = numDims * 2; % In 3D: x position, y position, z position, x displacement, y displacement, z displacement

if pers < 3
    error('The Kalman filter needs at least 3 points in the track')
end

vMMStoPixelDispPerFrame = timePerFrame / 1e3 * [xpixelsPerM, ypixelsPerM, zpixelsPerM];

% Define the matrices that map the state transition, and the state ->
% observation transformation
Fk = [1, 0, 0, 1, 0, 0; ...
      0, 1, 0, 0, 1, 0; ...
      0, 0, 1, 0, 0, 1; ...
      0, 0, 0, 1, 0, 0; ...
      0, 0, 0, 0, 1, 0; ...
      0, 0, 0, 0, 0, 1];
Hk = [1, 0, 0, 0, 0, 0; ...
      0, 1, 0, 0, 0, 0; ...
      0, 0, 1, 0, 0, 0; ...
      0, 0, 0, 1, 0, 0; ...
      0, 0, 0, 0, 1, 0; ...
      0, 0, 0, 0, 0, 1];

%%%%%%%% these covariance matrices are from the Song et al. 2020 paper %%%%%%%%%
Qk = diag(ones(numStates, 1)) .* 0.5;      % Covariance matrix of the system/process noise
Rk = diag(ones(numStates, 1)) .* 4;   % Covariance matrix of the observation noise

tic
for ti = 1:length(tracksVS_MMS) % Go through each track
% for ti = 1:10000
    trackTemp = tracksVS_MMS{ti};
    trackTemp = [trackTemp; [0, 0, trackTemp(end, 6:8), 0, 0, 0, 0, 0, 0]]; % Append the last position onto the column the KF looks at, and fill the rest with zeros (not used except velocities)

    % Initialize variables for the state vector and covariance matrix
    xk = NaN(numStates, size(trackTemp, 1)); % xk has dimensions (6 from xyz pos and velocity, # frames in the track)
    Pk = NaN(numStates, numStates, size(trackTemp, 1));

    % Initial values               
    xk(:, 1) = [trackTemp(1, 3:5), trackTemp(1, 9:11) .* vMMStoPixelDispPerFrame]'; % Initial state vector

    Pk(:, :, 1) = [1, 0, 0, 0, 0, 0; ... % Initial covariance matrix
                   0, 1, 0, 0, 0, 0; ...
                   0, 0, 1, 0, 0, 0; ...
                   0, 0, 0, 10, 0, 0; ...
                   0, 0, 0, 0, 10, 0; ...
                   0, 0, 0, 0, 0, 10];

    % Go through each step of the track
    for k = 2:size(trackTemp, 1)
        % Prediction
        xkp = Fk * xk(:, k - 1); 
        Pkp = Fk * Pk(:, :, k - 1) * Fk + Qk;

        % Observation
        yk = [trackTemp(k, 3:5), trackTemp(k, 9:11) .* vMMStoPixelDispPerFrame]';

        % Update
        Kku = Pkp * Hk' * inv(Hk * Pkp * Hk' + Rk); % Kalman gain matrix
        Iku = yk - Hk * xkp; % Innovation vector (difference between the observed state and the predicted state transformed into an observation at step k)
        xku = xkp + Kku * Iku; % Updated (weighted) estimate for the state at step k
        Pku = (eye(length(xku)) - Kku * Hk) * Pkp; % Updated covariance matrix at step k

        % Store the updated state and covariance
        xk(:, k) = xku;
        Pk(:, :, k) = Pku;

        % Adjustments in case the Kalman filter puts some points outside the original region
        xktemp = xk(1:3, :);
        xktemp(xktemp < 1) = 1;
        xk(1:3, :) = xktemp;
        xk(1, xk(1, :) > img_size(1)) = img_size(1);
        xk(2, xk(2, :) > img_size(2)) = img_size(2);
        xk(3, xk(3, :) > img_size(3)) = img_size(3);

    end
    % Store the post-KF values
    xkt = xk'; % Temporarily store the transpose of xk
    tracksVS_KF_MMS{ti}(:, 3:5) = round(xkt(1:end-1, 1:3)); % "Start" positions
    tracksVS_KF_MMS{ti}(:, 6:8) = round(xkt(2:end, 1:3)); % "End" positions
%     tracksVS_KF_MMS{ti}(:, 9:11) = xkt(2:end, 4:6) ./ vMMStoPixelDispPerFrame;
    tracksVS_KF_MMS{ti}(:, 9:11) = xkt(1:end - 1, 4:6) ./ vMMStoPixelDispPerFrame; % Velocities
end

toc
disp('Kalman filter applied')
clear n tn tln k xk Pk yk Kku Iku xku Pku track

%% 11.5. Plot some stuff post-KF
track3DPlot(tracksVS_KF_MMS, 42000:42100)

trackSpeedSpaghettiPlot(tracksVS_KF_MMS, 30000:40010)

%% Filter low speed tracks post-KF
low_prctile = 25;
low_speed_threshold = 2; % mm/s
high_speed_threshold = 50; % mm/s

[tracksVS_KF_MMS_filtered] = refineLowSpeedTracks(tracksVS_KF_MMS, low_prctile, low_speed_threshold, high_speed_threshold);

%% Look at low speed tracks post-KF
% track3DPlot(tracksVS_KF_MMS_filtered, 100 + [1:10])
track3DPlot(tracksVS_KF_MMS_filtered, 1:length(tracksVS_KF_MMS_filtered))
% trackSpeedSpaghettiPlot(tracksVS_KF_MMS_filtered, 200 + [1:10])

%% Interpolated speed map for the non-constrained KF, low speed filtered data
[SM_SmoothedKF_LS_LI, SM_SmoothedKF_LS_LI_counter] = interpolatedSpeedMap(tracksVS_KF_MMS_filtered, img_size, maxPixelDistPerFrame); % flow speed map, linearly interpolated
SM_SmoothedKF_LS_LI_Rfn = SM_SmoothedKF_LS_LI;
SM_SmoothedKF_LS_LI_Rfn = thresholdMaps(SM_SmoothedKF_LS_LI_Rfn, SM_SmoothedKF_LS_LI_counter, 1, 300);

cmap = colormap_ULM;
figure; imagesc(squeeze(max(SM_SmoothedKF_LS_LI_Rfn(:, :, :), [], 1))'); colormap(cmap); clim([0, maxSpeedExpectedMMPerS])
figure; imagesc(squeeze(max(SM_SmoothedKF_LS_LI_Rfn(:, :, :), [], 3))'); colormap(cmap); clim([0, maxSpeedExpectedMMPerS]); axis square

%% 12. Acceleration and direction constraints

tic
% bVelocityConstrained = applyConstraints(bVelocityM, vTrimmedMeanPercentage, aThresholdFactor, angleChangeThreshold, timePerFrame);

% bVelocityMSmoothedMMSConstrained = applyConstraints(bVelocityMSmoothedMMS, vTrimmedMeanPercentage, aThresholdFactor, angleChangeThreshold, timePerFrame);

bVelocityMSmoothedKFConstrainedMMS = applyConstraints(bVelocityMSmoothedKFMMS, vTrimmedMeanPercentage, aThresholdFactor, angleChangeThreshold, timePerFrame);
toc
disp('Acceleration and direction constraints applied')
clear n tln tn trackAlreadyDeleted track vTrack vTrackTrimmedMean aThresholdMag aTrackMag angleTrack angelTrackChanges

%%
% addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Processing\Jerman Enhancement Filter\')
% test = vesselness3D(bSum .^ 0.5, 1:5, [xpix_spacing; ypix_spacing; zpix_spacing], 0.5, true);
% volumeViewer(test)

%% Plot the interpolated density map with video
actualSize = [lateral_width, lateral_width, axial_depth];
BDM_video = interpolatedDensityMapWithVideo(bVelocityMSmoothedMMSConstrained, img_size, startFrame, maxPixelDistPerFrame, actualSize); % test density map interpolated
% BDM_video = interpolatedDensityMapWithVideo(bVelocityM, img_size, startFrame, maxPixelDistPerFrame, actualSize); % test density map interpolated

% %% No KF test
% [SM_Smoothed_LI, SM_Smoothed_LI_counter] = interpolatedSpeedMap(tracksVS_MMS, img_size, maxPixelDistPerFrame); % flow speed map, linearly interpolated
% SM_Smoothed_LI_Rfn = SM_Smoothed_LI;
% SM_Smoothed_LI_Rfn = thresholdMaps(SM_Smoothed_LI_Rfn, SM_Smoothed_LI_counter, 2, 300);


%% Get the speed maps
% Non-constrained no KF
% [SM_LI, SM_LI_counter] = interpolatedSpeedMap(tracksVS_MMS, img_size, maxPixelDistPerFrame); % flow speed map, linearly interpolated
% SM_LI_Rfn = SM_LI;
% SM_LI_Rfn = thresholdMaps(SM_LI_Rfn, SM_LI_counter, 1, 300);

% Non-constrained KF
[SM_SmoothedKF_LI, SM_SmoothedKF_LI_counter] = interpolatedSpeedMap(tracksVS_KF_MMS, img_size, maxPixelDistPerFrame); % flow speed map, linearly interpolated
SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI;
SM_SmoothedKF_LI_Rfn = thresholdMaps(SM_SmoothedKF_LI_Rfn, SM_SmoothedKF_LI_counter, 1, 300);

% Constrained KF
% [SM_SmoothedKFConstrained_LI, SM_SmoothedKFConstrained_LI_counter] = interpolatedSpeedMap(bVelocityMSmoothedKFConstrainedMMS, img_size, startFrame, maxPixelDistPerFrame); % flow speed map, linearly interpolated
% SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI;
% SM_SmoothedKFConstrained_LI_Rfn = thresholdMaps(SM_SmoothedKFConstrained_LI_Rfn, SM_SmoothedKFConstrained_LI_counter, 2, 300);

% Look at the smoothed constrained no KF data %%%%%%
% [SM_SC_LI, SM_SC_LI_counter] = interpolatedSpeedMap(bVelocityMSmoothedMMSConstrained, img_size, startFrame, maxPixelDistPerFrame); % flow speed map, linearly interpolated
% SM_SC_LI_Rfn = SM_SC_LI;
% SM_SC_LI_Rfn = thresholdMaps(SM_SC_LI_Rfn, SM_SC_LI_counter, 2, 300);

%% Generate a file for other software to read
% writematrix(SM_SmoothedKFConstrained_LI_Rfn, 'D:\Allen\Data\AZ02 Stroke ULM RC15gV\04-22-2025 7d left eye\ULM subpixel processing results\Speed maps\vol.dat')

%% Plot speed map

% volumeViewer(SM_LI_Rfn)
% plotMIPs(SM_LI_RSC, 1)

% volumeViewer(SM_SmoothedKFConstrained_LI_Rfn)

cmap = colormap_ULM;
% figure; imagesc(squeeze(max(SM_LI_Rfn(:, :, :), [], 1))'); colormap(cmap); clim([0, maxSpeedExpectedMMPerS])
% figure; imagesc(squeeze(max(SM_LI_Rfn(:, :, :), [], 3))'); colormap(cmap); clim([0, maxSpeedExpectedMMPerS]); axis square
% figure; imagesc(squeeze(max(SM_LI(:, :, :), [], 1))'); colormap(cmap); clim([0, maxSpeedExpectedMMPerS])
% figure; imagesc(squeeze(max(SM_LI(:, :, :), [], 3))'); colormap(cmap); clim([0, maxSpeedExpectedMMPerS]); axis square

figure; imagesc(squeeze(max(SM_SmoothedKF_LI_Rfn(:, :, :), [], 1))'); colormap(cmap); clim([0, maxSpeedExpectedMMPerS])
figure; imagesc(squeeze(max(SM_SmoothedKF_LI_Rfn(:, :, :), [], 3))'); colormap(cmap); clim([0, maxSpeedExpectedMMPerS]); axis square
% % figure; imagesc(squeeze(max(SM_SmoothedKF_LI_counter(300:500, :, :), [], 1) .^ 0.7)'); colormap hot

% figure; imagesc(squeeze(max(SM_LI_Rfn(300:500, :, :), [], 1))'); colormap(cmap);

%%
test = SM_SmoothedKFConstrained_LI_Rfn;
testlim = 10;
test(test >= testlim) = 0;
figure; imagesc(squeeze(max(test(400:600, :, :), [], 1))'); colormap(cmap); clim([0, testlim])

%% Make speed map MIPs
% [cmap, ~, ~, ~, ~] = Colormaps_fUS;
% [~, ~, cmap, ~, ~] = Colormaps_fUS;
% cmap = 'turbo';
cmap = colormap_ULM;
% plotSpeedMIPs(SM_LI_RSC, 1)
% generateTiffStack_multi([{SM_LI_RSC}], [8.8, 8.8, 8], cmap, 10)
% generateTiffStack_multi([{SM_SC_LI_Rfn}], [8.8, 8.8, 8], cmap, 10, [0, 50])
% generateTiffStack_multi([{SM_SmoothedKFConstrained_LI_Rfn}], [8.8, 8.8, 8], cmap, 50, [0, 40])
generateTiffStack_multi([{SM_SmoothedKF_LI_Rfn}], [8.8, 8.8, 8], cmap, 50, [0, 40])
% generateTiffStack_multi([{SM_SmoothedKF_LI_Rfn}], [8.8, 8.8, 8], cmap, 1, [0, 40])
% generateTiffStack_multi([{SMs_AZ04_baseline.SM_SmoothedKF_LI_Rfn}], [8.8, 8.8, 8], cmap, 1, [0, 40])
% generateTiffStack_multi([{test}], [8.8, 8.8, 8], cmap, 50, [0, testlim])

% generateTiffStack_multi([{SM_SmoothedKFConstrained_LI_Rfn}], [8.8, 8.8, 8], cmap, 1)

%% Convert the speed maps to a structure
% SMs_AZ02_baseline.SM_LI = SM_LI;
% SMs_AZ02_baseline.SM_LI_counter = SM_LI_counter;
% SMs_AZ02_baseline.SM_LI_RSC = SM_LI_RSC;
% SMs_AZ02_baseline.SM_SC_LI_counter = SM_SC_LI_counter;
% SMs_AZ02_baseline.SM_SC_LI_RSC = SM_SC_LI_RSC;
% SMs_AZ02_baseline.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ02_baseline.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
% SMs_AZ02_baseline.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;

% % SMs_AZ02_hour1.SM = SM;
% SMs_AZ02_hour1.SM_LI = SM_LI;
% SMs_AZ02_hour1.SM_LI_counter = SM_LI_counter;
% SMs_AZ02_hour1.SM_LI_RSC = SM_LI_RSC;
% SMs_AZ02_hour1.SM_SC_LI_counter = SM_SC_LI_counter;
% SMs_AZ02_hour1.SM_SC_LI_RSC = SM_SC_LI_RSC;
% SMs_AZ02_hour1.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ02_hour1.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
% SMs_AZ02_hour1.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;

% SMs_AZ02_day3.SM_LI = SM_LI;
% SMs_AZ02_day3.SM_LI_counter = SM_LI_counter;
% SMs_AZ02_day3.SM_LI_RSC = SM_LI_RSC;
% SMs_AZ02_day3.SM_SC_LI_counter = SM_SC_LI_counter;
% SMs_AZ02_day3.SM_SC_LI_RSC = SM_SC_LI_RSC;
% SMs_AZ02_day3.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ02_day3.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
% SMs_AZ02_day3.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;

% SMs_AZ02_day7.SM_LI = SM_LI;
% SMs_AZ02_day7.SM_LI_counter = SM_LI_counter;
% SMs_AZ02_day7.SM_LI_RSC = SM_LI_RSC;
% SMs_AZ02_day7.SM_SC_LI_counter = SM_SC_LI_counter;
% SMs_AZ02_day7.SM_SC_LI_RSC = SM_SC_LI_RSC;
% SMs_AZ02_day7.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ02_day7.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
% SMs_AZ02_day7.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;
SMs_AZ02_day7.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
SMs_AZ02_day7.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
SMs_AZ02_day7.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

% SMs_AZ03_baseline.SM_LI = SM_LI;
% SMs_AZ03_baseline.SM_LI_counter = SM_LI_counter;
% SMs_AZ03_baseline.SM_LI_Rfn = SM_LI_Rfn;
% SMs_AZ03_baseline.SM_SC_LI_counter = SM_SC_LI_counter;
% SMs_AZ03_baseline.SM_SC_LI_Rfn = SM_SC_LI_Rfn;
% SMs_AZ03_baseline.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ03_baseline.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
% SMs_AZ03_baseline.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;
SMs_AZ03_baseline.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
SMs_AZ03_baseline.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
SMs_AZ03_baseline.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

% SMs_AZ03_hour1.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ03_hour1.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
% SMs_AZ03_hour1.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;

% SMs_AZ03_day1.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ03_day1.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
% SMs_AZ03_day1.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;
% SMs_AZ03_day1.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
% SMs_AZ03_day1.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
% SMs_AZ03_day1.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

% SMs_AZ03_day4.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ03_day4.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
% SMs_AZ03_day4.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;
% SMs_AZ03_day4.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
% SMs_AZ03_day4.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
% SMs_AZ03_day4.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

% SMs_AZ03_day7.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ03_day7.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
% SMs_AZ03_day7.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;
% SMs_AZ03_day7.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
% SMs_AZ03_day7.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
% SMs_AZ03_day7.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

% SMs_AZ04_baseline.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
% SMs_AZ04_baseline.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
% SMs_AZ04_baseline.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;
SMs_AZ04_hour1.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
SMs_AZ04_hour1.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
SMs_AZ04_hour1.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

% SMs_AZ04_day1.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
% SMs_AZ04_day1.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
% SMs_AZ04_day1.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

% SMs_AZ04_day3.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
% SMs_AZ04_day3.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
% SMs_AZ04_day3.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

% SMs_AZ04_day7.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
% SMs_AZ04_day7.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
% SMs_AZ04_day7.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

% SMs_AZ04_day14.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
% SMs_AZ04_day14.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
% SMs_AZ04_day14.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

% SMs_AZ06_day1.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
% SMs_AZ06_day1.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
% SMs_AZ06_day1.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

% SMs_AZ06_day4.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
% SMs_AZ06_day4.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
% SMs_AZ06_day4.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

SMs_AZ06_day7.SM_SmoothedKF_LI = SM_SmoothedKF_LI;
SMs_AZ06_day7.SM_SmoothedKF_LI_Rfn = SM_SmoothedKF_LI_Rfn;
SMs_AZ06_day7.SM_SmoothedKF_counter = SM_SmoothedKF_LI_counter;

%% Helper functions

% function MIPvideo(bSum, xws, yws, zws, framerate, power) % Define x, y, z window sizes for a MIP flythrough video of the bubble density map
%     img_size = size(bSum);
%     savepath = [uigetdir('', 'Select the save path'), '\'];
% 
%     filename = datestr(now, 0);
% 
%     vo = VideoWriter([savepath, filename]);
%     vo.Quality = 100;
%     vo.FrameRate = framerate;
%     open(vo);
% 
%     vf = figure;
%     
%     for f = 1:img_size(1) - xws
% %     for f = 1:100
%         v1f = max(bSum(f:xws, :, :), [], 1); % MIP across x
%         
%         imagesc(squeeze(v1f)')
%         
%         cv = getframe(vf);     % get the current volume
%         rgb = frame2im(cv);      % convert the frame to rgb data
%     
%         writeVideo(vo, rgb);
%     end
%     close(vo);
% 
% end

function bVelocityConstrained = applyConstraints(bVelocity, vTrimmedMeanPercentage, aThresholdFactor, angleChangeThreshold, timePerFrame)
    bVelocityConstrained = bVelocity; % Initialize the variable for the smoothed, Kalman filtered, constrained velocity map
    parfor n = 1:size(bVelocityConstrained, 1) % Changed to parfor 5/8/25
        tln = bVelocityConstrained{n};    % Track list n
        for tn = size(tln, 1):-1:1                      % Go through each track number tn
            trackAlreadyDeleted = false;                % Reset the flag
            track = squeeze(tln(tn, :, :))';            % Get the track
            vTrack = track(:, 7:9);                     % Velocities of the track
    
            % Acceleration constraint
            vTrackTrimmedMean = trimmean(vTrack, vTrimmedMeanPercentage); % Exclude some percentage of the values when taking the mean. It goes across the first non-singleton dimension.
            aThresholdMag = abs(aThresholdFactor .* vTrackTrimmedMean ./ timePerFrame);
            aTrackMag = abs(diff(vTrack, 1) ./ timePerFrame); % Accelerations of the track
            if any(aTrackMag > aThresholdMag, 'all') % Remove the track if the acceleration constraint is violated
                bVelocityConstrained{n}(tn, :, :) = [];
                trackAlreadyDeleted = true;
            end
    
            % Direction constraint
            if ~trackAlreadyDeleted % Don't need to go through the direction calculation if the track was already deleted for the acceleration constraint
                angleTrack = atan2(vTrack(:, 2), vTrack(:, 1));         % Angle of each segment on the track
                angleTrackChanges = diff(angleTrack);                   % Change in angle between segments on the track
                if any(abs(angleTrackChanges) > angleChangeThreshold)   % Apply the threshold
                    bVelocityConstrained{n}(tn, :, :) = [];
                end
            end
        end
    end
end

function [densityMapInterpolated] = interpolatedDensityMap(bVelocityM, img_size, startFrame, maxPixelDistPerFrame)
    densityMapInterpolated = zeros(img_size(1), img_size(2), img_size(3));
%     plotPower = 1;

    % Counters for proper averaging if there are overlapped pixels from
    % different tracks
    densityMapInterpolatedCounter = zeros(size(densityMapInterpolated));
    
    tic
    for ti = startFrame:size(bVelocityM, 1)
%     for ti = startFrame:startFrame+100
%     for ti = 15000:15100
%     for ti = 12000:size(bVelocityM, 1)
        bvTemp = bVelocityM{ti}; % get the ti-th entry
        pers = size(bvTemp, 3);
        if ~isempty(bvTemp) % only do stuff if the bubble velocity cell array entry is not empty
            for bpi = 1:size(bvTemp, 1) % bubble pair index
    
                % Initialize temporary start and end coordinate matrices.
                % Each have dimensions [# persistence frames, 3] where each row is [x coord, y coord, z coord].
                coordsStart = NaN(pers, 3);
                coordsEnd = NaN(pers, 3);
    
                % Go through the # of persistence frames and get the
                % coordinates at each frame pfi for the bubble track bpi.
                %   Each row corresponds to a persistence frame, and contains
                %   [x, y, z] velocity.
                for pfi = 1:pers % persistence frame index
                    coordsStart(pfi, :) = bvTemp(bpi, 1:3, pfi);
                    coordsEnd(pfi, :) = bvTemp(bpi, 4:6, pfi);
                end

                vecDist = coordsEnd - coordsStart;
%                 totalDist = sqrt(sum(vecDist .^ 2, 2));

                % Only interpolate if the distance between points in a
                % track is less than the max pixel dist per frame as
                % calculated before
                if all(abs(vecDist) - maxPixelDistPerFrame <= 0, 'all')
                    vTemp = squeeze(bvTemp(bpi, 7:9, :)); % Velocity components
%                 speedTemp = sqrt(vTemp(:, 1).^2 + vTemp(:, 2).^2 + vTemp(:, 3).^2); % Speed vector: one value per persistence frame index
                    speedTemp = sqrt(sum(vTemp.^2, 1))';
                    roundOrNot = true;
                    interpPts = ULM_interp3D_linear(coordsStart, coordsEnd, speedTemp, roundOrNot); % Get interpolated points with the corresponding z velocity value. each row is [z coord, x coord, z velocity]
    
    %                 figure; scatter3(interpPts(:, 1), interpPts(:, 2), interpPts(:, 3))
                    for ipi = 1:size(interpPts, 1) % interpolated point index
                        interpPtsTemp = interpPts(ipi, 1:3);
                        speedValTemp = interpPts(ipi, 4);
                        
    %                     speedMap(smoothedPtsTemp(1), smoothedPtsTemp(2), smoothedPtsTemp(3)) = speedMap(smoothedPtsTemp(1), smoothedPtsTemp(2), smoothedPtsTemp(3)) + speedValTemp;
                        densityMapInterpolatedCounter(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) = densityMapInterpolatedCounter(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) + 1;
                    end
                end
                
%                 if any(totalDist > pixel)
    
                
    
    %             for ipi = 1:size(interpPts, 1) % interpolated point index
    %                 interpPtsTemp = interpPts(ipi, :);
    %                 speedValTemp = interpPtsTemp(4);
    %                 
    %                 speedMap(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) = speedMap(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) + speedValTemp;
    %                 speedMapCounter(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) = speedMapCounter(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) + 1;
    %             end
            end
    %     else
    %         interpPts = [];
    %         coordsStart = [];
    %         coordsEnd = [];
    %         zvTemp = [];
    %         bvTemp = [];
    %         zVelTemp = [];
    %         interpPtsTemp = [];
        end
    %     disp(strcat("Frame ", num2str(ti), " stored."))
    %     clear bvTemp
    end
    disp('Density map interpolated')
    toc
%     clear speedValTemp ti bpi pfi ipi interpPtsTemp speedTemp coordsStart coordsEnd bvTemp
    densityMapInterpolated = densityMapInterpolatedCounter;
    % Take the average for pixels with overlapping tracks
%     speedMask = densityMapInterpolatedCounter > 0;
%     speedMap(speedMask) = speedMap(speedMask) ./ densityMapInterpolatedCounter(speedMask);
end

function [speedMapInterpolated, speedMapInterpolatedCounter] = interpolatedSpeedMap(tracksVS_KF_MMS, img_size, maxPixelDistPerFrame)
    speedMapInterpolated = zeros(img_size(1), img_size(2), img_size(3));

    % Counters for proper averaging if there are overlapped pixels from
    % different tracks
    speedMapInterpolatedCounter = zeros(size(speedMapInterpolated));
    
    tic
    for ti = 1:length(tracksVS_KF_MMS)
%     for ti = startFrame:startFrame+100
%     for ti = 15000:15100
%     for ti = 12000:size(bVelocityM, 1)
        trackTemp = tracksVS_KF_MMS{ti}; % get the ti-th entry
        nfit = size(trackTemp, 1); % # frames in track

        for fi = 1:nfit % frame index (within the track)

            % Initialize temporary start and end coordinate matrices.
            % Each have dimensions [# frames in the track, 3] where each row is [x coord, y coord, z coord].
%             coordsStart = NaN(nfit, 3);
%             coordsEnd = NaN(nfit, 3);

            % Go through the # of persistence frames and get the
            % coordinates at each frame pfi for the bubble track bpi.
            %   Each row corresponds to a persistence frame, and contains
            %   [x, y, z] velocity.
%             for pfi = 1:nfit % persistence frame index
%                 coordsStart(pfi, :) = trackTemp(fi, 1:3, pfi);
%                 coordsEnd(pfi, :) = trackTemp(fi, 4:6, pfi);
%             end
            coordsStart = trackTemp(fi, 3:5);
            coordsEnd = trackTemp(fi, 6:8);

            vecDist = coordsEnd - coordsStart;
%                 totalDist = sqrt(sum(vecDist .^ 2, 2));

            % Only interpolate if the distance between points in a
            % track is less than the max pixel dist per frame as
            % calculated before
            if all(abs(vecDist) - maxPixelDistPerFrame <= 0, 'all')
                vTemp = squeeze(trackTemp(fi, 9:11)); % Velocity components for the track # bpi
%                 speedTemp = sqrt(vTemp(:, 1).^2 + vTemp(:, 2).^2 + vTemp(:, 3).^2); % Speed vector: one value per persistence frame index
                speedTemp = sqrt(sum(vTemp.^2, 2))';
                roundOrNot = true;
                interpPts = ULM_interp3D_linear(coordsStart, coordsEnd, speedTemp, roundOrNot); % Get interpolated points with the corresponding z velocity value. each row is [z coord, x coord, z velocity]

%                 figure; scatter3(interpPts(:, 1), interpPts(:, 2), interpPts(:, 3))
                for ipi = 1:size(interpPts, 1) % interpolated point index
                    interpPtsTemp = interpPts(ipi, 1:3);
                    speedValTemp = interpPts(ipi, 4);
                    
%                     speedMap(smoothedPtsTemp(1), smoothedPtsTemp(2), smoothedPtsTemp(3)) = speedMap(smoothedPtsTemp(1), smoothedPtsTemp(2), smoothedPtsTemp(3)) + speedValTemp;
                    speedMapInterpolated(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) = speedMapInterpolated(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) + speedValTemp; % Accumulate the speed so we can take the average for overlapping tracks at a voxel
                    speedMapInterpolatedCounter(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) = speedMapInterpolatedCounter(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) + 1; % Increment the counter for taking the average
                end
            end
            
        end

    end
    disp('Speed map interpolated')
    toc
%     clear speedValTemp ti bpi pfi ipi interpPtsTemp speedTemp coordsStart coordsEnd bvTemp

    % Take the average speed    for pixels with overlapping tracks
    speedMask = speedMapInterpolatedCounter > 0;
    speedMapInterpolated(speedMask) = speedMapInterpolated(speedMask) ./ speedMapInterpolatedCounter(speedMask);
end

function [densityMapInterpolated] = interpolatedDensityMapWithVideo(bVelocityM, img_size, startFrame, maxPixelDistPerFrame, actualSize)
    densityMapInterpolated = zeros(img_size(1), img_size(2), img_size(3));

    % Set up the figure for the video
    savepath = uigetdir('D:\Allen\Data\', 'Select the save path');
    savepath = [savepath, '\'];

    vo = VideoWriter([savepath, 'volumeVideo']);

    vo.Quality = 100;
    vo.FrameRate = 30;
    open(vo);

    vf = figure;
    V = volshow(densityMapInterpolated);
    V_old = V;
%     V.Alphamap(1:100) = 0;          % Change transparency
    V.BackgroundColor = [1, 1, 1];  % Make background white
    V.ScaleFactors(3) = size(densityMapInterpolated, 1) / size(densityMapInterpolated, 3) * actualSize(1) / actualSize(3); % scale with # pixels and region size
%     V.CameraPosition = V.CameraPosition ./ 2;
%     V.CameraPosition = [2.1161 -3.7332 -0.1764];
%     V.CameraUpVector = [0.1853 -0.3160 -0.9305];
%     V.CameraViewAngle = 15;

%     plotPower = 1;

    % Counters for proper averaging if there are overlapped pixels from
    % different tracks
    densityMapInterpolatedCounter = zeros(size(densityMapInterpolated));
    
    tic
    for ti = startFrame:size(bVelocityM, 1)
%     for ti = startFrame:startFrame+100
%     for ti = 15000:16100
%     for ti = 12000:size(bVelocityM, 1)
        bvTemp = bVelocityM{ti}; % get the ti-th entry
        pers = size(bvTemp, 3);
        if ~isempty(bvTemp) % only do stuff if the bubble velocity cell array entry is not empty
            for bpi = 1:size(bvTemp, 1) % bubble pair index
    
                % Initialize temporary start and end coordinate matrices.
                % Each have dimensions [# persistence frames, 3] where each row is [x coord, y coord, z coord].
                coordsStart = NaN(pers, 3);
                coordsEnd = NaN(pers, 3);
    
                % Go through the # of persistence frames and get the
                % coordinates at each frame pfi for the bubble track bpi.
                %   Each row corresponds to a persistence frame, and contains
                %   [x, y, z] velocity.
                for pfi = 1:pers % persistence frame index
                    coordsStart(pfi, :) = bvTemp(bpi, 1:3, pfi);
                    coordsEnd(pfi, :) = bvTemp(bpi, 4:6, pfi);
                end

                vecDist = coordsEnd - coordsStart;
%                 totalDist = sqrt(sum(vecDist .^ 2, 2));

                % Only interpolate if the distance between points in a
                % track is less than the max pixel dist per frame as
                % calculated before
                if all(abs(vecDist) - maxPixelDistPerFrame <= 0, 'all')
                    vTemp = squeeze(bvTemp(bpi, 7:9, :)); % Velocity components
%                 speedTemp = sqrt(vTemp(:, 1).^2 + vTemp(:, 2).^2 + vTemp(:, 3).^2); % Speed vector: one value per persistence frame index
                    speedTemp = sqrt(sum(vTemp.^2, 1))';
                    roundOrNot = true;
                    interpPts = ULM_interp3D_linear(coordsStart, coordsEnd, speedTemp, roundOrNot); % Get interpolated points with the corresponding z velocity value. each row is [z coord, x coord, z velocity]
    
    %                 figure; scatter3(interpPts(:, 1), interpPts(:, 2), interpPts(:, 3))
                    for ipi = 1:size(interpPts, 1) % interpolated point index
                        interpPtsTemp = interpPts(ipi, 1:3);
                        speedValTemp = interpPts(ipi, 4);
                        
    %                     speedMap(smoothedPtsTemp(1), smoothedPtsTemp(2), smoothedPtsTemp(3)) = speedMap(smoothedPtsTemp(1), smoothedPtsTemp(2), smoothedPtsTemp(3)) + speedValTemp;
                        densityMapInterpolatedCounter(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) = densityMapInterpolatedCounter(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) + 1;
                    end
                end
                
%                 if any(totalDist > pixel)
    
                
    
    %             for ipi = 1:size(interpPts, 1) % interpolated point index
    %                 interpPtsTemp = interpPts(ipi, :);
    %                 speedValTemp = interpPtsTemp(4);
    %                 
    %                 speedMap(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) = speedMap(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) + speedValTemp;
    %                 speedMapCounter(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) = speedMapCounter(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) + 1;
    %             end
            end
    %     else
    %         interpPts = [];
    %         coordsStart = [];
    %         coordsEnd = [];
    %         zvTemp = [];
    %         bvTemp = [];
    %         zVelTemp = [];
    %         interpPtsTemp = [];
        end
    %     disp(strcat("Frame ", num2str(ti), " stored."))
    %     clear bvTemp
        
        %%%%%%%%%%%%%%%%
        if mod(ti, 100) == 0 % only get a video frame every N frames
            disp(ti)
            V = volshow(densityMapInterpolatedCounter .^ 0.5, 'Renderer', 'MaximumIntensityProjection');

%             V.Alphamap(1:100) = 0;          % Change transparency
            V.BackgroundColor = [1, 1, 1];  % Make background white
            V.ScaleFactors(3) = size(densityMapInterpolated, 1) / size(densityMapInterpolated, 3) * actualSize(1) / actualSize(3); % scale with # pixels and region size
        %     V.CameraPosition = V.CameraPosition ./ 2;
%             V.CameraPosition = [2.1161 -3.7332 -0.1764];
%             V.CameraUpVector = [0.1853 -0.3160 -0.9305];
            V.CameraUpVector = [0, 0, -1]; % Flip the z axis
            V.CameraViewAngle = 15;
            V.ScaleFactors(3) = size(densityMapInterpolated, 1) / size(densityMapInterpolated, 3) * actualSize(1) / actualSize(3); % scale with # pixels and region size
            V.ScaleFactors = V.ScaleFactors .* 1.5; % zoom

            cv = getframe(vf);     % get the current volume
            rgb = frame2im(cv);      % convert the frame to rgb data
            writeVideo(vo, rgb);
        end
    end
    disp('Density map interpolated')
    close(vo)
    toc
%     clear speedValTemp ti bpi pfi ipi interpPtsTemp speedTemp coordsStart coordsEnd bvTemp
    densityMapInterpolated = densityMapInterpolatedCounter;
    % Take the average for pixels with overlapping tracks
%     speedMask = densityMapInterpolatedCounter > 0;
%     speedMap(speedMask) = speedMap(speedMask) ./ densityMapInterpolatedCounter(speedMask);
end

function plotSpeedMIPs(data, gamcp) % expects 4D input (x, y, z, frames)
    % gamcp = gamma compression power
    
    figure; imagesc(squeeze(max(data, [], 1))' .^ gamcp); colormap jet; colorbar
    title('xz MIP')
    xlabel('y pixels')
    ylabel('z pixels')

    figure; imagesc(squeeze(max(data, [], 2))' .^ gamcp); colormap jet; colorbar
    title('yz MIP')
    xlabel('x pixels')
    ylabel('z pixels')

    figure; imagesc(squeeze(max(data, [], 3))' .^ gamcp); colormap jet; colorbar
    title('xy MIP')
    xlabel('x pixels')
    ylabel('y pixels')

end

function [Tmap] = thresholdMaps(map, counter, lowerCutoff, upperCutoff) % threshold a bubble density map or speed map to remove low and high counts (noise and/or false positives)
    Tmap = map;
    Tmap(counter <= lowerCutoff) = 0;
    Tmap(counter >= upperCutoff) = 0;
end

function track3DPlot(tracksV, indices) % Plot the individual combined tracks in 3D space, for some range of indices
    figure; hold on
    for ind = indices
    %     plot3(tracksIndividualCombined{ test_ind }(:, 2), tracksIndividualCombined{ test_ind }(:, 3), tracksIndividualCombined{ test_ind }(:, 4), '-o')
        plot3([tracksV{ind}(1, 3); tracksV{ind}(:, 6)], [tracksV{ind}(1, 4); tracksV{ind}(:, 7)], [tracksV{ind}(1, 5); tracksV{ind}(:, 8)], '.')
    end
    hold off
end

function trackSpeedSpaghettiPlot(tracksV, indices)
%     figure; title("Speed along tracks " + num2str(min(indices)) + " to " num2str(max(indices)))
    figure; title("Speed along some tracks")
    hold on
    for ind = indices
        temp_v = tracksV{ind}(:, 9:11);
        plot(sqrt(sum(temp_v .^ 2, 2)), '.-')
    end
    hold off
end

function tracksVS_KF_MMS = applyKF(tracksVS_MMS, vMMStoPixelDispPerFrame, img_size, pers)
    tic

    tracksVS_KF_MMS = tracksVS_MMS;
    
    numDims = 3;
    numStates = numDims * 2; % In 3D: x position, y position, z position, x displacement, y displacement, z displacement
    
    if pers < 3
        error('The Kalman filter needs at least 3 points in the track')
    end
    
%     vMMStoPixelDispPerFrame = timePerFrame / 1e3 * [xpixelsPerM, ypixelsPerM, zpixelsPerM];
    
    % Define the matrices that map the state transition, and the state ->
    % observation transformation
    Fk = [1, 0, 0, 1, 0, 0; ...
          0, 1, 0, 0, 1, 0; ...
          0, 0, 1, 0, 0, 1; ...
          0, 0, 0, 1, 0, 0; ...
          0, 0, 0, 0, 1, 0; ...
          0, 0, 0, 0, 0, 1];
    Hk = [1, 0, 0, 0, 0, 0; ...
          0, 1, 0, 0, 0, 0; ...
          0, 0, 1, 0, 0, 0; ...
          0, 0, 0, 1, 0, 0; ...
          0, 0, 0, 0, 1, 0; ...
          0, 0, 0, 0, 0, 1];
    
    %%%%%%%% these covariance matrices are from the Song et al. 2020 paper %%%%%%%%%
    Qk = diag(ones(numStates, 1)) .* 0.5;      % Covariance matrix of the system/process noise
    Rk = diag(ones(numStates, 1)) .* 4;   % Covariance matrix of the observation noise
    
    for ti = 1:length(tracksVS_MMS) % Go through each track
        trackTemp = tracksVS_MMS{ti};
        trackTemp = [trackTemp; [0, 0, trackTemp(end, 6:8), 0, 0, 0, 0, 0, 0]]; % Append the last position onto the column the KF looks at, and fill the rest with zeros (not used except velocities)
    
        % Initialize variables for the state vector and covariance matrix
        xk = NaN(numStates, size(trackTemp, 1)); % xk has dimensions (6 from xyz pos and velocity, # frames in the track)
        Pk = NaN(numStates, numStates, size(trackTemp, 1));
    
        % Initial values               
        xk(:, 1) = [trackTemp(1, 3:5), trackTemp(1, 9:11) .* vMMStoPixelDispPerFrame]'; % Initial state vector
    
        Pk(:, :, 1) = [1, 0, 0, 0, 0, 0; ... % Initial covariance matrix
                       0, 1, 0, 0, 0, 0; ...
                       0, 0, 1, 0, 0, 0; ...
                       0, 0, 0, 10, 0, 0; ...
                       0, 0, 0, 0, 10, 0; ...
                       0, 0, 0, 0, 0, 10];
    
        % Go through each step of the track
        for k = 2:size(trackTemp, 1)
            % Prediction
            xkp = Fk * xk(:, k - 1); 
            Pkp = Fk * Pk(:, :, k - 1) * Fk + Qk;
    
            % Observation
            yk = [trackTemp(k, 3:5), trackTemp(k, 9:11) .* vMMStoPixelDispPerFrame]';
    
            % Update
            Kku = Pkp * Hk' * inv(Hk * Pkp * Hk' + Rk); % Kalman gain matrix
            Iku = yk - Hk * xkp; % Innovation vector (difference between the observed state and the predicted state transformed into an observation at step k)
            xku = xkp + Kku * Iku; % Updated (weighted) estimate for the state at step k
            Pku = (eye(length(xku)) - Kku * Hk) * Pkp; % Updated covariance matrix at step k
    
            % Store the updated state and covariance
            xk(:, k) = xku;
            Pk(:, :, k) = Pku;
    
            % Adjustments in case the Kalman filter puts some points outside the original region
            xktemp = xk(1:3, :);
            xktemp(xktemp < 1) = 1;
            xk(1:3, :) = xktemp;
            xk(1, xk(1, :) > img_size(1)) = img_size(1);
            xk(2, xk(2, :) > img_size(2)) = img_size(2);
            xk(3, xk(3, :) > img_size(3)) = img_size(3);
    
        end
        % Store the post-KF values
        xkt = xk'; % Temporarily store the transpose of xk
        tracksVS_KF_MMS{ti}(:, 3:5) = round(xkt(1:end-1, 1:3)); % "Start" positions
        tracksVS_KF_MMS{ti}(:, 6:8) = round(xkt(2:end, 1:3)); % "End" positions
    %     tracksVS_KF_MMS{ti}(:, 9:11) = xkt(2:end, 4:6) ./ vMMStoPixelDispPerFrame;
        tracksVS_KF_MMS{ti}(:, 9:11) = xkt(1:end - 1, 4:6) ./ vMMStoPixelDispPerFrame; % Velocities
    end
    
    toc
    disp('Kalman filter applied')
end

function [tracks_filtered] = refineLowSpeedTracks(tracksVS_MMS, low_prctile, low_speed_threshold, high_speed_threshold)

    % Make a new variable that keeps only tracks with some velocity characteristics
%     low_prctile = 25;
%     low_speed_threshold = 2; % mm/s
%     
%     high_speed_threshold = 50;

    tracksVS_MMS_thresholded_mask = false(size(tracksVS_MMS));
    tic
    for ti = 1:length(tracksVS_MMS)     % track index - go through each track
        trackTemp = tracksVS_MMS{ti};
        speedsTemp = sqrt(sum(trackTemp(:, 9:11) .^ 2, 2));
        
        speed_low_prctile = prctile(speedsTemp, low_prctile);
        if speed_low_prctile < low_speed_threshold & max(speedsTemp) <= high_speed_threshold
            tracksVS_MMS_thresholded_mask(ti) = true;
        end
    
    end
    toc
    tracks_filtered = tracksVS_MMS(tracksVS_MMS_thresholded_mask);


end