% Description: Part 2 of the ULM processing for a 2D array. This script
%   takes the result of the bubble localization from 'RCA_analyze_ULM.m' ('centers-filenum.mat') and
%   performs pairing and tracking, plus further image refinement for
%   plotting purposes.

%   centers is a a logical matrix of dimensions (# refined x pixels, # refined y
%   pixels, # refined z pixels, # frames per file) with ones at the pixels where a bubble
%   center was localized. There is one 'centers' for each file that the
%   localization was performed on.

%   This script turns that logical matrix into a cell array with their
%   corresponding coordinates, then does pairing with the Hungarian method
%   (assignmunkres). Tracks are created with a frame persistence condition (only keep a
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
%   Colormaps_fUS (From Jianbo's code)
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
load([datapath, 'proc_params.mat'])

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
totalFrames = numFiles * P.numFramesPerBuffer; % Total frames to process

% Cell array with an entry for each frame. Each entry contains (# bubbles) of coordinate pairs (z, x) of the detected bubble centers
% centerCoords = cell(totalFrames, 1); 
centerCoords = {};

% Concatenate all the centers- files
for n = startFile:endFile   % Go through each center file (for each buffer)
% for n = startFile
    tic
    load([datapath, 'centers-', num2str(n)])
    centerCoords = [centerCoords; centersRC];
    disp(strcat("Center coordinates for file ", num2str(n), " stored."))
    toc
end

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

% for cci = 1:length(centerCoords_corrected) % centerCoords index
%     cc = centerCoords_corrected{cci};
% %     cc = cc + pad_dims;
% %     cc = round(cc); %%%%%%%%% FOR TESTING %%%%%%%%%
%     for nbcci = 1:size(cc, 1) % # bubbles in centerCoords_corrected at index cci
%         bubbleDensityMapRaw(cc(nbcci, 1), cc(nbcci, 2), cc(nbcci, 3)) = bubbleDensityMapRaw(cc(nbcci, 1), cc(nbcci, 2), cc(nbcci, 3)) + 1;
%     end
% end
% 
% % volumeViewer(bubbleDensityMapRaw .^ 0.5)
% % figure; imagesc(squeeze(sum(bubbleDensityMapRaw, 1))' .^ 0.5); colormap hot; title('Raw bubble density, sum across y'); colorbar
% figure; imagesc(squeeze(max(bubbleDensityMapRaw, [], 1))' .^ 0.5); colormap hot; title('Raw bubble density, MIP across y'); colorbar
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

% %% 5. Calculate pairing
% bubblePairs = cell(totalFrames - 1, 1);   % Initialize cell vector of paired bubble indices
% 
% ubS = cell(totalFrames - 1, 1);             % unassigned bubbles from the source frames
% ubT = cell(totalFrames - 1, 1);             % unassigned bubbles from the target frames
% 
% tic
% % parfor f = startFrame:totalFrames - 1       % Go through frames
% parfor f = startFrame:startFrame+50
% % parfor f = startFrame:5000
% % for f = startFrame:startFrame+0
%     sourceFrame = centerCoords_corrected{f};     % Get the coordinates for the source frame (f)
%     targetFrame = centerCoords_corrected{f + 1}; % Get the coordinates for the target frame (f + 1)
% 
%     nbS = size(sourceFrame, 1);             % number of bubbles in the source frame
%     nbT = size(targetFrame, 1);             % number of bubbles in the target frame
%     D = NaN(nbS, nbT);                      % Initialize the distance matrix comparing source and target points' distances
% 
%     if nbS > 1 && nbT > 1                   % Only go through the pairing if the source and target frames have at least one bubble
%         % Go through each source frame point and get the distance to each
%         % point in the target frame. Store in the distance matrix D.
%         for spi = 1:nbS % source point index
%             sourcePoint = sourceFrame(spi, :);      % x, y, and z coords of the source frame's point "spi"
%             d = targetFrame - sourcePoint;          % vectorized difference between the x, y, and z coords of all the points in the target frame and point spi from the source frame
%             d(:, 1) = d(:, 1) ./ xpixelsPerM;       % Convert the x distance differences into natural distance units [m]
%             d(:, 2) = d(:, 2) ./ ypixelsPerM;       % Convert the y distance differences into natural distance units [m]
%             d(:, 3) = d(:, 3) ./ zpixelsPerM;       % Convert the z distance differences into natural distance units [m]
%             D(spi, :) = sqrt(sum((d .^ 2), 2));     % distance formula on the above
%     
%         end
%         
%         D(D > maxDistPerFrameM) = Inf;      % Set the elements above the distance per frame threshold to Inf so they aren't considered for pairing
%         [assignment, unassignedrows, unassignedcolumns] = assignmunkres(D, 100); % Pair with the Munkres algorithm, which minimizes the total cost (total paired distance)
%           bubblePairs{f} = assignment; % assignment will always be sorted to make the second column in order
% 
%         % Use Yi Cao's munkres for efficiency, which might break if there
%         % are Infs?
% %         [assignment, cost] = munkres(D);
% %         indTemp = find(assignment);
% %         [sourceInd, targetInd] = ind2sub(size(assignment), indTemp);
% %         bubblePairs{f} = [sourceInd, targetInd]; % assignment will always be sorted to make the second column in order
% 
% %         ubS{f} = unassignedrows;
% %         ubT{f} = unassignedcolumns;
%         
%     end
% end
% toc
% disp('Pairing done')
% 
% clear f nbS nbT D spi assignment unassignedrows unassignedcolumns

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

%% 6. Create tracks with persistence

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

%% 8. Create the velocity map

% bVelocityC = cell(size(tracksClean, 1), pers); % Optional cell array storage of the bubble data
bVelocityM = cell(size(tracksClean, 1), 1); % Matrix storage of the bubble data. Each element in the cell is a [# bubbles per frame in the track, 9, # persistence frames] matrix. Each row corresponds to [bubble f x coord, bubble f y coord, bubble f z coord, bubble f+1 x coord, bubble f+1 y coord, bubble f+1 z coord, x velocity, y velocity, z velocity]
tic
for ti = startFrame:length(tracksClean)     % track index - go through each track
    tracksTemp = tracksClean{ti};           % Get track ti
    nbiti = nbitAll(ti);                    % # of bubbles in the tracks starting in index ti
    if nbiti > 0
        for fn = 1:pers                     % Go through all the frames in the tracks with origin frame ti
            startPoints = tracksTemp((fn - 1) * nbiti + 1 : fn * nbiti, 2:4);
            endPoints = tracksTemp((fn) * nbiti + 1 : (fn + 1) * nbiti, 2:4);
            vfn = (endPoints - startPoints) ./ timePerFrame; % velocity = displacement/time
%             bVelocityC{ti, fn} = [startPoints, endPoints, vfn];  % each row is [x start coord, y start coord, z start coord, x end coord, y end coord, z end coord, x velocity, y velocity, z velocity]
            bVelocityM{ti}(:, :, fn) = [startPoints, endPoints, vfn];
        end
    end
end
toc
disp('Velocity map created')
clear ti fn tracksTemp startPoints endPoints vfn nbiti

%% 9. Refine the velocity map
bVelocityMSmoothed = bVelocityM;          % Initialize the velocity data, which will be smoothed across frames with a moving mean

for n = startFrame:size(bVelocityM, 1) % Go through each track collection n
    vt = bVelocityM{n};
    vtSmoothed = vt;        % Temporary variable with the smoothed velocities for each track collection n
    for tn = 1:size(vt, 1)  % track number (go through the path of each bubble over time)
        tnxVel = squeeze(vt(tn, 7, :));         % get the x velocity at each point in track tn
        tnxVelSmoothed = movmean(tnxVel, mmws); % moving mean
        vtSmoothed(tn, 7, :) = tnxVelSmoothed;  % Store the smoothed velocity

        tnyVel = squeeze(vt(tn, 8, :));         % get the y velocity at each point in track tn
        tnyVelSmoothed = movmean(tnyVel, mmws);
        vtSmoothed(tn, 8, :) = tnyVelSmoothed;

        tnzVel = squeeze(vt(tn, 9, :));         % get the z velocity at each point in track tn
        tnzVelSmoothed = movmean(tnzVel, mmws);
        vtSmoothed(tn, 9, :) = tnzVelSmoothed;
        
    end
    bVelocityMSmoothed{n} = vtSmoothed; % Store the smoothed velocity data for track collection n
end
clear n tn vt vtSmoothed tnzVel tnzVelSmoothed vtSmoothed tnxVel tnxVelSmoothed tnyVel tnyVelSmoothed 

% Scale the smoothed velocity into [mm/s]
bVelocityMSmoothedMMS = bVelocityMSmoothed;
for n = startFrame:size(bVelocityMSmoothedMMS, 1)
    if ~isempty(bVelocityMSmoothedMMS{n})
        bVelocityMSmoothedMMS{n}(:, 7, :) = bVelocityMSmoothed{n}(:, 7, :) ./ xpixelsPerM * 1e3;
        bVelocityMSmoothedMMS{n}(:, 8, :) = bVelocityMSmoothed{n}(:, 8, :) ./ ypixelsPerM * 1e3;
        bVelocityMSmoothedMMS{n}(:, 9, :) = bVelocityMSmoothed{n}(:, 9, :) ./ zpixelsPerM * 1e3;
    end
end

disp('Velocity map smoothed')

%% 10a. Kalman filter with velocity in the state
numDims = 3;
numStates = numDims * 2; % In 3D: x position, y position, z position, x displacement, y displacement, z displacement

if pers < 3
    error('The Kalman filter needs at least 3 points in the track')
end

% Initialize the filtered output variable
bVelocityMSmoothedKFMMS = bVelocityMSmoothedMMS;
% vMMStoPixelDispPerFrame = timePerFrame / 1e3 * pixelsPerM;
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
for n = startFrame:size(bVelocityMSmoothedKFMMS, 1)
% for n = 7961
    tln = bVelocityMSmoothedKFMMS{n}; % Track list n
    %%%% MAYBE MOVE THE BELOW LOOP INTO A FUNCTION WITH A PARFOR %%%
    for tn = 1:size(tln, 1) % Track number tn
        track = squeeze(tln(tn, :, :))'; % Get the track
        track = [track; [track(end, 4:6), 0, 0, 0, 0, 0, 0]];

        % Initialize variables for the state vector and covariance matrix
        xk = NaN(numStates, size(track, 1));
        Pk = NaN(numStates, numStates, size(track, 1));

        % Initial values                 
        xk(:, 1) = [track(1, 1:3), track(1, 7:9) .* vMMStoPixelDispPerFrame]'; % Initial state vector

        Pk(:, :, 1) = [1, 0, 0, 0, 0, 0; ... % Initial covariance matrix
                       0, 1, 0, 0, 0, 0; ...
                       0, 0, 1, 0, 0, 0; ...
                       0, 0, 0, 10, 0, 0; ...
                       0, 0, 0, 0, 10, 0; ...
                       0, 0, 0, 0, 0, 10];
        %%%%%% WHAT INITIAL COVARIANCES SHOULD I USE??? %%%%%%

%         xk(:, 2) = [track(2, 1:2), track(1, 5:6) .* vMMStoPixelDispPerFrame]';
% %         xk(:, 2) = [track(2, 1:2), track(1, 5:6)]';
%         Pk(:, :, 2) = [1, 0, 0, 0; ...
%                  0, 1, 0, 0; ...
%                  0, 0, 10^4, 0; ...
%                  0, 0, 0, 10^4];

        for k = 2:size(track, 1) % Go through each step of the track
            % Prediction
            xkp = Fk * xk(:, k - 1); 
            Pkp = Fk * Pk(:, :, k - 1) * Fk + Qk;

            % Observation
            yk = [track(k, 1:3), track(k, 7:9) .* vMMStoPixelDispPerFrame]';

            % Update
            Kku = Pkp * Hk' * inv(Hk * Pkp * Hk' + Rk); % Kalman gain matrix
            Iku = yk - Hk * xkp; % Innovation vector (difference between the observed state and the predicted state transformed into an observation at step k)
            xku = xkp + Kku * Iku; % Updated (weighted) estimate for the state at step k
            Pku = (eye(length(xku)) - Kku * Hk) * Pkp; % Updated covariance matrix at step k

            % Store the updated state and covariance
            xk(:, k) = xku;
            Pk(:, :, k) = Pku;
        end

        % Adjustments in case the Kalman filter puts some points outside the original region
        xktemp = xk(1:3, :);
        xktemp(xktemp < 1) = 1;
        xk(1:3, :) = xktemp;
        xk(1, xk(1, :) > img_size(1)) = img_size(1);
        xk(2, xk(2, :) > img_size(2)) = img_size(2);
        xk(3, xk(3, :) > img_size(3)) = img_size(3);
        
        bVelocityMSmoothedKFMMS{n}(tn, 1:3, :) = round(xk(1:3, 1:end-1));
        bVelocityMSmoothedKFMMS{n}(tn, 4:6, :) = round(xk(1:3, 2:end));
%         bVelocityMSmoothedKFMMS{n}(tn, 7:9, :) = xk(4:6, 2:end) ./ vMMStoPixelDispPerFrame;
        bVelocityMSmoothedKFMMS{n}(tn, 7:9, :) = xk(4:6, 2:end) ./ repmat(vMMStoPixelDispPerFrame', 1, pers);
%         bVelocityTestMSmoothedKFMMS{n}(tn, 5:6, :) = (xk(1:2, 2:end) - xk(1:2, 1:end-1)); % ./ timePerFrame;

        % plot test to compare for a single track
%         figure
%         hold on
%         plot(xk(1, :), xk(2, :), '-o')
%         plot(track(:, 1), track(:, 2), '--x')
%         clear i xki
%         legend('Kalman filtered', 'Original')
%         hold off
    end
end
toc
disp('Kalman filter applied')
clear n tn tln k xk Pk yk Kku Iku xku Pku track

%% 10b. Kalman filter without velocity in the observation
numDims = 3;
numStates = numDims * 2; % In 3D: x position, y position, z position, x displacement, y displacement, z displacement

if pers < 3
    error('The Kalman filter needs at least 3 points in the track')
end

% Initialize the filtered output variable
bVelocityMSmoothedKFMMS = bVelocityMSmoothedMMS;
% vMMStoPixelDispPerFrame = timePerFrame / 1e3 * pixelsPerM;
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
      0, 0, 1, 0, 0, 0];

%%%%%%%% these covariance matrices are from the Song et al. 2020 paper %%%%%%%%%
Qk = diag(ones(numStates, 1)) .* 0.5;      % Covariance matrix of the system/process noise
Rk = diag(ones(numStates, 1)) .* 4;   % Covariance matrix of the observation noise

tic
for n = startFrame:size(bVelocityMSmoothedKFMMS, 1)
% for n = 7961
    tln = bVelocityMSmoothedKFMMS{n}; % Track list n
    %%%% MAYBE MOVE THE BELOW LOOP INTO A FUNCTION WITH A PARFOR %%%
    for tn = 1:size(tln, 1) % Track number tn
        track = squeeze(tln(tn, :, :))'; % Get the track
        track = [track; [track(end, 4:6), 0, 0, 0, 0, 0, 0]];

        % Initialize variables for the state vector and covariance matrix
        xk = NaN(numStates, size(track, 1));
        Pk = NaN(numStates, numStates, size(track, 1));

        % Initial values                 
        xk(:, 1) = [track(1, 1:3), track(1, 7:9) .* vMMStoPixelDispPerFrame]'; % Initial state vector

        Pk(:, :, 1) = [1, 0, 0, 0, 0, 0; ... % Initial covariance matrix
                       0, 1, 0, 0, 0, 0; ...
                       0, 0, 1, 0, 0, 0; ...
                       0, 0, 0, 10, 0, 0; ...
                       0, 0, 0, 0, 10, 0; ...
                       0, 0, 0, 0, 0, 10];
        %%%%%% WHAT INITIAL COVARIANCES SHOULD I USE??? %%%%%%

%         xk(:, 2) = [track(2, 1:2), track(1, 5:6) .* vMMStoPixelDispPerFrame]';
% %         xk(:, 2) = [track(2, 1:2), track(1, 5:6)]';
%         Pk(:, :, 2) = [1, 0, 0, 0; ...
%                  0, 1, 0, 0; ...
%                  0, 0, 10^4, 0; ...
%                  0, 0, 0, 10^4];

        for k = 2:size(track, 1) % Go through each step of the track
            % Prediction
            xkp = Fk * xk(:, k - 1); 
            Pkp = Fk * Pk(:, :, k - 1) * Fk + Qk;

            % Observation
            yk = [track(k, 1:3), track(k, 7:9) .* vMMStoPixelDispPerFrame]';

            % Update
            Kku = Pkp * Hk' * inv(Hk * Pkp * Hk' + Rk); % Kalman gain matrix
            Iku = yk - Hk * xkp; % Innovation vector (difference between the observed state and the predicted state transformed into an observation at step k)
            xku = xkp + Kku * Iku; % Updated (weighted) estimate for the state at step k
            Pku = (eye(length(xku)) - Kku * Hk) * Pkp; % Updated covariance matrix at step k

            % Store the updated state and covariance
            xk(:, k) = xku;
            Pk(:, :, k) = Pku;
        end

        % Adjustments in case the Kalman filter puts some points outside the original region
        xktemp = xk(1:3, :);
        xktemp(xktemp < 1) = 1;
        xk(1:3, :) = xktemp;
        xk(1, xk(1, :) > img_size(1)) = img_size(1);
        xk(2, xk(2, :) > img_size(2)) = img_size(2);
        xk(3, xk(3, :) > img_size(3)) = img_size(3);
        
        bVelocityMSmoothedKFMMS{n}(tn, 1:3, :) = round(xk(1:3, 1:end-1));
        bVelocityMSmoothedKFMMS{n}(tn, 4:6, :) = round(xk(1:3, 2:end));
%         bVelocityMSmoothedKFMMS{n}(tn, 7:9, :) = xk(4:6, 2:end) ./ vMMStoPixelDispPerFrame;
        bVelocityMSmoothedKFMMS{n}(tn, 7:9, :) = xk(4:6, 2:end) ./ repmat(vMMStoPixelDispPerFrame', 1, pers);
%         bVelocityTestMSmoothedKFMMS{n}(tn, 5:6, :) = (xk(1:2, 2:end) - xk(1:2, 1:end-1)); % ./ timePerFrame;

        % plot test to compare for a single track
%         figure
%         hold on
%         plot(xk(1, :), xk(2, :), '-o')
%         plot(track(:, 1), track(:, 2), '--x')
%         clear i xki
%         legend('Kalman filtered', 'Original')
%         hold off
    end
end
toc
disp('Kalman filter applied')
clear n tn tln k xk Pk yk Kku Iku xku Pku track

%% 11. Acceleration and direction constraints

tic
% bVelocityConstrained = applyConstraints(bVelocityM, vTrimmedMeanPercentage, aThresholdFactor, angleChangeThreshold, timePerFrame);

% bVelocityMSmoothedMMSConstrained = applyConstraints(bVelocityMSmoothedMMS, vTrimmedMeanPercentage, aThresholdFactor, angleChangeThreshold, timePerFrame);

bVelocityMSmoothedKFConstrainedMMS = applyConstraints(bVelocityMSmoothedKFMMS, vTrimmedMeanPercentage, aThresholdFactor, angleChangeThreshold, timePerFrame);
toc
disp('Acceleration and direction constraints applied')
clear n tln tn trackAlreadyDeleted track vTrack vTrackTrimmedMean aThresholdMag aTrackMag angleTrack angelTrackChanges

%% 12. Plot bubble density map(s) with the paired bubbles after persistence

[BDM] = densityMap3D(bVelocityM, img_size, startFrame);
[BDM_Constrained] = densityMap3D(bVelocityConstrained, img_size, startFrame);
[BDM_SmoothedMMS] = densityMap3D(bVelocityMSmoothedMMS, img_size, startFrame);

[BDM_SmoothedKF] = densityMap3D(bVelocityMSmoothedKFMMS, img_size, startFrame);
[BDM_SmoothedKFConstrained] = densityMap3D(bVelocityMSmoothedKFConstrainedMMS, img_size, startFrame);

clear n tn iti trackTemp tempBuf

% volumeViewer(BDM .^ 0.3)
%%
BDM_Constrained_Rfn = thresholdMaps(BDM_Constrained, BDM_Constrained, 2, 500);
volumeViewer(BDM_Constrained .^ 0.3)
%%
volumeViewer(BDM_SmoothedMMS .^ 0.3)

% volumeViewer(bSumSmoothedKF .^ 0.5)
volumeViewer(BDM_SmoothedKFConstrained .^ 0.3)

% imgRefinementFactor_map = [2, 2, 2];
% refineMap = @(map, imgRefinementFactor_map) imresize3(map, [size(map, 1) * imgRefinementFactor_map(1), size(map, 2) * imgRefinementFactor_map(2), size(map, 3) * imgRefinementFactor_map(3)]);
% bSumRefined = refineMap(bSum, imgRefinementFactor_map);
% volumeViewer(bSumRefined .^ 1)

% yrange_plot_MIP = 140:180;
% figure; imagesc(abs(squeeze(max(bSumRefined(yrange_plot_MIP, :, :), [], 1) .^ 0.3)')); colormap hot; colorbar
% title("bSumRefined Maximum Intensity Projection from y = " + num2str(yrange_plot_MIP(1)) + " to " + num2str(yrange_plot_MIP(end)) + " \^ 0.3")

% figure; imagesc(sum(bSumSmoothedKFConstrained .^ 0.5, 3))
% holeMaskThreshold = 10;
% maskHole = sum(bSumSmoothedKFConstrained .^ 0.5, 3) < holeMaskThreshold;
% figure; imagesc(maskHole)

% Set voxels with a bubble count of 2 or less = 0
% bSumSmoothedKFConstrainedTest = bSumSmoothedKFConstrained;
% bSumSmoothedKFConstrainedTest(bSumSmoothedKFConstrainedTest <= 2) = 0;
% volumeViewer(bSumSmoothedKFConstrainedTest .^ 0.3)

% 2D sum plots
% plotPower2D = 0.3;
% figure; imagesc(squeeze(sum(BDM, 1).^ plotPower2D)'); colormap hot; title('BDM sum across y'); colorbar
% figure; imagesc(squeeze(sum(bSum, 2).^ plotPower2D)'); colormap hot; title('bSum sum across x'); colorbar
% figure; imagesc(squeeze(sum(bSum, 3).^ plotPower2D)'); colormap hot; title('bSum sum across z'); colorbar

% 2D selected slice plots or max intensity projection from a small range of
% slices
% figure; imagesc(squeeze(bSum(80, :, :) .^ 0.5)'); colormap hot
% figure; imagesc(squeeze(max(bSum, [], 1).^ 0.5)'); colormap hot; title('bSum MIP across y')
% yrange_plot_MIP = 70:90;
% yrange_plot_MIP = 45:55;
yrange_plot_MIP = 75:85;
figure; imagesc(squeeze(max(BDM(yrange_plot_MIP, :, :), [], 1) .^ 0.3)'); colormap hot; colorbar
title("bSum Maximum Intensity Projection from y = " + num2str(yrange_plot_MIP(1)) + " to " + num2str(yrange_plot_MIP(end)) + " \^ 0.3")

xrange_plot_MIP = 70:90;
figure; imagesc(squeeze(max(BDM(:, xrange_plot_MIP, :), [], 2) .^ 0.3)'); colormap hot; colorbar
title("bSum Maximum Intensity Projection from x = " + num2str(xrange_plot_MIP(1)) + " to " + num2str(xrange_plot_MIP(end)) + " \^ 0.3")

zrange_plot_MIP = 100:150;
figure; imagesc(squeeze(max(BDM(:, :, zrange_plot_MIP), [], 3) .^ 0.5)'); colormap hot; colorbar
title("bSum Maximum Intensity Projection from z = " + num2str(zrange_plot_MIP(1)) + " to " + num2str(zrange_plot_MIP(end)) + " \^ 0.3")

% figure; imagesc(squeeze(max(bSumSmoothedKFConstrainedTest(yrange_plot_MIP, :, :), [], 1) .^ 0.3)'); colormap hot
% title("bSumSmoothedKFConstrained count<=2 removed Maximum Intensity Projection from y = " + num2str(yrange_plot_MIP(1)) + " to " + num2str(yrange_plot_MIP(end)) + " \^0.3")

% addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Processing\Jerman Enhancement Filter\')
% test = vesselness3D(bSum .^ 0.5, 1:5, [xpix_spacing; ypix_spacing; zpix_spacing], 0.5, true);
% volumeViewer(test)
%%
% test_dmi = interpolatedDensityMap(bVelocityMSmoothedKFConstrainedMMS, img_size, startFrame); % test density map interpolated
% test_dmi = interpolatedDensityMap(bVelocityM, img_size, startFrame); % test density map interpolated

BDM_LI = interpolatedDensityMap(bVelocityM, img_size, startFrame, maxPixelDistPerFrame); % density map, linearly interpolated
% BDM_LI = interpolatedDensityMap(bVelocityM, img_size, startFrame, maxPixelDistPerFrame ./ 2); % density map, linearly interpolated
% Probably should only do this if we're confident about the max speed input
% volumeViewer(test_dmi .^ 0.3)
% volumeViewer(BDM_LI .^ 0.3)
%%
plotMIPs(BDM_LI , 0.4)
%%
% BDM_SmoothedKFConstrained_LI = interpolatedDensityMap(bVelocityMSmoothedKFConstrainedMMS, img_size, startFrame, maxPixelDistPerFrame); % density map, linearly interpolated
BDM_SmoothedKFConstrained_LI_Rfn = BDM_SmoothedKFConstrained_LI;
BDM_SmoothedKFConstrained_LI_Rfn = thresholdMaps(BDM_SmoothedKFConstrained_LI_Rfn, BDM_SmoothedKFConstrained_LI_Rfn, 2, 100);
% volumeViewer(BDM_SmoothedKFConstrained_LI_RSC .^ 0.4)
plotMIPs(BDM_SmoothedKFConstrained_LI_Rfn, 0.7)

generateTiffStack_multi([{BDM_SmoothedKFConstrained_LI_Rfn .^ 0.5}], [8.8, 8.8, 8], 'hot', 5)

%% Plot the interpolated density map with video
actualSize = [lateral_width, lateral_width, axial_depth];
BDM_video = interpolatedDensityMapWithVideo(bVelocityMSmoothedMMSConstrained, img_size, startFrame, maxPixelDistPerFrame, actualSize); % test density map interpolated
% BDM_video = interpolatedDensityMapWithVideo(bVelocityM, img_size, startFrame, maxPixelDistPerFrame, actualSize); % test density map interpolated

%% MIP plots for the interpolation test
% yrange_plot_MIP = 75:85;
yrange_plot_MIP = 70:90;
figure; imagesc(squeeze(max(test_dmi_2(yrange_plot_MIP, :, :), [], 1) .^ 0.3)'); colormap hot; colorbar
% figure; imagesc(squeeze(max(test_dmi(yrange_plot_MIP, :, :), [], 1) .^ 1)'); colormap hot; colorbar
title("test_dmi_2 Maximum Intensity Projection from y = " + num2str(yrange_plot_MIP(1)) + " to " + num2str(yrange_plot_MIP(end)) + " \^ 0.3")

figure; imagesc(squeeze(max(BDM(yrange_plot_MIP, :, :), [], 1) .^ 0.3)'); colormap hot; colorbar
% figure; imagesc(squeeze(max(test_dmi(yrange_plot_MIP, :, :), [], 1) .^ 1)'); colormap hot; colorbar
title("bSum Maximum Intensity Projection from y = " + num2str(yrange_plot_MIP(1)) + " to " + num2str(yrange_plot_MIP(end)) + " \^ 0.3")

xrange_plot_MIP = 75:85;
figure; imagesc(squeeze(max(test_dmi_2(:, xrange_plot_MIP, :), [], 2) .^ 0.3)'); colormap hot; colorbar
title("test_dmi_2 Maximum Intensity Projection from x = " + num2str(xrange_plot_MIP(1)) + " to " + num2str(xrange_plot_MIP(end)) + " \^ 0.3")

% zrange_plot_MIP = 100:150;
zrange_plot_MIP = 110:130;
figure; imagesc(squeeze(max(test_dmi_2(:, :, zrange_plot_MIP), [], 3) .^ 0.5)'); colormap hot; colorbar
title("test_dmi_2 Maximum Intensity Projection from z = " + num2str(zrange_plot_MIP(1)) + " to " + num2str(zrange_plot_MIP(end)) + " \^ 0.3")

figure; imagesc(squeeze(max(BDM(:, :, zrange_plot_MIP), [], 3) .^ 0.4)'); colormap hot; colorbar
title("bSum Maximum Intensity Projection from z = " + num2str(zrange_plot_MIP(1)) + " to " + num2str(zrange_plot_MIP(end)) + " \^ 0.3")

%%
BDM_LI_Rfn = BDM_LI;
BDM_LI_Rfn = thresholdMaps(BDM_LI_Rfn, BDM_LI_Rfn, 2, 300);


volumeViewer(BDM_LI_Rfn .^ 0.4)

%% generate Tiff stacks
% generateTiffStack(BDM_LI_Rfn .^ 1, [lateral_width, lateral_width, axial_depth], 'gray')
generateTiffStack_multi([{BDM_LI_Rfn .^ 0.5}], [8.8, 8.8, 8], 'hot', 10)

%%
yrange_plot_MIP = 70:90;
figure; imagesc(squeeze(max(BDM_LI_Rfn(yrange_plot_MIP, :, :), [], 1) .^ 0.4)'); colormap hot; colorbar
% figure; imagesc(squeeze(max(test_dmi(yrange_plot_MIP, :, :), [], 1) .^ 1)'); colormap hot; colorbar
title("test_dmi_remove_smallcounts Maximum Intensity Projection from y = " + num2str(yrange_plot_MIP(1)) + " to " + num2str(yrange_plot_MIP(end)) + " \^ 0.3")

xrange_plot_MIP = 70:90;
figure; imagesc(squeeze(max(BDM_LI_Rfn(:, xrange_plot_MIP, :), [], 2) .^ 0.3)'); colormap hot; colorbar
title("test_dmi_remove_smallcounts Maximum Intensity Projection from x = " + num2str(xrange_plot_MIP(1)) + " to " + num2str(xrange_plot_MIP(end)) + " \^ 0.3")

zrange_plot_MIP = 100:150;
figure; imagesc(squeeze(max(BDM_LI_Rfn(:, :, zrange_plot_MIP), [], 3) .^ 0.5)'); colormap hot; colorbar
title("test_dmi_remove_smallcounts Maximum Intensity Projection from z = " + num2str(zrange_plot_MIP(1)) + " to " + num2str(zrange_plot_MIP(end)) + " \^ 0.3")

% figure; imagesc(squeeze(test_dmi(80, :, :))' .^ plotPower2D); colormap hot
% figure; imagesc(squeeze(sum(test_dmi, 1))' .^ plotPower2D); colormap hot
% findfigs
% bw = imbinarize(bSum .^ 1);

%% convert outdated names to new ones
% BDM_LI = test_dmi_2;
% clear test_dmi_2
% BDM = bSum; clear bSum
% BDM_Constrained = bSumConstrained; clear bSumConstrained
% BDM_SmoothedMMS = bSumSmoothedMMS; clear bSumSmoothedMMS
% BDM_SmoothedKF = bSumSmoothedKF; clear bSumSmoothedKF
% BDM_SmoothedKFConstrained = bSumSmoothedKFConstrained; clear bSumSmoothedKFConstrained

% Structure to make it easier
% BDMs_AZ02_day3.BDM = BDM;
% BDMs_AZ02_day3.BDM_LI = BDM_LI;
% BDMs_AZ02_day3.BDM_LI_RSC = BDM_LI_RSC;
% BDMs_AZ02_day3.BDM_Constrained = BDM_Constrained;
% BDMs_AZ02_day3.BDM_SmoothedMMS = BDM_SmoothedMMS;
% BDMs_AZ02_day3.BDM_SmoothedKFConstrained = BDM_SmoothedKFConstrained;

% BDMs_AZ02_hour1.BDM = BDM;
% BDMs_AZ02_hour1.BDM_LI = BDM_LI;
% BDMs_AZ02_hour1.BDM_LI_RSC = BDM_LI_RSC;
% BDMs_AZ02_hour1.BDM_Constrained = BDM_Constrained;
% BDMs_AZ02_hour1.BDM_SmoothedMMS = BDM_SmoothedMMS;
% BDMs_AZ02_hour1.BDM_SmoothedKFConstrained = BDM_SmoothedKFConstrained;

% BDMs_AZ02_baseline.BDM = BDM;
% BDMs_AZ02_baseline.BDM_LI = BDM_LI;
% BDMs_AZ02_baseline.BDM_LI_RSC = BDM_LI_RSC;
% BDMs_AZ02_baseline.BDM_Constrained = BDM_Constrained;
% BDMs_AZ02_baseline.BDM_SmoothedMMS = BDM_SmoothedMMS;
% BDMs_AZ02_baseline.BDM_SmoothedKFConstrained = BDM_SmoothedKFConstrained;
% BDMs_AZ02_baseline.BDM_SmoothedKFConstrained_LI_RSC = BDM_SmoothedKFConstrained_LI_RSC;

% BDMs_AZ02_day7.BDM = BDM;
% BDMs_AZ02_day7.BDM_LI = BDM_LI;
% BDMs_AZ02_day7.BDM_LI_RSC = BDM_LI_Rfn;
% BDMs_AZ02_day7.BDM_Constrained = BDM_Constrained;
% BDMs_AZ02_day7.BDM_SmoothedMMS = BDM_SmoothedMMS;
% BDMs_AZ02_day7.BDM_SmoothedKFConstrained = BDM_SmoothedKFConstrained;
% BDMs_AZ02_day7.BDM_SmoothedKFConstrained_LI_RSC = BDM_SmoothedKFConstrained_LI_RSC;

BDMs_AZ03_baseline.BDM = BDM;
BDMs_AZ03_baseline.BDM_LI = BDM_LI;
BDMs_AZ03_baseline.BDM_LI_Rfn = BDM_LI_Rfn;
BDMs_AZ03_baseline.BDM_Constrained = BDM_Constrained;
BDMs_AZ03_baseline.BDM_SmoothedMMS = BDM_SmoothedMMS;
BDMs_AZ03_baseline.BDM_SmoothedKFConstrained = BDM_SmoothedKFConstrained;
BDMs_AZ03_baseline.BDM_SmoothedKFConstrained_LI_RSC = BDM_SmoothedKFConstrained_LI_Rfn;
%% Get the speed maps
% [SM_LI, SM_LI_counter] = interpolatedSpeedMap(bVelocityM, img_size, startFrame, maxPixelDistPerFrame); % flow speed map, linearly interpolated
% % Refine the speed map
% SM_LI_Rfn = SM_LI;
% SM_LI_Rfn = thresholdMaps(SM_LI_Rfn, SM_LI_counter, 2, 300);

% Constrained KF
[SM_SmoothedKFConstrained_LI, SM_SmoothedKFConstrained_LI_counter] = interpolatedSpeedMap(bVelocityMSmoothedKFConstrainedMMS, img_size, startFrame, maxPixelDistPerFrame); % flow speed map, linearly interpolated
SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI;
SM_SmoothedKFConstrained_LI_Rfn = thresholdMaps(SM_SmoothedKFConstrained_LI_Rfn, SM_SmoothedKFConstrained_LI_counter, 2, 300);

% Look at the smoothed constrained no KF data %%%%%%
% [SM_SC_LI, SM_SC_LI_counter] = interpolatedSpeedMap(bVelocityMSmoothedMMSConstrained, img_size, startFrame, maxPixelDistPerFrame); % flow speed map, linearly interpolated
% 
% SM_SC_LI_Rfn = SM_SC_LI;
% SM_SC_LI_Rfn = thresholdMaps(SM_SC_LI_Rfn, SM_SC_LI_counter, 2, 300);

%% Generate a file for other software to read
% writematrix(SM_SmoothedKFConstrained_LI_Rfn, 'D:\Allen\Data\AZ02 Stroke ULM RC15gV\04-22-2025 7d left eye\ULM subpixel processing results\Speed maps\vol.dat')

%% Plot speed map

% volumeViewer(SM_LI_Rfn)
% plotMIPs(SM_LI_RSC, 1)

% volumeViewer(SM_SmoothedKFConstrained_LI_Rfn)

cmap = colormap_ULM;
figure; imagesc(squeeze(max(SM_SmoothedKFConstrained_LI_Rfn(400:600, :, :), [], 1))'); colormap(cmap); clim([0, 40])

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
generateTiffStack_multi([{SM_SmoothedKFConstrained_LI_Rfn}], [8.8, 8.8, 8], cmap, 50, [0, 40])
% generateTiffStack_multi([{test}], [8.8, 8.8, 8], cmap, 50, [0, testlim])

% generateTiffStack_multi([{SM_SmoothedKFConstrained_LI_Rfn}], [8.8, 8.8, 8], cmap, 1)

%% Convert the speed maps to a structure
% SMs_AZ02_baseline.SM_LI = SM_LI;
% SMs_AZ02_baseline.SM_LI_counter = SM_LI_counter;
% SMs_AZ02_baseline.SM_LI_RSC = SM_LI_RSC;
% SMs_AZ02_baseline.SM_SC_LI_counter = SM_SC_LI_counter;
% SMs_AZ02_baseline.SM_SC_LI_RSC = SM_SC_LI_RSC;
% SMs_AZ02_baseline.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ02_baseline.SM_SmoothedKFConstrained_LI_RSC = SM_SmoothedKFConstrained_LI_RSC;
% SMs_AZ02_baseline.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;

% % SMs_AZ02_hour1.SM = SM;
% SMs_AZ02_hour1.SM_LI = SM_LI;
% SMs_AZ02_hour1.SM_LI_counter = SM_LI_counter;
% SMs_AZ02_hour1.SM_LI_RSC = SM_LI_RSC;
% SMs_AZ02_hour1.SM_SC_LI_counter = SM_SC_LI_counter;
% SMs_AZ02_hour1.SM_SC_LI_RSC = SM_SC_LI_RSC;
% SMs_AZ02_hour1.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ02_hour1.SM_SmoothedKFConstrained_LI_RSC = SM_SmoothedKFConstrained_LI_RSC;
% SMs_AZ02_hour1.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;

% SMs_AZ02_day3.SM_LI = SM_LI;
% SMs_AZ02_day3.SM_LI_counter = SM_LI_counter;
% SMs_AZ02_day3.SM_LI_RSC = SM_LI_RSC;
% SMs_AZ02_day3.SM_SC_LI_counter = SM_SC_LI_counter;
% SMs_AZ02_day3.SM_SC_LI_RSC = SM_SC_LI_RSC;
% SMs_AZ02_day3.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ02_day3.SM_SmoothedKFConstrained_LI_RSC = SM_SmoothedKFConstrained_LI_RSC;
% SMs_AZ02_day3.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;

% SMs_AZ02_day7.SM_LI = SM_LI;
% SMs_AZ02_day7.SM_LI_counter = SM_LI_counter;
% SMs_AZ02_day7.SM_LI_RSC = SM_LI_RSC;
% SMs_AZ02_day7.SM_SC_LI_counter = SM_SC_LI_counter;
% SMs_AZ02_day7.SM_SC_LI_RSC = SM_SC_LI_RSC;
% SMs_AZ02_day7.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ02_day7.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
% SMs_AZ02_day7.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;

% SMs_AZ03_baseline.SM_LI = SM_LI;
% SMs_AZ03_baseline.SM_LI_counter = SM_LI_counter;
% SMs_AZ03_baseline.SM_LI_Rfn = SM_LI_Rfn;
% SMs_AZ03_baseline.SM_SC_LI_counter = SM_SC_LI_counter;
% SMs_AZ03_baseline.SM_SC_LI_Rfn = SM_SC_LI_Rfn;
% SMs_AZ03_baseline.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ03_baseline.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
% SMs_AZ03_baseline.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;


% SMs_AZ03_baseline.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
% SMs_AZ03_baseline.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
% SMs_AZ03_baseline.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;

SMs_AZ03_hour1.SM_SmoothedKFConstrained_LI = SM_SmoothedKFConstrained_LI;
SMs_AZ03_hour1.SM_SmoothedKFConstrained_LI_Rfn = SM_SmoothedKFConstrained_LI_Rfn;
SMs_AZ03_hour1.SM_SmoothedKFConstrained_counter = SM_SmoothedKFConstrained_LI_counter;

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

function [speedMapInterpolated, speedMapInterpolatedCounter] = interpolatedSpeedMap(bVelocityM, img_size, startFrame, maxPixelDistPerFrame)
    speedMapInterpolated = zeros(img_size(1), img_size(2), img_size(3));

    % Counters for proper averaging if there are overlapped pixels from
    % different tracks
    speedMapInterpolatedCounter = zeros(size(speedMapInterpolated));
    
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
                    vTemp = squeeze(bvTemp(bpi, 7:9, :)); % Velocity components for the track # bpi
%                 speedTemp = sqrt(vTemp(:, 1).^2 + vTemp(:, 2).^2 + vTemp(:, 3).^2); % Speed vector: one value per persistence frame index
                    speedTemp = sqrt(sum(vTemp.^2, 1))';
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
    %     clear bvTemp
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