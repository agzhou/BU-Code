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

% Acknowledgement: using Jianbo Tang's ULM code, the Song group's ULM papers, and Jean-Yves Tinevez's simpletracker as references
clearvars

%% 1. Add dependencies and load parameters
% Get data path of the localized bubble centers
datapath = uigetdir('F:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\', 'Select the data path');
datapath = [datapath, '\'];

% Load localization processing parameters: proc_params.mat
load([datapath, 'proc_params.mat'])

% Choose and load the params.mat file (from the acquisition)
[params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', 'G:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\params.mat');
load([params_pathname, params_filename])

% Choose the path and add Jianbo's functions to the Matlab path for the Colormaps_fUS function
oldcodePath = uigetdir('C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Previous lab code\A-US-ULM\SubFunctions\', 'Select the old SubFunctions path');
oldcodePath = [oldcodePath, '\'];
addpath(oldcodePath)

% Load localization processing parameters: proc_params.mat
load([datapath, 'proc_params.mat'])

% [add line to load the recon pixel spacing .mat file]
% For now, hard coding and assuming equal x, y, and z spacing
pix_spacing = P.wl/2;

% Prompt for parameter user input
parameterPrompt = {'Number of files', 'x pixel spacing [um]', 'y pixel spacing [um]', 'z pixel spacing [um]', 'Maximum expected flow speed [mm/s]', 'Persistence frames', 'Moving window size [frames]', 'Acceleration constraint factor', 'Trimmed mean percentage', 'Direction constraint'};
parameterDefaults = {'', num2str(P.Trans.spacingMm * 1e3), num2str(P.Trans.spacingMm * 1e3), num2str(P.wl/2 * 1e6), '50', '3', '3', '2', '20', 'pi/2'};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

xpix_spacing = str2double(parameterUserInput{2});
ypix_spacing = str2double(parameterUserInput{3});
zpix_spacing = str2double(parameterUserInput{4});

%% 2. Turn the individual center files into one cell array and turn the logical matrices into coordinates
numFiles = str2double(parameterUserInput{1});  % # of files (superframes/buffers) to process
totalFrames = numFiles * P.numFramesPerBuffer; % Total frames to process

% Cell array with an entry for each frame. Each entry contains (# bubbles) of coordinate pairs (z, x) of the detected bubble centers
centerCoords = cell(totalFrames, 1); 

% tic
for n = 1:numFiles   % Go through each center file (for each buffer)
    tic
    load([datapath, 'centers-', num2str(n)])
    tsl = size(centers, 1) * size(centers, 2) * size(centers, 3); % troubleshooting length to account for all ones in the centers matrix

    for bfi = 1:size(centers, 4) % buffer frame index
        centersTemp = squeeze(centers(:, :, :, bfi));
        indTemp = find(centersTemp);
        [xc, yc, zc] = ind2sub(size(centersTemp), indTemp);

        if ~((length(xc) == tsl) & (length(yc) == tsl) & (length(zc) == tsl))
            centerCoords{(n - 1) * P.numFramesPerBuffer + bfi} = [xc, yc, zc];
        end
    end
    disp(strcat("Center coordinates for file ", num2str(n), " stored."))
    toc
end
% toc

img_size = size(centers); % Save image size if we want to clear allCenters
% clear bfi cai caiGlobal xc zc centers

%% 3. Calculate bubble count and correct for frames that have bubbles at every voxel
bubbleCount = zeros(length(centerCoords), 1);   % numFiles/# buffers x # frames per buffer. Count of bubbles in each frame
centerCoords_corrected = centerCoords;          % Correct the centerCoords because some frames have every pixel identified as        a bubble

parfor fi = 1:length(centerCoords_corrected) % frame index - go through every frame
    bufTemp = centerCoords{fi};

    % Correct for some error that makes every pixel a bubble
    if size(bufTemp, 1) >= img_size(1) * img_size(2) * img_size(3)
        centerCoords_corrected{fi} = NaN;
        bubbleCount(fi) = 0;
    else
        bubbleCount(fi) = size(bufTemp, 1);
    end
end

totalCount = sum(bubbleCount, 'all');

% Plot the bubble count
% figure; plot(1:length(bubbleCount), bubbleCount, '.')
figure; plot(1:length(bubbleCount), bubbleCount)
title('Bubble count')
xlabel('Frame number')
ylabel('Bubble count')

clear fi bufTemp

%% 4. Define max speed (distance per frame) threshold and initialize variables
startFramePrompt = {'Start frame'}; % User input for the frame to start processing at
startFrameDefault = {'1'};
startFrameUserInput = inputdlg(startFramePrompt, 'Choose start frame #', 1, startFrameDefault);
startFrame = str2double(startFrameUserInput{1});                      % Frame to start processing at
maxSpeedExpectedMMPerS = str2double(parameterUserInput{5});           % max expected flow speed [mm/s]
timePerFrame = 1 / P.frameRate;                                       % time elapsed per frame [s]
maxDistPerFrameM = (maxSpeedExpectedMMPerS / 1000) * timePerFrame;    % max distance traveled per frame [m], according to the max expected flow speed and frame rate
% pixelsPerM = 1 / pix_spacing * imgRefinementFactor(1);              % # of pixels per meter, which depends on the pixel spacing from reconstruction and the image refinement factor from the localization
% maxPixelDistPerFrame = maxDistPerFrameM * pixelsPerM;               % max distance traveled per frame in units of pixels
xpixelsPerM = 1 / (xpix_spacing / 1e6) * imgRefinementFactor(1);      % # of x pixels per meter, which depends on the pixel spacing from reconstruction and the image refinement factor from the localization
ypixelsPerM = 1 / (ypix_spacing / 1e6) * imgRefinementFactor(2);      % # of y pixels per meter, which depends on the pixel spacing from reconstruction and the image refinement factor from the localization
zpixelsPerM = 1 / (zpix_spacing / 1e6) * imgRefinementFactor(3);      % # of z pixels per meter, which depends on the pixel spacing from reconstruction and the image refinement factor from the localization
maxxPixelDistPerFrame = maxDistPerFrameM * xpixelsPerM;               % max x distance traveled per frame in units of pixels
maxyPixelDistPerFrame = maxDistPerFrameM * ypixelsPerM;               % max y distance traveled per frame in units of pixels
maxzPixelDistPerFrame = maxDistPerFrameM * zpixelsPerM;               % max z distance traveled per frame in units of pixels

bubblePairs = cell(totalFrames - 1, 1);   % Initialize cell vector of paired bubble indices

%% 5. Calculate pairing

ubS = cell(totalFrames - 1, 1);             % unassigned bubbles from the source frames
ubT = cell(totalFrames - 1, 1);             % unassigned bubbles from the target frames

tic
parfor f = startFrame:totalFrames - 1       % Go through frames
    sourceFrame = centerCoords_corrected{f};     % Get the coordinates for the source frame (f)
    targetFrame = centerCoords_corrected{f + 1}; % Get the coordinates for the target frame (f + 1)

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
    
            D(D > maxDistPerFrameM) = Inf;      % Set the elements above the distance per frame threshold to Inf so they aren't considered for pairing
        end

        [assignment, unassignedrows, unassignedcolumns] = assignmunkres(D, 100000000000); % Pair with the Munkres algorithm, which minimizes the total cost (total paired distance)
        bubblePairs{f} = assignment; % assignment will always be sorted to make the second column in order
        ubS{f} = unassignedrows;
        ubT{f} = unassignedcolumns;
        
    end
end
toc
disp('Pairing done')

clear f nbS nbT D spi assignment unassignedrows unassignedcolumns

%% 6. Create tracks with persistence
pers = str2double(parameterUserInput{6}); % # of frames a track needs to persist through to keep it

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
mmws = str2double(parameterUserInput{7}); % Moving mean window size [frames]
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
    bVelocityMSmoothedMMS{n}(:, 7, :) = bVelocityMSmoothed{n}(:, 7, :) ./ xpixelsPerM * 1e3;
    bVelocityMSmoothedMMS{n}(:, 8, :) = bVelocityMSmoothed{n}(:, 8, :) ./ ypixelsPerM * 1e3;
    bVelocityMSmoothedMMS{n}(:, 9, :) = bVelocityMSmoothed{n}(:, 9, :) ./ zpixelsPerM * 1e3;
end

disp('Velocity map smoothed')

%% 10. Kalman filter
numDims = 3;
numStates = numDims * 2; % In 3D: x position, y position, z position, x displacement, y displacement, z displacement

if pers < 3
    error('The Kalman filter needs at least 3 points in the track')
end

% Initialize the filtered output variable
bVelocityMSmoothedKFMMS = bVelocityMSmoothedMMS;
vMMStoPixelDispPerFrame = timePerFrame / 1e3 * pixelsPerM;

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
        bVelocityMSmoothedKFMMS{n}(tn, 7:9, :) = xk(4:6, 2:end) ./ vMMStoPixelDispPerFrame;
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
aThresholdFactor = str2double(parameterUserInput{8});         % Acceleration change threshold factor
vTrimmedMeanPercentage = str2double(parameterUserInput{9});   % Trimmed mean percentage for the acceleration change threshold [%]
angleChangeThreshold = str2double(parameterUserInput{10});    % Angle change threshold [radians]

bVelocityMSmoothedKFConstrainedMMS = bVelocityMSmoothedKFMMS; % Initialize the variable for the smoothed, Kalman filtered, constrained velocity map
tic
for n = 1:size(bVelocityMSmoothedKFConstrainedMMS, 1)
    tln = bVelocityMSmoothedKFConstrainedMMS{n};    % Track list n
    for tn = size(tln, 1):-1:1                      % Go through each track number tn
        trackAlreadyDeleted = false;                % Reset the flag
        track = squeeze(tln(tn, :, :))';            % Get the track
        vTrack = track(:, 7:9);                     % Velocities of the track

        % Acceleration constraint
        vTrackTrimmedMean = trimmean(vTrack, vTrimmedMeanPercentage); % Exclude some percentage of the values when taking the mean. It goes across the first non-singleton dimension.
        aThresholdMag = abs(aThresholdFactor .* vTrackTrimmedMean ./ timePerFrame);
        aTrackMag = abs(diff(vTrack, 1) ./ timePerFrame); % Accelerations of the track
        if any(aTrackMag > aThresholdMag, 'all') % Remove the track if the acceleration constraint is violated
            bVelocityMSmoothedKFConstrainedMMS{n}(tn, :, :) = [];
            trackAlreadyDeleted = true;
        end

        % Direction constraint
        if ~trackAlreadyDeleted % Don't need to go through the direction calculation if the track was already deleted for the acceleration constraint
            angleTrack = atan2(vTrack(:, 2), vTrack(:, 1));         % Angle of each segment on the track
            angleTrackChanges = diff(angleTrack);                   % Change in angle between segments on the track
            if any(abs(angleTrackChanges) > angleChangeThreshold)   % Apply the threshold
                bVelocityMSmoothedKFConstrainedMMS{n}(tn, :, :) = [];
            end
        end
    end
end

toc
disp('Acceleration and direction constraints applied')
clear n tln tn trackAlreadyDeleted track vTrack vTrackTrimmedMean aThresholdMag aTrackMag angleTrack angelTrackChanges

%% 12. Plot density map with the paired bubbles after persistence
bSum = zeros(img_size(1), img_size(2), img_size(3)); % Initialize the bubble density map variable
% bSum = padarray(bSum, 50); % Pad array in case the Kalman filter puts some points outside the original region

for n = startFrame:length(bVelocityM) % Go through each track collection and plot
    tempBuf = bVelocityM{n};
%     tempBuf = bVelocityMSmoothedKFConstrainedMMS{n};
%     tempBuf = bVelocityMSmoothedKFMMS{n};
    for tn = 1:size(tempBuf, 1) % Go through each track number tn
        trackTemp = squeeze(tempBuf(tn, :, :))';

        % Add 1 to a pixel's count if the current track intersects it
        bSum(trackTemp(1, 1), trackTemp(1, 2), trackTemp(1, 3)) = bSum(trackTemp(1, 1), trackTemp(1, 2), trackTemp(1, 3)) + 1;
        for iti = 1:size(trackTemp, 1) % inside track index
            bSum(trackTemp(iti, 4), trackTemp(iti, 5), trackTemp(iti, 6)) = bSum(trackTemp(iti, 4), trackTemp(iti, 5), trackTemp(iti, 6)) + 1;
        end
    end
end

clear n tn iti trackTemp tempBuf

volumeViewer(bSum .^ 0.5)

addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Processing\Jerman Enhancement Filter\')

%% Plot speed map after persistence with linear interpolation, on the cleaned and refined velocity data
speedMap = zeros(img_size(1), img_size(2), img_size(3));
plotPower = 1;

% Counters for proper averaging if there are overlapped pixels from
% different tracks
speedMapCounter = zeros(img_size(1), img_size(2), img_size(3));

for ti = startFrame:size(bVelocityMSmoothed, 1)
% for ti = startFrame
%     bvTemp = bVelocityTestMSmoothedKFConstrainedMMS{ti};
    bvTemp = bVelocityMSmoothedKFMMS{ti};
%     bvTemp = bVelocityTestMSmoothedMMS{ti}; % get the ti-th entry
%     bvTemp = bVelocityTestMSmoothed{ti}; % get the ti-th entry
%     bvTemp = bVelocityTestM{ti}; % get the ti-th entry
    if ~isempty(bvTemp) % only do stuff if the bubble velocity cell array entry is not empty
        for bpi = 1:size(bvTemp, 1) % bubble pair index
            % Initialize temporary start and end coordinate matrices.
            % Each have dimensions [# persistence frames, 2] where each row is [z coord, x coord].
            coordsStart = NaN(pers, 3);
            coordsEnd = NaN(pers, 3);

            % Go through the # of persistence frames and get the
            % coordinates at each frame pfi for the bubble track bpi
            for pfi = 1:pers % persistence frame index
                coordsStart(pfi, :) = bvTemp(bpi, 1:3, pfi);
                coordsEnd(pfi, :) = bvTemp(bpi, 4:6, pfi);
            end

            vTemp = squeeze(bvTemp(bpi, 7:9, :)); % Velocity components
            speedTemp = sqrt(vTemp(:, 1).^2 + vTemp(:, 2).^2 + vTemp(:, 3).^2); % Speed

            interpPts = ULM_interp3D_linear(coordsStart, coordsEnd, speedTemp); % Get interpolated points with the corresponding z velocity value. each row is [z coord, x coord, z velocity]
%             interpPts = ULM_interp2D_linear(coordsStart, coordsEnd, zvTemp, ti); % Get interpolated points with the corresponding z velocity value. each row is [z coord, x coord, z velocity]
%             interpPts = ULM_interp2D_spline(coordsStart, coordsEnd, zvTemp, ti);

            for ipi = 1:size(interpPts, 1) % interpolated point index
                interpPtsTemp = interpPts(ipi, :);
                speedValTemp = interpPtsTemp(4);
                
                speedMap(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) = speedMap(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) + speedValTemp;
                speedMapCounter(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) = speedMapCounter(interpPtsTemp(1), interpPtsTemp(2), interpPtsTemp(3)) + 1;
            end
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
clear speedValTemp ti bpi pfi ipi interpPtsTemp speedTemp coordsStart coordsEnd bvTemp

% Take the average for pixels with overlapping tracks
speedMask = speedMapCounter > 0;
speedMap(speedMask) = speedMap(speedMask) ./ speedMapCounter(speedMask);
%%
% Compression test
% zvUpMap = zvUpMap .^ 1/3;
% zvDownMap = -1 .* abs(zvDownMap).^1/3;

% Plot bubble density
bubbleDensityMap = (zvUpMapCounter + zvDownMapCounter);
% figure; imagesc(bubbleDensityMap .^ 0.2); colormap hot

% Make a small PSF, adapting Jianbo's code
FWHM_X=10; % x resolution, FWHM-Amplitude, um
FWHM_Y=10; % y resolution, FWHM-Amplitude, um
FWHM_Z=10;  % z resolution, FWHM-Amplitude, um
Sigma_X=FWHM_X/(2*sqrt(2*log(2)));
Sigma_Y=FWHM_Y/(2*sqrt(2*log(2)));
Sigma_Z=FWHM_Z/(2*sqrt(2*log(2)));
xPSF0=-30:30; yPSF0 = xPSF0; zPSF0=xPSF0; % pixels
[xPSF,yPSF, zPSF]=meshgrid(xPSF0,yPSF0,zPSF0);
% PRSSinfo.sysPSF=exp(-((xPSF/(Sigma_X/PRSSinfo.lPix)).^2+(zPSF/(Sigma_Z/PRSSinfo.lPix)).^2)/2);
smallPSF=exp(-((xPSF/(Sigma_X * (pixelsPerM / 1e6))).^2 + (yPSF/(Sigma_Y * (pixelsPerM / 1e6))).^2 + (zPSF/(Sigma_Z * (pixelsPerM / 1e6))).^2)/2);
% volumeViewer(smallPSF)

bubbleDensityMapConv = conv2(bubbleDensityMap, smallPSF);
% figure; imagesc(bubbleDensityMapConv .^ 0.3); colormap hot

% Velocity map convolution with small PSF
speedMapConv = convn(speedMap, smallPSF);

% Load Jianbo's colormaps
[VzCmap, VzCmapDn, VzCmapUp, pdiCmapUp, PhtmCmap] = Colormaps_fUS;
zvMapFig = figure;

% Plot with two linked axes (one for up Z, other for down Z)
% % hold on
vCrange = [-maxSpeedExpectedMMPerS, maxSpeedExpectedMMPerS];
% vCrange = [-maxSpeedExpectedMMPerS, maxSpeedExpectedMMPerS] .* 0.6;
% vCrange = [-5, 5];
% vCrange = [-8800, 8800];
% vCrange = [-15000, 15000];
figure(zvMapFig)
h1 = axes;
zvUpMap = zvUpMap .^ plotPower;
% imagesc(zvUpMap .^ plotPower)
imagesc(zvUpMap)
% alpha(h1, double(abs(zvUpMap) > 1))
colormap(zvMapFig, VzCmap)
caxis(vCrange);
axis tight
colorbar
hold on

h2 = axes;
zvDownMap = -1 .* (abs(zvDownMap) .^ plotPower);
imagesc(zvDownMap)
alpha(h2, double(abs(zvDownMap) > 1))
colormap(zvMapFig, VzCmap)
caxis(vCrange);
axis tight
cb = colorbar;
axis off
linkaxes([h1, h2]);
ylabel(cb, 'z velocity (mm/s)') % label colorbar

% clim([])
