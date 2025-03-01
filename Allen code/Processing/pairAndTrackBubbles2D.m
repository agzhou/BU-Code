% Description: Part 2 of the ULM processing for a linear array. This script
%   takes the result of the bubble localization from 'LA_analyze_ULM.m' (allCenters) and
%   performs pairing and tracking, plus further image refinement for
%   plotting purposes.

%   allCenters is a cell array with dimensions (1, # of files) where each
%   element contains a logical matrix of dimensions (# refined z pixels, # refined x
%   pixels, # frames per file) with ones at the pixels where a bubble
%   center was localized.

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
%   ULM_interp2D.m (I should update the code to use the below function instead)
%   ULM_interp2D_linear.m
%   Colormaps_fUS (From Jianbo's code)

% Acknowledgement: using Jianbo Tang's ULM code, the Song group's ULM papers, and Jean-Yves Tinevez's simpletracker as references
clearvars
%% 1. Add dependencies and load parameters + bubble center locations
% Get data path of the localized bubble centers
datapath = uigetdir('D:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processed Data\with NLM\', 'Select the data path');
datapath = [datapath, '\'];

% Load localization processing parameters: proc_params.mat
load([datapath, 'proc_params.mat'])

% Choose and load the params file (from the acquisition)
[params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', 'D:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\params.mat');
load([params_pathname, params_filename])

% Add Jianbo's functions to the path for the Colormaps_fUS function
addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Previous lab code\A-US-ULM\SubFunctions\')

% [add line to load the recon pixel spacing .mat file]
% For now, hard coding and assuming equal z and x spacing

% Choose and load the localized centers file
if ~exist('allCenters', 'var') % Only load if the variable doesn't already exist in the workspace
%     allCenters_generalpath = fullfile(datapath, 'allCenters*.mat');
%     allCentersDir = dir(allCenters_generalpath);
%     load([allCentersDir.folder, '\', allCentersDir.name])
    [allCenters_filename, allCenters_pathname, ~] = uigetfile('*.mat', 'Select the allCenters file', 'D:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processed Data\with NLM\allCenters_02_10_2025_SVs10to80_0p4XCThreshold_10refinementfactor.mat');
    load([allCenters_filename, allCenters_pathname])
end

% Prompt for parameter user input
parameterPrompt = {'z pixel spacing [um]', 'x pixel spacing [um]', 'Maximum expected flow speed [mm/s]', 'Persistence frames', 'Moving window size [frames]', 'Acceleration constraint factor', 'Trimmed mean percentage', 'Direction constraint'};
parameterDefaults = {num2str(P.wl/2 * 1e6), num2str(P.Trans.spacingMm * 1e3), '50', '5', '3', '2', '20', 'pi/2'};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

zpix_spacing = str2double(parameterUserInput{1});
xpix_spacing = str2double(parameterUserInput{2});

clear allCenters_generalpath allCentersDir

%% 2. Turn the logical matrices into coordinates
totalFrames = length(allCenters) * size(allCenters{1}, 3); % Get the total number of frames to consider

% Cell array with an entry for each frame. Each entry contains (# bubbles) of coordinate pairs (z, x) of the detected bubble centers
centerCoords = cell(totalFrames, 1); 

% Go through each file/frame and extract the coordinates of the bubble centers
caiGlobal = 1;
for cai = 1:length(allCenters)                          % cell array index - go through each file (buffer/superframe)
    for bfi = 1:size(allCenters{cai}, 3)                % buffer frame index - go through all the frames within each file
        [zc, xc] = find(allCenters{cai}(:, :, bfi));    % Get the coordinates from the logical matrix
        centerCoords{caiGlobal} = [zc, xc];             % Store the coordinates
        caiGlobal = caiGlobal + 1;                      % increment the overall frame index
    end
end

img_size = size(allCenters{1}); % Save image size (helps if we want to clear allCenters)
clear bfi cai caiGlobal xc zc

%% 3. Calculate bubble count and correct for frames that have bubbles at every pixel
bubbleCount = zeros(length(centerCoords), 1);   % numFiles/# buffers x # frames per buffer. Count of bubbles in each frame
centerCoords_corrected = centerCoords;          % Correct the centerCoords because some frames have every pixel identified as a bubble

parfor fi = 1:length(centerCoords_corrected)    % frame index - go through every frame
    bufTemp = centerCoords{fi};                 % Temporary variable - the coordinates of the detected bubbles at frame fi

    % Correct for some error that makes every pixel a bubble
    if size(bufTemp, 1) >= img_size(1) * img_size(2)
        centerCoords_corrected{fi} = NaN;
        bubbleCount(fi) = 0;
    else
        bubbleCount(fi) = size(bufTemp, 1);
    end
end

totalCount = sum(bubbleCount, 'all');           % Total bubble count across the experiment

% Plot the bubble count
% figure; plot(1:length(bubbleCount), bubbleCount, '.')
figure; plot(1:length(bubbleCount), bubbleCount)
title('Bubble count')
xlabel('Frame number')
ylabel('Bubble count')

clear fi bufTemp

%% 4. Define max speed (distance per frame) threshold and initialize variables
maxSpeedExpectedMMPerS = str2double(parameterUserInput{3});           % max expected flow speed [mm/s]
timePerFrame = 1 / P.frameRate;                                       % time elapsed per frame [s]
% totalFrames = size(centerCoords, 1);                                % total number of frames
maxDistPerFrameM = (maxSpeedExpectedMMPerS / 1000) * timePerFrame;    % max distance traveled per frame [m], according to the max expected flow speed and frame rate
zpixelsPerM = 1 / (zpix_spacing / 1e6) * imgRefinementFactor(1);      % # of z pixels per meter, which depends on the pixel spacing from reconstruction and the image refinement factor from the localization
xpixelsPerM = 1 / (xpix_spacing / 1e6) * imgRefinementFactor(2);      % # of x pixels per meter, which depends on the pixel spacing from reconstruction and the image refinement factor from the localization
maxzPixelDistPerFrame = maxDistPerFrameM * zpixelsPerM;               % max z distance traveled per frame in units of pixels
maxxPixelDistPerFrame = maxDistPerFrameM * xpixelsPerM;               % max x distance traveled per frame in units of pixels

bubblePairs = cell(totalFrames - 1, 1);   % Initialize cell vector of paired bubble indices

%% 5. Calculate pairing
% testFig = figure;

ubS = cell(totalFrames - 1, 1);             % unassigned bubbles from the source frames
ubT = cell(totalFrames - 1, 1);             % unassigned bubbles from the target frames

tic
parfor f = 1:totalFrames - 1                % Go through every frame
% for f = 7
    sourceFrame = centerCoords_corrected{f};     % Get the coordinates for the source frame (f)
    targetFrame = centerCoords_corrected{f + 1}; % Get the coordinates for the target frame (f + 1)

    nbS = size(sourceFrame, 1);             % number of bubbles in the source frame
    nbT = size(targetFrame, 1);             % number of bubbles in the target frame
    D = NaN(nbS, nbT);                      % Initialize the distance matrix comparing source and target points' distances

    if nbS > 1 && nbT > 1                   % Only go through the pairing if the source and target frames have at least one bubble
        % Go through each source frame point and get the distance to each
        % point in the target frame. Store in the distance matrix D.
        for spi = 1:nbS % source point index
            sourcePoint = sourceFrame(spi, :);      % z and x coords of the source frame's point "spi"
            d = targetFrame - sourcePoint;          % vectorized difference between the z and x coords of all the points in the target frame and point spi from the source frame
            d(:, 1) = d(:, 1) ./ zpixelsPerM;       % Convert the z distance differences into natural distance units [m]
            d(:, 2) = d(:, 2) ./ xpixelsPerM;       % Convert the x distance differences into natural distance units [m]
            D(spi, :) = sqrt(sum((d .^ 2), 2));     % distance formula on the above
    
            D(D > maxDistPerFrameM) = Inf;          % Set the elements above the distance per frame threshold to Inf so they aren't considered for pairing
        end

        [assignment, unassignedrows, unassignedcolumns] = assignmunkres(D, 100000000000); % Pair with the Munkres algorithm, which minimizes the total cost (total paired distance)
        bubblePairs{f} = assignment; % assignment will always be sorted to make the second column in order
        ubS{f} = unassignedrows;
        ubT{f} = unassignedcolumns;
    
        % Plot a line for each pair
%         figure; scatter(sourceFrame(:, 1), sourceFrame(:, 2), 'ro')
%         hold on; scatter(targetFrame(:, 1), targetFrame(:, 2), 'b*');
%         title(strcat("Frame ", num2str(f)))
%         for npb = 1:size(bubblePairs{f}, 1)
%             pairInd = bubblePairs{f}(npb, :);
%             ZVecTemp = [sourceFrame(pairInd(1), 1), targetFrame(pairInd(2), 1)];
%             XVecTemp = [sourceFrame(pairInd(1), 2), targetFrame(pairInd(2), 2)];
%             line(ZVecTemp, XVecTemp, 'LineWidth', 1.5)
%         end
%         hold off
        
    end
%     hold off
end
toc
disp('Pairing done')

clear f nbS nbT D spi assignment unassignedrows unassignedcolumns

%% 6. Create tracks with the frame persistence condition
pers = str2double(parameterUserInput{4}); % # of frames a track needs to persist through to keep it

tic
% Separate the pairs of coordinates so we can change their sizes independently
bubblePairsPers = cell(length(bubblePairs), 2); % bubble pairs with the pairs separated into another cell dimension
for n = 1:length(bubblePairs)
    if ~isempty(bubblePairs{n})
        bubblePairsPers{n, 1} = bubblePairs{n}(:, 1); % Store the info for the source frame
        bubblePairsPers{n, 2} = bubblePairs{n}(:, 2); % Store the info for the target frame
    end
end

% Store the bubble indices that contribute to each track of at least (pers) # of frames
tracks = cell(size(bubblePairsPers, 1) - pers, 1);
for n = 1:length(bubblePairsPers) - pers
    bubblePairsPersTemp = bubblePairsPers;
    for pfc = 1:pers - 1 % persistence frame count
        startIndex = bubblePairsPersTemp{n + pfc - 1, 2};    % indices of the paired "target" bubbles in frame n + 1 (pair n), which will be sorted in ascending order
        endIndex = bubblePairsPersTemp{n + pfc, 1};  % indices of the paired "source" bubbles in frame n + 2 (pair n + 1), not necessarily in ascending order because it's aligned with the "target" indices in frame n + 1
        [trackContinuesIndices, is, ie] = intersect(startIndex, endIndex, 'stable'); % find the common values in the start and end vectors, and the corresponding indices for each
        
        % Update the paired and tracked list for this immediate set of pairs
        bubblePairsPersTemp{n + pfc - 1, 1} = bubblePairsPersTemp{n + pfc - 1, 1}(is);
        bubblePairsPersTemp{n + pfc - 1, 2} = bubblePairsPersTemp{n + pfc - 1, 2}(is);
        bubblePairsPersTemp{n + pfc, 1} = bubblePairsPersTemp{n + pfc, 1}(ie); % I think this preserves the order
        bubblePairsPersTemp{n + pfc, 2} = bubblePairsPersTemp{n + pfc, 2}(ie); % Need to preserve the order, so set the next starting point????
        

    end
    % Go back and update the previous pairs too
%         for recpfc = pfc - 1 : -1 : 1   % recursive persistence frame count
    for recpfc = 1:pfc - 1   % recursive persistence frame count
        recStartIndex = bubblePairsPersTemp{n + pfc - recpfc - 1, 2};    % indices of the paired "target" bubbles in frame n + 1 (pair n), which will be sorted in ascending order
        recEndIndex = bubblePairsPersTemp{n + pfc - recpfc, 1};  % indices of the paired "source" bubbles in frame n + 2 (pair n + 1), not necessarily in ascending order because it's aligned with the "target" indices in frame n + 1
        [recTrackContinuesIndices, ris, rie] = intersect(recStartIndex, recEndIndex, 'stable'); % find the common values in the start and end vectors, and the corresponding indices for each
    
        bubblePairsPersTemp{n + pfc - recpfc - 1, 1} = bubblePairsPersTemp{n + pfc - recpfc - 1, 1}(ris);
        bubblePairsPersTemp{n + pfc - recpfc - 1, 2} = bubblePairsPersTemp{n + pfc - recpfc - 1, 2}(ris);
    end
    tracks{n} = bubblePairsPersTemp(n : n + pfc, :);
end
toc
clear bubblePairsPersTemp n pfc recpfc

%% Clean tracks

% turn tracks into a proper link of coordinates and indices - remove the
% redundant cross-frame stuff
tracksClean = cell(size(tracks));
nbifAll = zeros(size(tracks));
for n = 1:size(tracks, 1)
% for n = 1:2
    tracksTemp = tracks{n};
    if ~isempty(tracksTemp)
        stt = size(tracksTemp);
        nbif = length(tracksTemp{1}); % # of paired bubbles in each frame of the track
        nbifAll(n) = nbif;
        tracksClean{n} = zeros((stt(1) + 1) * nbif, 4); % There are stt(1) + 1 frames represented in each entry of tracksTemp
        for fn = 1:stt(1)
            tracksClean{n}((fn) * nbif + 1 : (fn + 1) * nbif, 1) = tracksTemp{fn, 2};
            tracksClean{n}((fn) * nbif + 1 : (fn + 1) * nbif, 2:3) = centerCoords_corrected{n + fn}(tracksTemp{fn, 2}, :);
            tracksClean{n}((fn) * nbif + 1 : (fn + 1) * nbif, 4) = repmat(n + fn, nbif, 1); % add frame number
        end
        tracksClean{n}((0) * nbif + 1 : (1) * nbif, 1) = tracksTemp{1, 1};
        tracksClean{n}((0) * nbif + 1 : (1) * nbif, 2:3) = centerCoords_corrected{n}(tracksTemp{1, 1}, :);
        tracksClean{n}((0) * nbif + 1 : (1) * nbif, 4) = repmat(n, nbif, 1); % add frame number
    else
        tracksClean{n} = NaN;
    end

end
clear n tracksTemp nbif stt

%% Velocity map after doing persistence, on the cleaned tracks

bVelocity = cell(size(tracksClean, 1), 1);  % bubble velocity - for each entry in the cell array, [# points x 4] where it has[z position, x position, z velocity, x velocity] in units of pixels and pixels/s
bVelocityTest = cell(size(tracksClean, 1), pers);
bVelocityTestM = cell(size(tracksClean, 1), 1);
for ti = 1:length(tracksClean)              % track index
% for ti = 1:2
        
    tracksTemp = tracksClean{ti};
    nbiti = nbifAll(ti);                      % # of bubbles in the tracks starting in index ti
    vmapTemp = [];
    for fn = 1:pers % Go through all the frames in the tracks with origin frame ti
%     for fn = 1:2
        startPoints = tracksTemp((fn - 1) * nbiti + 1 : fn * nbiti, 2:3);
        endPoints = tracksTemp((fn) * nbiti + 1 : (fn + 1) * nbiti, 2:3);
        vfn = (endPoints - startPoints) ./ timePerFrame;        % velocity = displacement/time
        bVelocityTest{ti, fn} = [startPoints, endPoints, vfn];  % each row is [z start coord, x start coord, z end coord, x end coord, z velocity, x velocity]
        bVelocityTestM{ti}(:, :, fn) = [startPoints, endPoints, vfn];

        for i = 1:nbiti % Go through each pair of coordinates used to calculate the velocity and add the interpolated points to the velocity + interpolated points at which to plot that velocity's matrix
%         for i = 1
            [zcInterp, xcInterp] = ULM_interp2D(startPoints(i, :), endPoints(i, :));
            vmapTempi = [zcInterp, xcInterp, repmat(vfn(i, :), length(zcInterp), 1)];
            vmapTemp = [vmapTemp; vmapTempi];
        end
    end
    bVelocity{ti} = vmapTemp;
end

clear ti fn tracksTemp startPoints endPoints vfn zcInterp xcInterp i vmapTempi vmapTemp

%% Refine the velocity map
mmws = str2double(parameterUserInput{5}); % moving mean window size [frames]
bVelocityTestMSmoothed = bVelocityTestM;
for n = 1:size(bVelocityTestM, 1)
% for n = 1
    vt = bVelocityTestM{n};
    vtSmoothed = vt;
    for tn = 1:size(vt, 1) % track number
%     for tn = 1
        tnzVel = squeeze(vt(tn, 5, :)); % get the z velocity at each point in track tn
        tnzVelSmoothed = movmean(tnzVel, mmws);
        vtSmoothed(tn, 5, :) = tnzVelSmoothed;

        tnxVel = squeeze(vt(tn, 6, :)); % get the x velocity at each point in track tn
        tnxVelSmoothed = movmean(tnxVel, mmws);
        vtSmoothed(tn, 6, :) = tnxVelSmoothed;
        
    end
    bVelocityTestMSmoothed{n} = vtSmoothed;
end
clear n tn vt vtSmoothed tnzVel tnzVelSmoothed vtSmoothed tnxVel tnxVelSmoothed

% scale the smoothed velocity into mm/s
bVelocityTestMSmoothedMMS = bVelocityTestMSmoothed;
for n = 1:size(bVelocityTestMSmoothedMMS, 1)
    bVelocityTestMSmoothedMMS{n}(:, 5, :) = bVelocityTestMSmoothed{n}(:, 5, :) ./ zpixelsPerM * 1e3;
    bVelocityTestMSmoothedMMS{n}(:, 6, :) = bVelocityTestMSmoothed{n}(:, 6, :) ./ xpixelsPerM * 1e3;
end

%% Kalman filter version 1 (observation is only the positions)
numDims = 2;
numStates = numDims * 2; % In 2D: axial position (z), lateral position (x), axial displacement (dz), lateral displacement (dx)

if pers < 3
    error('The Kalman filter needs at least 3 points in the track')
end

% Initialize the filtered output variable
bVelocityTestMSmoothedKFMMS = bVelocityTestMSmoothedMMS;
vMMStoPixelDispPerFrame = timePerFrame / 1e3 * pixelsPerM;

% Define the matrices that map the state transition, and the state ->
% observation transformation
Fk = [1, 0, 1, 0; ...
      0, 1, 0, 1; ...
      0, 0, 1, 0; ...
      0, 0, 0, 1];
Hk = [1, 0, 0, 0; ...
      0, 1, 0, 0];

%%%%%%%% these covariance matrices are from the Song et al. 2020 paper %%%%%%%%%
Qk = diag(ones(numStates, 1));      % Covariance matrix of the system noise
Rk = diag(ones(numDims, 1)) .* 2;   % Covariance matrix of the observation noise

for n = 1:size(bVelocityTestMSmoothedMMS, 1)
    tln = bVelocityTestMSmoothedMMS{n}; % Track list n
    for tn = 1:size(tln, 1) % Track number tn
        track = squeeze(tln(tn, :, :))'; % Get the track
        track = [track; [track(end, 3:4), 0, 0, 0, 0]];
        % Initialize variables for the state vector and covariance matrix
        xk = NaN(numStates, size(track, 1));
        Pk = NaN(numStates, numStates, size(track, 1));

        % Initial values
%         xk(:, 1) = [track(1, 1:2), 0, 0]';    % Initial state vector
        xk(:, 1) = [track(1, 1:2), track(1, 5:6) .* vMMStoPixelDispPerFrame]'; % Initial state vector
        
        Pk(:, :, 1) = [1, 0, 0, 0; ...          % Initial covariance matrix
                 0, 1, 0, 0; ...
                 0, 0, 10, 0; ...
                 0, 0, 0, 10];

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
            yk = [track(k, 1:2)]';

            % Update
            Kku = Pkp * Hk' * inv(Hk * Pkp * Hk' + Rk); % Kalman gain matrix
            Iku = yk - Hk * xkp; % Innovation vector (difference between the observed state and the predicted state transformed into an observation at step k)
            xku = xkp + Kku * Iku; % Updated (weighted) estimate for the state at step k
            Pku = (eye(length(xku)) - Kku * Hk) * Pkp; % Updated covariance matrix at step k

            % Store the updated state and covariance
            xk(:, k) = xku;
            Pk(:, :, k) = Pku;
        end
        bVelocityTestMSmoothedKFMMS{n}(tn, 1:2, :) = xk(1:2, 1:end-1);
        bVelocityTestMSmoothedKFMMS{n}(tn, 3:4, :) = xk(1:2, 2:end);
%         bVelocityTestMSmoothedKFMMS{n}(tn, 5:6, :) = xk(3:4, 2:end);
        bVelocityTestMSmoothedKFMMS{n}(tn, 5:6, :) = (xk(1:2, 2:end) - xk(1:2, 1:end-1)); % ./ timePerFrame;

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

%% Kalman filter version 2 (observation includes velocities)
numDims = 2;
numStates = numDims * 2; % In 2D: axial position (z), lateral position (x), axial displacement (dz), lateral displacement (dx)

if pers < 3
    error('The Kalman filter needs at least 3 points in the track')
end

% Initialize the filtered output variable
bVelocityTestMSmoothedKFMMS = bVelocityTestMSmoothedMMS;
vMMStoPixelDispPerFrame = timePerFrame / 1e3 * pixelsPerM;

% Define the matrices that map the state transition, and the state ->
% observation transformation
Fk = [1, 0, 1, 0; ...
      0, 1, 0, 1; ...
      0, 0, 1, 0; ...
      0, 0, 0, 1];
Hk = [1, 0, 0, 0; ...
      0, 1, 0, 0; ...
      0, 0, 1, 0; ...
      0, 0, 0, 1];

Qk = diag(ones(numStates, 1)) .* 0.5;       % Covariance matrix of the system/process noise
Rk = diag(ones(numStates, 1)) .* 4;         % Covariance matrix of the observation noise

for n = 1:size(bVelocityTestMSmoothedKFMMS, 1)
% for n = 1
    tln = bVelocityTestMSmoothedKFMMS{n};   % Track list n
    for tn = 1:size(tln, 1)                 % Track number tn
%     for tn = 1
%         track = tln(tn, :, :);
        track = squeeze(tln(tn, :, :))';    % Get the track
        track = [track; [track(end, 3:4), 0, 0, 0, 0]];
        % Initialize variables for the state vector and covariance matrix
        xk = NaN(numStates, size(track, 1));
        Pk = NaN(numStates, numStates, size(track, 1));

        % Initial values
%         xk(:, 1) = [track(1, 1:2), 0, 0]'; % Initial state vector
        %%%%%% WHAT INITIAL DISPLACEMENT/VELOCITY SHOULD I USE??? %%%%%%
%         Pk{1} = eye(length(x0));                 % Initial covariance matrix
        xk(:, 1) = [track(1, 1:2), track(1, 5:6) .* vMMStoPixelDispPerFrame]'; % Initial state vector

        Pk(:, :, 1) = [1, 0, 0, 0; ...
                 0, 1, 0, 0; ...
                 0, 0, 10, 0; ...
                 0, 0, 0, 10];
        %%%%%% WHAT INITIAL COVARIANCES SHOULD I USE??? %%%%%%

%         xk(:, 2) = [track(2, 1:2), track(1, 5:6) .* vMMStoPixelDispPerFrame]';
% %         xk(:, 2) = [track(2, 1:2), track(1, 5:6)]';
%         Pk(:, :, 2) = [1, 0, 0, 0; ...
%                  0, 1, 0, 0; ...
%                  0, 0, 10^4, 0; ...
%                  0, 0, 0, 10^4];

        for k = 2:size(track, 1) % Go through each step of the track
%             xk{k} = Fk * xk{k - 1} + Qk;
            % Prediction
            xkp = Fk * xk(:, k - 1); 
            Pkp = Fk * Pk(:, :, k - 1) * Fk + Qk;

            % Observation
            yk = [track(k, 1:2), track(k, 5:6) .* vMMStoPixelDispPerFrame]';

            % Update
            Kku = Pkp * Hk' * inv(Hk * Pkp * Hk' + Rk); % Kalman gain matrix
            Iku = yk - Hk * xkp; % Innovation vector (difference between the observed state and the predicted state transformed into an observation at step k)
            xku = xkp + Kku * Iku; % Updated (weighted) estimate for the state at step k
            Pku = (eye(length(xku)) - Kku * Hk) * Pkp; % Updated covariance matrix at step k

            % Store the updated state and covariance
            xk(:, k) = xku;
            Pk(:, :, k) = Pku;
        end
        bVelocityTestMSmoothedKFMMS{n}(tn, 1:2, :) = xk(1:2, 1:end-1);
        bVelocityTestMSmoothedKFMMS{n}(tn, 3:4, :) = xk(1:2, 2:end);
        bVelocityTestMSmoothedKFMMS{n}(tn, 5:6, :) = xk(3:4, 2:end) ./ vMMStoPixelDispPerFrame;
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

clear n tn tln k xk Pk yk Kku Iku xku Pku track
%% Acceleration and direction constraints
% aThresholdFactor = 0.5;
aThresholdFactor = str2double(parameterUserInput{6});        % Acceleration change threshold factor
vTrimmedMeanPercentage = str2double(parameterUserInput{7});  % Trimmed mean percentage for the acceleration change threshold [%]
angleChangeThreshold = str2double(parameterUserInput{8});    % Angle change threshold [radians]

bVelocityTestMSmoothedKFConstrainedMMS = bVelocityTestMSmoothedKFMMS;
for n = 1:size(bVelocityTestMSmoothedKFConstrainedMMS, 1)
    tln = bVelocityTestMSmoothedKFConstrainedMMS{n}; % Track list n
    for tn = size(tln, 1):-1:1              % Go through each track number tn
        trackAlreadyDeleted = false;        % Reset the flag
        track = squeeze(tln(tn, :, :))';    % Get the track
        vTrack = track(:, 5:6);             % Velocities of the track

        % Acceleration constraint
        vTrackTrimmedMean = trimmean(vTrack, vTrimmedMeanPercentage); % Exclude some percentage of the values when taking the mean. It goes across the first non-singleton dimension.
        aThresholdMag = abs(aThresholdFactor .* vTrackTrimmedMean ./ timePerFrame);
        aTrackMag = abs(diff(vTrack, 1) ./ timePerFrame); % Accelerations of the track
        if any(aTrackMag > aThresholdMag, 'all') % Remove the track if the acceleration constraint is violated
            bVelocityTestMSmoothedKFConstrainedMMS{n}(tn, :, :) = [];
            trackAlreadyDeleted = true;
        end

        % Direction constraint
        if ~trackAlreadyDeleted % Don't need to go through the direction calculation if the track was already deleted for the acceleration constraint
            angleTrack = atan2(vTrack(:, 2), vTrack(:, 1));         % Angle of each segment on the track
            angleTrackChanges = diff(angleTrack);                   % Change in angle between segments on the track
            if any(abs(angleTrackChanges) > angleChangeThreshold)   % Apply the threshold
                bVelocityTestMSmoothedKFConstrainedMMS{n}(tn, :, :) = [];
            end
        end
    end
end

%% Plot density map with the paired bubbles only
% % bSum = zeros(size(allCenters{1}, 1), size(allCenters{1}, 2));
% bSum = zeros(img_size(1), img_size(2));
% 
% bSumFig = figure; colormap turbo
% % hold on
% figure(bSumFig)
% 
% % bSumFig.XDataSource = 
% 
% for f = 1:length(bubblePairs)
% % for f = 1:10000
%     disp(num2str(f))
%     if ~isempty(bubblePairs{f})
%         coordsCF = centerCoords_corrected{f}; % coords for the current frame
%         coordsNF = centerCoords_corrected{f + 1}; % coords for the current frame
%     %     for npb = 1:size(pairedBubbles{f}, 1)
% 
%         keepBubblesCF = bubblePairs{f}(:, 1); % which bubbles (indices) to keep for the current frame. Could do it more efficiently or correctly
%         for kbcfi = keepBubblesCF' % index for each kept bubble in the current frame
%             kbcfCoord = coordsCF(kbcfi, :);
%             bSum(kbcfCoord(1), kbcfCoord(2)) = bSum(kbcfCoord(1), kbcfCoord(2)) + 1; % Update the count for that pixel
%         end
% 
%         keepBubblesNF = bubblePairs{f}(:, 2); % which bubbles (indices) to keep for the next frame. Could do it more efficiently or correctly
%         for kbnfi = keepBubblesNF' % index for each kept bubble in the next frame
%             kbnfCoord = coordsNF(kbnfi, :);
%             bSum(kbnfCoord(1), kbnfCoord(2)) = bSum(kbnfCoord(1), kbnfCoord(2)) + 1; % Update the count for that pixel
%         end
% 
%     end
% end
% hold off
% bsIm = imagesc(bSum);

%% Plot z velocity map after persistence with linear interpolation, on the cleaned and refined velocity data
% zVelocityPersMap = zeros(img_size(1), img_size(2));
zvUpMap = zeros(img_size(1), img_size(2));
zvDownMap = zeros(img_size(1), img_size(2));
plotPower = 1;

% Counters for proper averaging if there are overlapped pixels from
% different tracks
zvUpMapCounter = zeros(img_size(1), img_size(2));
zvDownMapCounter = zeros(img_size(1), img_size(2));

for ti = 1:size(bVelocityTestMSmoothed, 1)
% for ti = 7
%     bvTemp = bVelocityTestMSmoothedKFConstrainedMMS{ti};
%     bvTemp = bVelocityTestMSmoothedKFMMS{ti};
    bvTemp = bVelocityTestMSmoothedMMS{ti}; % get the ti-th entry
%     bvTemp = bVelocityTestMSmoothed{ti}; % get the ti-th entry
%     bvTemp = bVelocityTestM{ti}; % get the ti-th entry
    if ~isempty(bvTemp) % only do stuff if the bubble velocity cell array entry is not empty
        for bpi = 1:size(bvTemp, 1) % bubble pair index
%             clear interpPts coordsStart coordsEnd zvTemp
%             interpPts = [];
%         for bpi = 2
            % Initialize temporary start and end coordinate matrices.
            % Each have dimensions [# persistence frames, 2] where each row is [z coord, x coord].
            coordsStart = NaN(pers, 2);
            coordsEnd = NaN(pers, 2);

            % Go through the # of persistence frames and get the
            % coordinates at each frame pfi for the bubble track bpi
            for pfi = 1:pers % persistence frame index
                coordsStart(pfi, :) = bvTemp(bpi, 1:2, pfi);
                coordsEnd(pfi, :) = bvTemp(bpi, 3:4, pfi);
            end

            zvTemp = squeeze(bvTemp(bpi, 5, :));
            interpPts = ULM_interp2D_linear(coordsStart, coordsEnd, zvTemp); % Get interpolated points with the corresponding z velocity value. each row is [z coord, x coord, z velocity]
%             interpPts = ULM_interp2D_linear(coordsStart, coordsEnd, zvTemp, ti); % Get interpolated points with the corresponding z velocity value. each row is [z coord, x coord, z velocity]
%             interpPts = ULM_interp2D_spline(coordsStart, coordsEnd, zvTemp, ti);

            for ipi = 1:size(interpPts, 1) % interpolated point index
                interpPtsTemp = interpPts(ipi, :);
                zVelTemp = interpPtsTemp(3);
                if zVelTemp > 0         % up z flow
                    zvUpMap(interpPtsTemp(1), interpPtsTemp(2)) = zvUpMap(interpPtsTemp(1), interpPtsTemp(2)) + zVelTemp;
                    zvUpMapCounter(interpPtsTemp(1), interpPtsTemp(2)) = zvUpMapCounter(interpPtsTemp(1), interpPtsTemp(2)) + 1;
                else
                    zvDownMap(interpPtsTemp(1), interpPtsTemp(2)) = zvDownMap(interpPtsTemp(1), interpPtsTemp(2)) + zVelTemp;
                    zvDownMapCounter(interpPtsTemp(1), interpPtsTemp(2)) = zvDownMapCounter(interpPtsTemp(1), interpPtsTemp(2)) + 1;
                end
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
clear zVelTemp zvUpMask zvDownMask ti bpi pfi ipi interpPtsTemp zvTemp coordsStart coordsEnd bvTemp

% Take the average for pixels with overlapping tracks
zvUpMask = zvUpMapCounter > 0;
zvDownMask = zvDownMapCounter > 0;
zvUpMap(zvUpMask) = zvUpMap(zvUpMask) ./ zvUpMapCounter(zvUpMask);
zvDownMap(zvDownMask) = zvDownMap(zvDownMask) ./ zvDownMapCounter(zvDownMask);

% Compression test
% zvUpMap = zvUpMap .^ 1/3;
% zvDownMap = -1 .* abs(zvDownMap).^1/3;

% Plot bubble density
bubbleDensityMap = (zvUpMapCounter + zvDownMapCounter);
% figure; imagesc(bubbleDensityMap .^ 0.2); colormap hot

% Make a small PSF, adapting Jianbo's code
FWHM_X=10; % lateral resolution, FWHM-Amplitude, um
FWHM_Z=10;  % axial resolution, FWHM-Amplitude, um
Sigma_X=FWHM_X/(2*sqrt(2*log(2)));
Sigma_Z=FWHM_Z/(2*sqrt(2*log(2)));
xPSF0=-30:30; zPSF0=xPSF0; % pixels
[xPSF,zPSF]=meshgrid(xPSF0,zPSF0);
% PRSSinfo.sysPSF=exp(-((xPSF/(Sigma_X/PRSSinfo.lPix)).^2+(zPSF/(Sigma_Z/PRSSinfo.lPix)).^2)/2);
smallPSF=exp(-((xPSF/(Sigma_X * (xpixelsPerM / 1e6))).^2+(zPSF/(Sigma_Z * (zpixelsPerM / 1e6))).^2)/2);
% figure; imagesc(smallPSF)

bubbleDensityMapConv = conv2(bubbleDensityMap, smallPSF);
% figure; imagesc(bubbleDensityMapConv .^ 0.3); colormap hot

% Velocity map convolution with small PSF
zvUpMap = conv2(zvUpMap, smallPSF);
zvDownMap = conv2(zvDownMap, smallPSF);

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
