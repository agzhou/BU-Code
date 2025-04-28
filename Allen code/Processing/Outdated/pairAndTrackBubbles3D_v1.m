% Acknowledgement: using Jean-Yves Tinevez's simpletracker as a reference
%% Add dependencies and load parameters + bubble center locations
% datapath = 'D:\Allen\Data\01-29-2025 AZ001 ULM RC15gV\run 1 left eye\Processed Data 02-21-2025\';
% datapath = 'F:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\Processed Data 02-24-2025\';
datapath = 'F:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\Processed Data 02-27-2025 in the hole\';
% addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Previous lab code\A-US-ULM\SubFunctions\')
addpath('C:\Users\BOAS-US\Documents\Allen\GitHub\BU-Code\Previous lab code\A-US-ULM\SubFunctions')
% load('D:\Allen\Data\01-29-2025 AZ001 ULM RC15gV\run 1 left eye\params.mat')
load('G:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\params.mat')

% Load localization processing parameters: proc_params.mat
load([datapath, 'proc_params.mat'])

% [add line to load the recon pixel spacing .mat file]
% For now, hard coding and assuming equal z and x spacing
pix_spacing = P.wl/2;

% Find and load the localized centers file, which starts with 'allCenters'
% if ~exist('allCenters', 'var') % Only load if the variable doesn't already exist in the workspace
%     allCenters_generalpath = fullfile(datapath, 'allCenters*.mat');
%     allCentersDir = dir(allCenters_generalpath);
%     load([allCentersDir.folder, '\', allCentersDir.name])
% end
% clear allCenters_generalpath allCentersDir
%% Turn the individual center files into one cell array and turn the logical matrices into coordinates
numFiles = 96;
totalFrames = numFiles * P.numFramesPerBuffer;
% allCenters = cell(numFiles, 1);

% Cell array with an entry for each frame. Each entry contains (# bubbles) of coordinate pairs (z, x) of the detected bubble centers
centerCoords = cell(totalFrames, 1); 

% tic
% caiGlobal = 1;
for n = 1:numFiles   % Go through each center file (for each buffer)
% for n = 30
    tic
    load([datapath, 'centers-', num2str(n)])
    tsl = size(centers, 1) * size(centers, 2) * size(centers, 3); % troubleshooting length to account for all ones in the centers matrix

    for bfi = 1:size(centers, 4) % buffer frame index
%     for bfi = 1
        centersTemp = squeeze(centers(:, :, :, bfi));
        indTemp = find(centersTemp);
        [xc, yc, zc] = ind2sub(size(centersTemp), indTemp);

        if ~((length(xc) == tsl) & (length(yc) == tsl) & (length(zc) == tsl))
            centerCoords{(n - 1) * P.numFramesPerBuffer + bfi} = [xc, yc, zc];
        end
%         caiGlobal = caiGlobal + 1;
    end
    disp(strcat("Center coords for file ", num2str(n), " stored."))
    toc
end
% toc

img_size = size(centers); % Save image size if we want to clear allCenters
% clear bfi cai caiGlobal xc zc centers

%% Calculate bubble count
bubbleCount = zeros(length(centerCoords), 1);   % numFiles/# buffers x # frames per buffer. Count of bubbles in each frame
centerCoords_corrected = centerCoords;          % Correct the centerCoords because some frames have every pixel identified as        a bubble
% mbci = 1; % microbubble count index
parfor fi = 1:length(centerCoords_corrected) % frame index
% for fi = 1
    bufTemp = centerCoords{fi};

    % Some error that makes every pixel a bubble
    if size(bufTemp, 1) >= img_size(1) * img_size(2) * img_size(3)
        centerCoords_corrected{fi} = NaN;
        bubbleCount(fi) = 0;
    else
        bubbleCount(fi) = size(bufTemp, 1);
    end
%     bubbleCount(mbci) = sum(bufTemp(:, :, f), 'all');
%     mbci = mbci + 1;
end

totalCount = sum(bubbleCount, 'all');

% Plot the bubble count
% figure; plot(1:length(bubbleCount), bubbleCount, '.')
figure; plot(1:length(bubbleCount), bubbleCount)
title('Bubble count')
xlabel('Frame number')
ylabel('Bubble count')

clear fi bufTemp

%% Max speed (distance per frame) threshold and initialize variables
startFrame = 4001;                                                  % Frame to start processing at
maxSpeedExpectedMMPerS = 50;                                       % max expected flow speed [mm/s]
timePerFrame = 1 / P.frameRate;                                     % time elapsed per frame [s]
totalFrames = size(centerCoords, 1);
maxDistPerFrameM = (maxSpeedExpectedMMPerS / 1000) * timePerFrame;  % max distance traveled per frame [m], according to the max expected flow speed and frame rate
pixelsPerM = 1 / pix_spacing * imgRefinementFactor(1);              % # of pixels per meter, which depends on the pixel spacing from reconstruction and the image refinement factor from the localization
maxPixelDistPerFrame = maxDistPerFrameM * pixelsPerM;               % max distance traveled per frame in units of pixels

bubblePairs = cell(totalFrames - 1, 1);   % Initialize cell vector of paired bubble indices

%% Calculate pairing
% testFig = figure;

ubS = cell(totalFrames - 1, 1);             % unassigned bubbles from the source frames
ubT = cell(totalFrames - 1, 1);             % unassigned bubbles from the target frames

tic
parfor f = startFrame:totalFrames - 1
% for f = 7
    sourceFrame = centerCoords_corrected{f};
    targetFrame = centerCoords_corrected{f + 1};

    nbS = size(sourceFrame, 1);             % number of bubbles in the source frame
    nbT = size(targetFrame, 1);             % number of bubbles in the target frame
    D = NaN(nbS, nbT);

    if nbS > 1 && nbT > 1
        % Go through each source frame point and get the distance to each
        % point in the target frame. Store in the distance matrix D.
        for spi = 1:nbS % source point index
    %     for cpi = 1
            sourcePoint = sourceFrame(spi, :);      % z and x coords of the source frame's point "spi"
            d = targetFrame - sourcePoint;          % vectorized difference between the z and x coords of all the points in the target frame and point spi from the source frame
            D(spi, :) = sqrt(sum((d .^ 2), 2));     % distance formula on the above
    
            D(D > maxPixelDistPerFrame) = Inf; 
        end

        [assignment, unassignedrows, unassignedcolumns] = assignmunkres(D, 100000000000);
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

%% Create tracks with persistence
pers = 2; % # of frames a track needs to persist through to keep it

tic
% Separate the pairs of coordinates so we can change their sizes independently
bubblePairsPers = cell(length(bubblePairs), 2); % bubble pairs with the pairs separated into another cell dimension
for n = startFrame:length(bubblePairs)
    if ~isempty(bubblePairs{n})
        bubblePairsPers{n, 1} = bubblePairs{n}(:, 1);
        bubblePairsPers{n, 2} = bubblePairs{n}(:, 2);
    end
end

% Store the bubble indices that contribute to each track of at least (pers) # of frames
tracks = cell(size(bubblePairsPers, 1) - pers, 1);
for n = 1:length(bubblePairsPers) - pers
% for n = 7
    bubblePairsPersTemp = bubblePairsPers;
    for pfc = 1:pers - 1 % persistence frame count
%     for pfc = 1
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
disp('Tracks with persistence created')
toc
clear bubblePairsPersTemp n pfc recpfc ris rie recTrackContinuesIndices recStartIndex recEndIndex trackContinuesIndices is ie startIndex endIndex

%% Clean tracks

% turn tracks into a proper link of coordinates and indices - remove the
% redundant cross-frame stuff
tracksClean = cell(size(tracks));
nbifAll = zeros(size(tracks)); % # bubbles in frame, all
for n = startFrame:size(tracks, 1)
% for n = 1:2
    tracksTemp = tracks{n};
    if ~isempty(tracksTemp)
        stt = size(tracksTemp);
        nbif = length(tracksTemp{1}); % # of paired bubbles in each frame of the track
        nbifAll(n) = nbif;
        tracksClean{n} = zeros((stt(1) + 1) * nbif, 5); % There are stt(1) + 1 frames represented in each entry of tracksTemp
        for fn = 1:stt(1)
            tracksClean{n}((fn) * nbif + 1 : (fn + 1) * nbif, 1) = tracksTemp{fn, 2};
            tracksClean{n}((fn) * nbif + 1 : (fn + 1) * nbif, 2:4) = centerCoords_corrected{n + fn}(tracksTemp{fn, 2}, :);
            tracksClean{n}((fn) * nbif + 1 : (fn + 1) * nbif, 5) = repmat(n + fn, nbif, 1); % add frame number
        end
        tracksClean{n}((0) * nbif + 1 : (1) * nbif, 1) = tracksTemp{1, 1};
        tracksClean{n}((0) * nbif + 1 : (1) * nbif, 2:4) = centerCoords_corrected{n}(tracksTemp{1, 1}, :);
        tracksClean{n}((0) * nbif + 1 : (1) * nbif, 5) = repmat(n, nbif, 1); % add frame number
    else
        tracksClean{n} = NaN;
    end

end
clear n tracksTemp nbif stt

%% Velocity map after doing persistence, on the cleaned tracks

% bVelocity = cell(size(tracksClean, 1), 1);  % bubble velocity - for each entry in the cell array, [# points x 4] where it has[z position, x position, z velocity, x velocity] in units of pixels and pixels/s
bVelocityTest = cell(size(tracksClean, 1), pers);
bVelocityTestM = cell(size(tracksClean, 1), 1);
tic
for ti = startFrame:length(tracksClean)              % track index
% for ti = 1:2
        
    tracksTemp = tracksClean{ti};
    nbiti = nbifAll(ti);                      % # of bubbles in the tracks starting in index ti
%     vmapTemp = [];
    if nbiti > 0
        for fn = 1:pers % Go through all the frames in the tracks with origin frame ti
    %     for fn = 1:2
            startPoints = tracksTemp((fn - 1) * nbiti + 1 : fn * nbiti, 2:4);
            endPoints = tracksTemp((fn) * nbiti + 1 : (fn + 1) * nbiti, 2:4);
            vfn = (endPoints - startPoints) ./ timePerFrame;        % velocity = displacement/time
            bVelocityTest{ti, fn} = [startPoints, endPoints, vfn];  % each row is [x start coord, y start coord, z start coord, x end coord, y end coord, z end coord, x velocity, y velocity, z velocity]
            bVelocityTestM{ti}(:, :, fn) = [startPoints, endPoints, vfn];
    
    %         for i = 1:nbiti % Go through each pair of coordinates used to calculate the velocity and add the interpolated points to the velocity + interpolated points at which to plot that velocity's matrix
    % %         for i = 1
    %             [zcInterp, xcInterp] = ULM_interp2D(startPoints(i, :), endPoints(i, :));
    %             vmapTempi = [zcInterp, xcInterp, repmat(vfn(i, :), length(zcInterp), 1)];
    %             vmapTemp = [vmapTemp; vmapTempi];
    %         end
        end
    end
%     bVelocity{ti} = vmapTemp;
end
toc
disp('Velocity map created')

clear ti fn tracksTemp startPoints endPoints vfn zcInterp xcInterp i vmapTempi vmapTemp nbiti

%% Combine cleaned tracks
tracksCombined = cell(size(tracksClean));
% for n = 1:size(tracksClean, 1) % go through frames
for n = 1
%     tracksntemp = tracksClean{n};            % tracks starting in frame n
%     nbifn = nbifAll(n);                      % # of bubbles in the tracks starting in frame n
%     tcTemp = zeros(1, 4);
    for fn = n:n + pers           % go through all the possible overlaps in frames
%     for fn = 1:20
        tcTemp = [];
%         tcTemp = tracksntemp((fn - 1) * nbifn + 1:fn * nbifn, :); % Accumulate the possible bubbles/coords for each global frame # fn
        rangefn = n:fn;  % range of tracks which include frame # fn

        for ti = rangefn % track index
            trackstiTemp = tracksClean{ti};  % tracks in index ti
            nbifti = nbifAll(ti);            % # of bubbles in the tracks starting in frame ti

            tracksfnTemp = trackstiTemp((fn - ti) * nbifti + 1:(fn - ti + 1) * nbifti, :); % tracks of global frame # fn within tracks{ti}
            tcTemp = [tcTemp; tracksfnTemp];
        end

        [uniqueTracksfnTemp, iat, ict] = unique(tcTemp, 'rows', 'stable');
        tracksCombined{fn} = uniqueTracksfnTemp; %%%%

%         for nfi = 1 % frame number inside tracknptemp
%             tracksntempfni = tracksntemp((nfi) * nbifn + 1 : (nfi + 1) * nbifn, :);
%             tracksnptempfni = tracksnptemp((nfi - 1) * nbifnp + 1 : nfi * nbifnp, :);
% %             tcTemp = [tcTemp; unique([tracksntempfni; tracksnptempfni], 'rows', 'stable')];
%             tcTempnp = [tcTempnp; unique([tracksntempfni; tracksnptempfni], 'rows', 'stable')];
%         end
    end
end
% clear trackntemp tracknptemp tcTemp n np

%% test on tracks
figure;
hold on
% for ti = 1:length(tracks)% track index
for ti = 1
    tracksTemp = tracks{ti};
    for tli = 1:size(tracksTemp, 1)% track length index
        startPointInd = tracksTemp{tli, 1};
        endPointInd = tracksTemp{tli, 2};

        coordsSFTemp = centerCoords_corrected{ti + tli - 1}(startPointInd, :); % coords for the source frame
        coordsTFTemp = centerCoords_corrected{ti + tli}(endPointInd, :); % coords for the target frame
        
        scatter(coordsSFTemp(:, 1), coordsSFTemp(:, 2), 'ro')
        scatter(coordsTFTemp(:, 1), coordsTFTemp(:, 2), 'bx')
%         for idk = size(startPointInd, 1)
%             line([coordsSFTemp(idk, 1); coordsTFTemp(idk, 1)], [coordsSFTemp(idk, 2); coordsTFTemp(idk, 2)])
%         end

    end
end
% clear startPointInd endPointInd coordsSFTemp coordsTFTemp

%% Refine the velocity map - (maybe remove tracks where the position is oscillating over time)
mmws = 3; % moving mean window size
bVelocityTestMSmoothed = bVelocityTestM;
for n = startFrame:size(bVelocityTestM, 1)
% for n = 1
    vt = bVelocityTestM{n};
    vtSmoothed = vt;
    for tn = 1:size(vt, 1) % track number
%     for tn = 1
        tnxVel = squeeze(vt(tn, 7, :)); % get the x velocity at each point in track tn
        tnxVelSmoothed = movmean(tnxVel, mmws);
        vtSmoothed(tn, 7, :) = tnxVelSmoothed;

        tnyVel = squeeze(vt(tn, 8, :)); % get the y velocity at each point in track tn
        tnyVelSmoothed = movmean(tnyVel, mmws);
        vtSmoothed(tn, 8, :) = tnyVelSmoothed;

        tnzVel = squeeze(vt(tn, 9, :)); % get the z velocity at each point in track tn
        tnzVelSmoothed = movmean(tnzVel, mmws);
        vtSmoothed(tn, 9, :) = tnzVelSmoothed;
        
    end
    bVelocityTestMSmoothed{n} = vtSmoothed;
end
clear n tn vt vtSmoothed tnzVel tnzVelSmoothed vtSmoothed tnxVel tnxVelSmoothed tnyVel tnyVelSmoothed 

% scale the smoothed velocity into mm/s
bVelocityTestMSmoothedMMS = bVelocityTestMSmoothed;
for n = startFrame:size(bVelocityTestMSmoothedMMS, 1)
    bVelocityTestMSmoothedMMS{n}(:, 7:9, :) = bVelocityTestMSmoothed{n}(:, 7:9, :) ./ pixelsPerM * 1e3;
end

%% Kalman filter attempt 1
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

% for n = 1:size(bVelocityTestMSmoothedMMS, 1)
for n = 1
    tln = bVelocityTestMSmoothedMMS{n}; % Track list n
    for tn = 1:size(tln, 1) % Track number tn
%     for tn = 1
%         track = tln(tn, :, :);
        track = squeeze(tln(tn, :, :))'; % Get the track
        track = [track; [track(end, 3:4), 0, 0, 0, 0]];
        % Initialize variables for the state vector and covariance matrix
        xk = NaN(numStates, size(track, 1));
        Pk = NaN(numStates, numStates, size(track, 1));

        % Initial values
        xk(:, 1) = [track(1, 1:2), 0, 0]'; % Initial state vector
%         xk(:, 1) = [track(1, 1:2), track(1, 5:6) .* vMMStoPixelDispPerFrame]'; % Initial state vector
        %%%%%% WHAT INITIAL DISPLACEMENT/VELOCITY SHOULD I USE??? %%%%%%
%         Pk{1} = eye(length(x0));                 % Initial covariance matrix
        Pk(:, :, 1) = [1, 0, 0, 0; ...
                 0, 1, 0, 0; ...
                 0, 0, 0, 0; ...
                 0, 0, 0, 0];
        %%%%%% WHAT INITIAL COVARIANCES SHOULD I USE??? %%%%%%

        xk(:, 2) = [track(2, 1:2), track(1, 5:6) .* vMMStoPixelDispPerFrame]';
%         xk(:, 2) = [track(2, 1:2), track(1, 5:6)]';
        Pk(:, :, 2) = [1, 0, 0, 0; ...
                 0, 1, 0, 0; ...
                 0, 0, 10^4, 0; ...
                 0, 0, 0, 10^4];

        for k = 3:size(track, 1) % Go through each step of the track
%             xk{k} = Fk * xk{k - 1} + Qk;
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
        figure
        hold on
        plot(xk(1, :), xk(2, :), '-o')
        plot(track(:, 1), track(:, 2), '--x')
        clear i xki
        legend('Kalman filtered', 'Original')
        hold off
    end
end

%% Kalman filter attempt 2
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
% for n = 1
    tln = bVelocityTestMSmoothedMMS{n}; % Track list n
    for tn = 1:size(tln, 1) % Track number tn
%     for tn = 1
%         track = tln(tn, :, :);
        track = squeeze(tln(tn, :, :))'; % Get the track
        track = [track; [track(end, 3:4), 0, 0, 0, 0]];
        % Initialize variables for the state vector and covariance matrix
        xk = NaN(numStates, size(track, 1));
        Pk = NaN(numStates, numStates, size(track, 1));

        % Initial values
%         xk(:, 1) = [track(1, 1:2), 0, 0]'; % Initial state vector
        xk(:, 1) = [track(1, 1:2), track(1, 5:6) .* vMMStoPixelDispPerFrame]'; % Initial state vector
        %%%%%% WHAT INITIAL DISPLACEMENT/VELOCITY SHOULD I USE??? %%%%%%
%         Pk{1} = eye(length(x0));                 % Initial covariance matrix
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

%% Kalman filter attempt 3
numDims = 3;
numStates = numDims * 2; % In 3D: x position, y position, z position, x displacement, y displacement, z displacement

if pers < 3
    error('The Kalman filter needs at least 3 points in the track')
end

% Initialize the filtered output variable
bVelocityTestMSmoothedKFMMS = bVelocityTestMSmoothedMMS;
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

for n = startFrame:size(bVelocityTestMSmoothedKFMMS, 1)
% for n = startFrame
    tln = bVelocityTestMSmoothedKFMMS{n}; % Track list n
    for tn = 1:size(tln, 1) % Track number tn
%     for tn = 1
%         track = tln(tn, :, :);
        track = squeeze(tln(tn, :, :))'; % Get the track
        track = [track; [track(end, 4:6), 0, 0, 0, 0, 0, 0]];
        % Initialize variables for the state vector and covariance matrix
        xk = NaN(numStates, size(track, 1));
        Pk = NaN(numStates, numStates, size(track, 1));

        % Initial values
%         xk(:, 1) = [track(1, 1:2), 0, 0]'; % Initial state vector
        %%%%%% WHAT INITIAL DISPLACEMENT/VELOCITY SHOULD I USE??? %%%%%%
%         Pk{1} = eye(length(x0));                 % Initial covariance matrix
        xk(:, 1) = [track(1, 1:3), track(1, 7:9) .* vMMStoPixelDispPerFrame]'; % Initial state vector

        Pk(:, :, 1) = [1, 0, 0, 0, 0, 0; ...
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
%             xk{k} = Fk * xk{k - 1} + Qk;
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
        bVelocityTestMSmoothedKFMMS{n}(tn, 1:3, :) = round(xk(1:3, 1:end-1));
        bVelocityTestMSmoothedKFMMS{n}(tn, 4:6, :) = round(xk(1:3, 2:end));
        bVelocityTestMSmoothedKFMMS{n}(tn, 7:9, :) = xk(4:6, 2:end) ./ vMMStoPixelDispPerFrame;
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
aThresholdFactor = 2;
vTrimmedMeanPercentage = 20;
angleChangeThreshold = pi/2; % Angle change threshold [radians]

bVelocityTestMSmoothedKFConstrainedMMS = bVelocityTestMSmoothedKFMMS;
for n = 1:size(bVelocityTestMSmoothedKFConstrainedMMS, 1)
% for n = 1
    tln = bVelocityTestMSmoothedKFConstrainedMMS{n}; % Track list n
    for tn = size(tln, 1):-1:1 % Track number tn
%     for tn = 1
        trackAlreadyDeleted = false;
        track = squeeze(tln(tn, :, :))'; % Get the track
        vTrack = track(:, 7:9); % Velocities of the track

        % Acceleration constraint
        vTrackTrimmedMean = trimmean(vTrack, vTrimmedMeanPercentage); % Exclude some percentage of the values when taking the mean. It goes across the first non-singleton dimension.
        aThresholdMag = abs(aThresholdFactor .* vTrackTrimmedMean ./ timePerFrame);
        aTrackMag = abs(diff(vTrack, 1) ./ timePerFrame); % Accelerations of the track
        if any(aTrackMag > aThresholdMag, 'all') % Remove the track if the acceleration constraint is violated
            bVelocityTestMSmoothedKFConstrainedMMS{n}(tn, :, :) = [];
            trackAlreadyDeleted = true;
        end

        % Direction constraint
        if ~trackAlreadyDeleted
            angleTrack = atan2(vTrack(:, 2), vTrack(:, 1));
            angleTrackChanges = diff(angleTrack);
            if any(abs(angleTrackChanges) > angleChangeThreshold)
                bVelocityTestMSmoothedKFConstrainedMMS{n}(tn, :, :) = [];
            end
        end

    end
end

%% Plot density map with the paired bubbles after persistence
bSum = zeros(img_size(1), img_size(2), img_size(3));

for n = startFrame:length(bVelocityTestM)
% for n = startFrame
    tempBuf = bVelocityTestM{n};
    for tn = 1:size(tempBuf, 1) % track number
%     for tn = 1
        trackTemp = squeeze(tempBuf(tn, :, :))';
        bSum(trackTemp(1, 1), trackTemp(1, 2), trackTemp(1, 3)) = bSum(trackTemp(1, 1), trackTemp(1, 2), trackTemp(1, 3)) + 1;
        for iti = 1:size(trackTemp, 1) % inside track index
            bSum(trackTemp(iti, 4), trackTemp(iti, 5), trackTemp(iti, 6)) = bSum(trackTemp(iti, 4), trackTemp(iti, 5), trackTemp(iti, 6)) + 1;
        end
    end
end

clear n tn iti trackTemp tempBuf

volumeViewer(bSum .^ 0.5)
%% Plot speed map after persistence with linear interpolation, on the cleaned and refined velocity data
speedMap = zeros(img_size(1), img_size(2), img_size(3));
plotPower = 1;

% Counters for proper averaging if there are overlapped pixels from
% different tracks
speedMapCounter = zeros(img_size(1), img_size(2), img_size(3));

for ti = startFrame:size(bVelocityTestMSmoothed, 1)
% for ti = startFrame
%     bvTemp = bVelocityTestMSmoothedKFConstrainedMMS{ti};
    bvTemp = bVelocityTestMSmoothedKFMMS{ti};
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
