% Acknowledgement: using Jean-Yves Tinevez's simpletracker as a reference
%% Add dependencies and load parameters + bubble center locations
datapath = 'D:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processed Data\with NLM\isotropic half wl pix spacing\';
addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Processing\munkres')
addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Previous lab code\A-US-ULM\SubFunctions\')
load('D:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\params.mat')

% Load localization processing parameters: proc_params.mat
load([datapath, 'proc_params.mat'])

% [add line to load the recon pixel spacing .mat file]
% For now, hard coding and assuming equal z and x spacing
pix_spacing = P.wl/2;

% Find and load the localized centers file, which starts with 'allCenters'
if ~exist('allCenters', 'var') % Only load if the variable doesn't already exist in the workspace
    allCenters_generalpath = fullfile(datapath, 'allCenters*.mat');
    allCentersDir = dir(allCenters_generalpath);
    load([allCentersDir.folder, '\', allCentersDir.name])
end

clear allCenters_generalpath allCentersDir
%% Turn the logical matrices into coordinates
totalFrames = length(allCenters) * size(allCenters{1}, 3);

% Cell array with an entry for each frame. Each entry contains (# bubbles) of coordinate pairs (z, x) of the detected bubble centers
centerCoords = cell(totalFrames, 1); 

caiGlobal = 1;
for cai = 1:length(allCenters)           % cell array index
    for bfi = 1:size(allCenters{cai}, 3) % buffer frame index
        [zc, xc] = find(allCenters{cai}(:, :, bfi));
        centerCoords{caiGlobal} = [zc, xc];
%         centerCoords{caiGlobal}  = find(allCenters{cai}(:, :, bfi));
        caiGlobal = caiGlobal + 1;
    end
    
end

img_size = size(allCenters{1}); % Save image size if we want to clear allCenters
clear bfi cai caiGlobal xc zc

%% Calculate bubble count
bubbleCount = zeros(length(centerCoords), 1);   % numFiles/# buffers x # frames per buffer. Count of bubbles in each frame
centerCoords_corrected = centerCoords;          % Correct the centerCoords because some frames have every pixel identified as        a bubble
% mbci = 1; % microbubble count index
parfor fi = 1:length(centerCoords_corrected) % frame index
% for fi = 1
    bufTemp = centerCoords{fi};

    % Some error that makes every pixel a bubble
    if size(bufTemp, 1) >= img_size(1) * img_size(2)
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
figure; plot(1:length(bubbleCount), bubbleCount, '.')
title('Bubble count')
xlabel('Frame number')
ylabel('Bubble count')

clear fi bufTemp
%% Max speed (distance per frame) threshold and initialize variables
maxSpeedExpectedMMPerS = 50;                                        % max expected flow speed [mm/s]
timePerFrame = 1 / P.frameRate;                                     % time elapsed per frame [s]
maxDistPerFrameM = (maxSpeedExpectedMMPerS / 1000) * timePerFrame;  % max distance traveled per frame [m], according to the max expected flow speed and frame rate
pixelsPerM = 1 / pix_spacing * imgRefinementFactor(1);              % # of pixels per meter, which depends on the pixel spacing from reconstruction and the image refinement factor from the localization
maxPixelDistPerFrame = maxDistPerFrameM * pixelsPerM;               % max distance traveled per frame in units of pixels

bubblePairs = cell(totalFrames - 1, 1);   % Initialize cell vector of paired bubble indices

%% Calculate pairing
% testFig = figure;

ubS = cell(totalFrames - 1, 1);             % unassigned bubbles from the source frames
ubT = cell(totalFrames - 1, 1);             % unassigned bubbles from the target frames

parfor f = 1:totalFrames - 1
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

clear f nbS nbT D spi assignment unassignedrows unassignedcolumns
%% Create tracks with persistence (might have a bug)
pers = 5; % # of frames a track needs to persist through to keep it

% Separate the pairs of coordinates so we can change their sizes independently
bubblePairsPers = cell(length(bubblePairs), 2); % bubble pairs with the pairs separated into another cell dimension
for n = 1:length(bubblePairs)
    if ~isempty(bubblePairs{n})
        bubblePairsPers{n, 1} = bubblePairs{n}(:, 1);
        bubblePairsPers{n, 2} = bubblePairs{n}(:, 2);
    end
end

% Store the bubble indices that contribute to each track of at least (pers) # of frames
tracks = cell(size(bubblePairsPers, 1) - pers, 1);
for n = 1:length(bubblePairsPers) - pers
% for n = 1:1
    bubblePairsPersTemp = bubblePairsPers;
    for pfc = 1:pers - 1 % persistence frame count
%     for pfc = 1
        startIndex = bubblePairsPersTemp{n + pfc - 1, 2};    % indices of the paired "target" bubbles in frame n + 1 (pair n), which will be sorted in ascending order
        endIndex = bubblePairsPersTemp{n + pfc, 1};  % indices of the paired "source" bubbles in frame n + 2 (pair n + 1), not necessarily in ascending order because it's aligned with the "target" indices in frame n + 1
        [trackContinuesIndices, is, ie] = intersect(startIndex, endIndex); % find the common values in the start and end vectors, and the corresponding indices for each
        
        % Update the paired and tracked list for this immediate set of pairs
        bubblePairsPersTemp{n + pfc - 1, 1} = bubblePairsPersTemp{n + pfc - 1, 1}(is);
        bubblePairsPersTemp{n + pfc - 1, 2} = bubblePairsPersTemp{n + pfc - 1, 2}(is);
        bubblePairsPersTemp{n + pfc, 1} = bubblePairsPersTemp{n + pfc, 1}(ie); % I think this preserves the order
        bubblePairsPersTemp{n + pfc, 2} = bubblePairsPersTemp{n + pfc, 2}(ie); % Need to preserve the order, so set the next starting point????
        
        % Go back and update the previous pairs too
%         for recpfc = pfc - 1 : -1 : 1   % recursive persistence frame count
        for recpfc = 1:pfc - 1   % recursive persistence frame count
            recStartIndex = bubblePairsPersTemp{n + pfc - recpfc - 1, 2};    % indices of the paired "target" bubbles in frame n + 1 (pair n), which will be sorted in ascending order
            recEndIndex = bubblePairsPersTemp{n + pfc - recpfc, 1};  % indices of the paired "source" bubbles in frame n + 2 (pair n + 1), not necessarily in ascending order because it's aligned with the "target" indices in frame n + 1
            [recTrackContinuesIndices, ris, rie] = intersect(recStartIndex, recEndIndex); % find the common values in the start and end vectors, and the corresponding indices for each
        
            bubblePairsPersTemp{n + pfc - recpfc - 1, 1} = bubblePairsPersTemp{n + pfc - recpfc - 1, 1}(ris);
            bubblePairsPersTemp{n + pfc - recpfc - 1, 2} = bubblePairsPersTemp{n + pfc - recpfc - 1, 2}(ris);
        end
    end
    tracks{n} = bubblePairsPersTemp(n : n + pfc, :);
end

clear bubblePairsPersTemp n pfc recpfc

%% Create tracks with persistence - hopefully fixing the bug
pers = 5; % # of frames a track needs to persist through to keep it

tic
% Separate the pairs of coordinates so we can change their sizes independently
bubblePairsPers = cell(length(bubblePairs), 2); % bubble pairs with the pairs separated into another cell dimension
for n = 1:length(bubblePairs)
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

% clear ti fn tracksTemp startPoints endPoints vfn zcInterp xcInterp i vmapTempi vmapTemp
%% Combine cleaned tracks old
tracksCombined = cell(size(tracksClean));
% for n = 1:size(tracksClean, 1) % go through frames
for n = 1
    tracksntemp = tracksClean{n};            % tracks starting in frame n
    nbifn = nbifAll(n);                      % # of bubbles in the tracks starting in frame n
%     tcTemp = zeros(1, 4);
    tcTemp = tracksntemp(1 : nbifn, :);
%     for np = 1:pers             % go through all the possible overlaps in frames
    for np = 1
        tracksnptemp = tracksClean{n + np};  % tracks starting in frame n + np
        nbifnp = nbifAll(n + np);            % # of bubbles in the tracks starting in frame n + np
%         for nfi = 1:pers
    
        for nfi = 1 % frame number inside tracknptemp
            tracksntempfni = tracksntemp((nfi) * nbifn + 1 : (nfi + 1) * nbifn, :);
            tracksnptempfni = tracksnptemp((nfi - 1) * nbifnp + 1 : nfi * nbifnp, :);
%             tcTemp = [tcTemp; unique([tracksntempfni; tracksnptempfni], 'rows', 'stable')];
            tcTempnp = [tcTempnp; unique([tracksntempfni; tracksnptempfni], 'rows', 'stable')];
        end
    end
end
% clear trackntemp tracknptemp tcTemp n np

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

%% Velocity map after doing persistence

bVelocityPers = cell(size(tracks, 1), pers); % bubble velocity - both components [z velocity, x velocity]
for ti = 1:length(tracks)% track index
% for ti = 1
        
    tracksTemp = tracks{ti};
    for tli = 1:size(tracksTemp, 1)% track length index
        if ~isempty(tracksTemp{tli})
            startPointInd = tracksTemp{tli, 1};
            endPointInd = tracksTemp{tli, 2};

            coordsSFTemp = centerCoords_corrected{ti + tli - 1}(startPointInd, :);  % coords for the source frame
            coordsTFTemp = centerCoords_corrected{ti + tli}(endPointInd, :);        % coords for the target frame
            vTemp = (coordsTFTemp - coordsSFTemp) ./ timePerFrame; % THIS IS IN PIXELS per second, NOT DISTANCE per second!!!!!!!!!!!!!!!!!
            bVelocityPers{ti, tli} = vTemp;
        else
            bVelocityPers{ti, tli} = NaN;
        end
    end
end
clear ti tli tracksTemp startPointInd endPointInd coordsSFTemp coordsTFTemp vTemp

%% Refine the velocity map - (maybe remove tracks where the position is oscillating over time)
bVelocityTestMSmoothed = bVelocityTestM;
for n = 1:size(bVelocityTestM, 1)
% for n = 1
    vt = bVelocityTestM{n};
    vtSmoothed = vt;
    for tn = 1:size(vt, 1) % track number
%     for tn = 1
        tnzVel = squeeze(vt(tn, 5, :)); % get the z velocity at each point in track tn
        tnzVelSmoothed = movmean(tnzVel, 3);
        vtSmoothed(tn, 5, :) = tnzVelSmoothed;

        tnxVel = squeeze(vt(tn, 6, :)); % get the x velocity at each point in track tn
        tnxVelSmoothed = movmean(tnxVel, 3);
        vtSmoothed(tn, 6, :) = tnxVelSmoothed;
        
    end
    bVelocityTestMSmoothed{n} = vtSmoothed;
end
clear n tn vt vtSmoothed tnzVel tnzVelSmoothed vtSmoothed tnxVel tnxVelSmoothed
% % test
% testMap = zeros(img_size(1), img_size(2));
% bv1 = bVelocity{1, 1};
% for i = size(bv1, 1)
%     bv1temp = bv1(i, :);
%     testMap(bv1temp(1), bv1temp(2)) = testMap(bv1temp(1), bv1temp(2)) + bv1temp(3);
% end
% figure; imagesc(testMap)
%% Plot density map with the paired bubbles only
% bSum = zeros(size(allCenters{1}, 1), size(allCenters{1}, 2));
bSum = zeros(img_size(1), img_size(2));

bSumFig = figure; colormap turbo
% hold on
figure(bSumFig)

% bSumFig.XDataSource = 

for f = 1:length(bubblePairs)
% for f = 1:10000
    disp(num2str(f))
    if ~isempty(bubblePairs{f})
        coordsCF = centerCoords_corrected{f}; % coords for the current frame
        coordsNF = centerCoords_corrected{f + 1}; % coords for the current frame
    %     for npb = 1:size(pairedBubbles{f}, 1)

        keepBubblesCF = bubblePairs{f}(:, 1); % which bubbles (indices) to keep for the current frame. Could do it more efficiently or correctly
        for kbcfi = keepBubblesCF' % index for each kept bubble in the current frame
            kbcfCoord = coordsCF(kbcfi, :);
            bSum(kbcfCoord(1), kbcfCoord(2)) = bSum(kbcfCoord(1), kbcfCoord(2)) + 1; % Update the count for that pixel
        end

        keepBubblesNF = bubblePairs{f}(:, 2); % which bubbles (indices) to keep for the next frame. Could do it more efficiently or correctly
        for kbnfi = keepBubblesNF' % index for each kept bubble in the next frame
            kbnfCoord = coordsNF(kbnfi, :);
            bSum(kbnfCoord(1), kbnfCoord(2)) = bSum(kbnfCoord(1), kbnfCoord(2)) + 1; % Update the count for that pixel
        end

    end
end
hold off
bsIm = imagesc(bSum);

%% Plot z velocity map without interpolation
velocityMap = zeros(img_size(1), img_size(2));

mymap = [zeros(256,2),linspace(1,0,256)';linspace(0,1,256)', zeros(256,2)];

vMapFig = figure; colormap(vMapFig, mymap)
% hold on
figure(vMapFig)
for n = 1:length(bubblePairs) - 1
% for n = 1
    zVel = bVelocity{n}(:, 1);
    % Use the source coordinates as the point to plot velocity?
    % OR FILL IN THE LINE BETWEEN????
    if ~isempty(bubblePairs{n})
        coords = centerCoords_corrected{n}(bubblePairs{n}(:, 1), :); % coords for the source frame
        for ci = 1:length(zVel) % coordinate index
            velocityMap(coords(ci, 1), coords(ci, 2)) = velocityMap(coords(ci, 1), coords(ci, 2)) + zVel(ci); % add the velocities to the existing values
        end
    end
end

imagesc(velocityMap)
% clim([])

%% Plot z velocity map after persistence with linear interpolation
% zVelocityPersMap = zeros(img_size(1), img_size(2));
zvUpPersMap = NaN(img_size(1), img_size(2));
zvDownPersMap = NaN(img_size(1), img_size(2));

% Counters for proper averaging if there are overlapped pixels from
% different tracks
zvUpPersMapCounter = ones(img_size(1), img_size(2));
zvDownPersMapCounter = ones(img_size(1), img_size(2));

% for ti = 1:length(tracks) % track index
for ti = 1000:15000
% for ti = 1:15517
    tracksTemp = tracks{ti};
    for pfi = 1:pers % persistence frame index
        zVel = bVelocityPers{ti, pfi}(:, 1); % z component (the :, 1) of velocity between frames ti + pfi and ti + pfi + 1 (I think?)
        if ~isnan(zVel)
%         if ~isnan(zVel) & (size(zVel, 1) > 1) % TEST
            % Get the coordinates that the velocity values correspond to
            startPointInd = tracksTemp{pfi, 1};
            endPointInd = tracksTemp{pfi, 2};
            coordsSFTemp = centerCoords_corrected{ti + pfi - 1}(startPointInd, :);  % coords for the source frame
            coordsTFTemp = centerCoords_corrected{ti + pfi}(endPointInd, :);        % coords for the target frame
            
            % Linearly interpolate to fill in the pixels between the
            % coordinates of the paired bubbles.
            % Go one by one through the zVel vector
            for pri = 1:length(zVel)    % pair index
                zVelTemp = zVel(pri);
                if zVelTemp > 0         % up z flow
                    [zcInterp, xcInterp] = ULM_interp2D(coordsSFTemp(pri, :), coordsTFTemp(pri, :));
    %                 zvUpMap([zcInterp, xcInterp]) = zvUpMap([zcInterp, xcInterp]) + zVelTemp;
                    for nip = 1:length(zcInterp) % number of interp points  
                        if isnan(zvUpPersMap(zcInterp(nip), xcInterp(nip)))
                            zvUpPersMap(zcInterp(nip), xcInterp(nip)) = zVelTemp; % If this is the first point entered at that pixel, use the value
                        else
    %                         zvUpPersMap(zcInterp(nip), xcInterp(nip)) = (zvUpPersMap(zcInterp(nip), xcInterp(nip)) + zVelTemp) / 2; % If it's not the first point entered at that pixel, average between the existing and new values
                            zvUpPersMap(zcInterp(nip), xcInterp(nip)) = zvUpPersMap(zcInterp(nip), xcInterp(nip)) + zVelTemp; % If it's not the first point entered at that pixel, add the new value
                            zvUpPersMapCounter(zcInterp(nip), xcInterp(nip)) = zvUpPersMapCounter(zcInterp(nip), xcInterp(nip)) + 1; % increment the counter so we can average later
                        end
                    end
                else                % down z flow
                    [zcInterp, xcInterp] = ULM_interp2D(coordsSFTemp(pri, :), coordsTFTemp(pri, :));
    %                 zvDownMap([zcInterp, xcInterp]) = zvDownMap([zcInterp, xcInterp]) + zVelTemp;
                    for nip = 1:length(zcInterp) % number of interp points  
                        if isnan(zvDownPersMap(zcInterp(nip), xcInterp(nip)))
                            zvDownPersMap(zcInterp(nip), xcInterp(nip)) = zVelTemp; % If this is the first point entered at that pixel, use the value
                        else
                            zvDownPersMap(zcInterp(nip), xcInterp(nip)) = zvDownPersMap(zcInterp(nip), xcInterp(nip)) + zVelTemp; % If it's not the first point entered at that pixel, add the new value
                            zvDownPersMapCounter(zcInterp(nip), xcInterp(nip)) = zvDownPersMapCounter(zcInterp(nip), xcInterp(nip)) + 1; % increment the counter so we can average later
    
                        end
                    end
                end
            end
        end
    end
end
zvUpPersMap = zvUpPersMap ./ zvUpPersMapCounter;
zvDownPersMap = zvDownPersMap ./ zvDownPersMapCounter;
% Plot z velocity map after persistence
clear zVel

zvDownPersMap(isnan(zvDownPersMap)) = 0;
zvUpPersMap(isnan(zvUpPersMap)) = 0;

%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%
% zvDownPersMap(zvDownPersMapCounter == 1) = 0;
% zvUpPersMap(zvUpPersMapCounter == 1) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Jianbo's colormaps
[VzCmap, VzCmapDn, VzCmapUp, pdiCmapUp, PhtmCmap] = Colormaps_fUS;
zvMapPersFig = figure;

% hold on
vCrange= [-8800, 8800];
figure(zvMapPersFig)
h1 = axes;
imagesc(zvUpPersMap)
% alpha(h1, double(abs(zvUpMap) > 1))
colormap(zvMapPersFig, VzCmap)
caxis(vCrange);
axis tight
colorbar
hold on

h2 = axes;
imagesc(zvDownPersMap)
alpha(h2, double(abs(zvDownPersMap) > 1))
colormap(zvMapPersFig, VzCmap)
caxis(vCrange);
axis tight
colorbar
axis off
linkaxes([h1, h2]);

% clim([])

%% Plot z velocity map after persistence with linear interpolation, on the cleaned and refined velocity data
% zVelocityPersMap = zeros(img_size(1), img_size(2));
zvUpMap = zeros(img_size(1), img_size(2));
zvDownMap = zeros(img_size(1), img_size(2));

% Counters for proper averaging if there are overlapped pixels from
% different tracks
zvUpMapCounter = zeros(img_size(1), img_size(2));
zvDownMapCounter = zeros(img_size(1), img_size(2));

for ti = 1:size(bVelocityTestMSmoothed, 1)
% for ti = 7
%     bvTemp = bVelocityTestMSmoothed{ti}; % get the ti-th entry
    bvTemp = bVelocityTestM{ti}; % get the ti-th entry
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
%             interpPts = ULM_interp2D_linear(coordsStart, coordsEnd, zvTemp); % Get interpolated points with the corresponding z velocity value. each row is [z coord, x coord, z velocity]
            interpPts = ULM_interp2D_linear(coordsStart, coordsEnd, zvTemp, ti); % Get interpolated points with the corresponding z velocity value. each row is [z coord, x coord, z velocity]
            
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
zvUpMask = zvUpMapCounter > 0;
zvDownMask = zvDownMapCounter > 0;
zvUpMap(zvUpMask) = zvUpMap(zvUpMask) ./ zvUpMapCounter(zvUpMask);
zvDownMap(zvDownMask) = zvDownMap(zvDownMask) ./ zvDownMapCounter(zvDownMask);
%

% Plot z velocity map after persistence
clear zVelTemp

% Load Jianbo's colormaps
[VzCmap, VzCmapDn, VzCmapUp, pdiCmapUp, PhtmCmap] = Colormaps_fUS;
zvMapFig = figure;

% hold on
vCrange= [-8800, 8800];
figure(zvMapFig)
h1 = axes;
imagesc(zvUpMap)
% alpha(h1, double(abs(zvUpMap) > 1))
colormap(zvMapFig, VzCmap)
caxis(vCrange);
axis tight
colorbar
hold on

h2 = axes;
imagesc(zvDownMap)
alpha(h2, double(abs(zvDownMap) > 1))
colormap(zvMapFig, VzCmap)
caxis(vCrange);
axis tight
colorbar
axis off
linkaxes([h1, h2]);

% clim([])
