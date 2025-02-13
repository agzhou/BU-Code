% Acknowledgement: using Jean-Yves Tinevez's simpletracker as a reference
%%
addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Processing\munkres')
load('D:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\params.mat')
% [add a line to load proc_params.mat]

if ~exist('allCenters', 'var')
    load('D:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processed Data\allCenters_02_10_2025_SVs10to80_0p4XCThreshold_10refinementfactor.mat')
end
%% Define some threshold parameters
% res_assumption = 10; % assumption of how good the localization is (um)
% maxDetectSpeed = (res_assumption / 1e6) / (1 / P.frameRate); % maximum bubble speed we can detect
% maxDistPerFrame = (res_assumption / 1e6); % maximum travel distance between frames (m)

%%
totalFrames = length(allCenters) * size(allCenters{1}, 3);

% allCentersE = cell(totalFrames, 1); % elongated (take the frames from each buffer out individually)
% 
% caiGlobal = 1;
% for cai = 1:length(allCenters) % cell array index
%     for bfi = 1:size(allCenters{cai}, 3) % buffer frame index
%         allCentersE{caiGlobal} = allCenters{cai}(:, :, bfi);
%         caiGlobal = caiGlobal + 1;
%     end
%     
% end
% 
% clear allCenters % save memory

% figure; spy(currentFrame, 'ro')
% hold on; spy(nextFrame, 'b*'); hold off
%% Turn the logical matrices into coordinates
% load('D:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processed Data\allCenters_02_05_2025_SVs10to80_files1to315.mat')
totalFrames = length(allCenters) * size(allCenters{1}, 3);

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

img_size = size(allCenters{1});
clear bfi cai caiGlobal allCenters xc zc

%% Calculate bubble count
bubbleCount = zeros(length(centerCoords), 1); % numFiles/# buffers x # frames per buffer. Count of bubbles in each frame
centerCoords_corrected = centerCoords;
% mbci = 1; % microbubble count index
parfor fi = 1:length(centerCoords) % frame index
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

%% Max speed (distance per frame) threshold and initialize variables
maxSpeedExpectedMMPerS = 100;   % max expected flow speed [mm/s]
timePerFrame = 1 / P.frameRate; % time elapsed per frame [s]
maxDistPerFrameM = (maxSpeedExpectedMMPerS / 1000) * timePerFrame;
pixelsPerM = 2 / P.wl * imgRefinementFactor(1);
% maxPixelDistPerFrame = maxDistPerFrameM / (P.wl/2/imgRefinementFactor(1)); % HARD CODE THIS FOR TESTING NOW
maxPixelDistPerFrame = maxDistPerFrameM * pixelsPerM;

pairedBubbles = cell(totalFrames - 1, 1);
costs = zeros(totalFrames - 1, 1); % cost vector for each frame pairing

%% Calculate pairing
testFig = figure;
% parfor f = 1:totalFrames - 1
for f = 17
    currentFrame = centerCoords_corrected{f};
    nextFrame = centerCoords_corrected{f + 1};



    nbC = size(currentFrame, 1); % number of bubbles in the current frame
    nbN = size(nextFrame, 1);    % number of bubbles in the next frame
    D = NaN(nbC, nbN);

    if nbC > 1 && nbN > 1
        % Go through each current frame point and get the distance to each
        % point in the next frame. Store in the distance matrix D.
        for cpi = 1:nbC % current point index
    %     for cpi = 1
            currentPoint = currentFrame(cpi, :);    % z and x coords of the current frame's point cpi
            d = nextFrame - currentPoint;           % vectorized difference between the z and x coords of all the points in the next frame and point cpi from the current frame
            D(cpi, :) = sqrt(sum((d .^ 2), 2));     % distance formula on the above
    
            D(D > maxPixelDistPerFrame) = Inf; 
        end
    
        [assignment, cost] = munkres(D);  % Use the munkres.m function from Yi Cao (Matlab File Exchange)
        costs(f) = cost;                  % update the overall cost vector
    
        [pzct, pxct] = find(assignment); % temporarily store the z and x coords from the paired point as stored in the "assignment" matrix
        pairedBubbles{f} = [pzct, pxct];
    
%         % Plot a line for each pair
%       figure(testFig); scatter(currentFrame(:, 1), currentFrame(:, 2), 'ro')
%       hold on; scatter(nextFrame(:, 1), nextFrame(:, 2), 'b*');
%       title(strcat("Frame ", num2str(f)))
%         for npb = 1:size(pairedBubbles{f}, 1)
%             pairInd = pairedBubbles{f}(npb, :);
%             ZVecTemp = [currentFrame(pairInd(1), 1), nextFrame(pairInd(2), 1)];
%             XVecTemp = [currentFrame(pairInd(1), 2), nextFrame(pairInd(2), 2)];
%             line(ZVecTemp, XVecTemp, 'LineWidth', 1.5)
%         end
        
    end
%     hold off
end

%% Plot density map with the paired bubbles only
% bSum = zeros(size(allCenters{1}, 1), size(allCenters{1}, 2));
bSum = zeros(img_size(1), img_size(2));

bSumFig = figure; colormap turbo
% hold on
figure(bSumFig)
bsIm = imagesc(bSum);
% bSumFig.XDataSource = 

for f = 1:length(pairedBubbles)
% for f = 1:10000
    disp(num2str(f))
    if ~isempty(pairedBubbles{f})
        coordsCF = centerCoords{f}; % coords for the current frame
        coordsNF = centerCoords{f + 1}; % coords for the current frame
    %     for npb = 1:size(pairedBubbles{f}, 1)

        keepBubblesCF = pairedBubbles{f}(:, 1); % which bubbles (indices) to keep for the current frame. Could do it more efficiently or correctly
        for kbcfi = keepBubblesCF' % index for each kept bubble in the current frame
            kbcfCoord = coordsCF(kbcfi, :);
            bSum(kbcfCoord(1), kbcfCoord(2)) = bSum(kbcfCoord(1), kbcfCoord(2)) + 1; % Update the count for that pixel
        end

        keepBubblesNF = pairedBubbles{f}(:, 2); % which bubbles (indices) to keep for the next frame. Could do it more efficiently or correctly
        for kbnfi = keepBubblesNF' % index for each kept bubble in the next frame
            kbnfCoord = coordsNF(kbnfi, :);
            bSum(kbnfCoord(1), kbnfCoord(2)) = bSum(kbnfCoord(1), kbnfCoord(2)) + 1; % Update the count for that pixel
        end

    end
end
hold off

%% Plot pairing
testFig = figure;
% for f = 1:totalFrames - 1
for f = 1
% for f = 3
    currentFrame = centerCoords{f};
    nextFrame = centerCoords{f + 1};

%     pairedBubbles{f} = [pzct, pxct];
    
    % Plot a line for each pair
    figure(testFig); scatter(currentFrame(:, 1), currentFrame(:, 2), 'ro')
    hold on; scatter(nextFrame(:, 1), nextFrame(:, 2), 'b*');
    title(strcat("Frame ", num2str(f)))
    for npb = 1:size(pairedBubbles{f}, 1)
        pairInd = pairedBubbles{f}(npb, :);
        ZVecTemp = [currentFrame(pairInd(1), 1), nextFrame(pairInd(2), 1)];
        XVecTemp = [currentFrame(pairInd(1), 2), nextFrame(pairInd(2), 2)];
        line(ZVecTemp, XVecTemp, 'LineWidth', 1.5)
    end
        
%     hold off
end
legend('Frame n', 'Frame n + 1')