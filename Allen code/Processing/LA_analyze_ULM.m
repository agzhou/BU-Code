%% Description
% ULM data analysis:
% Take reconstructed data from individual buffers/batches/superframes, each with
% some number of frames/subframes, each containing N acquisitions/angles.

% Re-sample with a rolling method to get more effective frames.

% SVD to separate the bubble signals from the tissue signal and other
% clutter

% ...
%% Use parallel processing for speed
% https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox

pp = gcp('nocreate');
if isempty(pp)
    % There is no parallel pool
    parpool LocalProfile1

end

%% Load parameters and make folder for saving the processed data
datapath = 'F:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\';
% datapath = 'D:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\';
load([datapath, 'params.mat'])
saveFolderName = 'Processed Data with NLM and isotropic half wl spacing\';
mkdir([datapath, saveFolderName])
savepath = [datapath, saveFolderName];
extHDsavepath = 'K:\Allen data backup\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processed Data\';

filename_structure = [P.Trans.name, '-IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.PRF), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-'];
%% Parameters for processing the data
% Define various processing parameters
% Singular value thresholds
sv_threshold_lower = 10;
sv_threshold_upper = 80;

% Region of interest (HARD CODING FOR NOW)
zrange = 40:120;
% xrange = 1:128;
xrange = 1:227;
% framerange = 1:size(IQf, 3);
% range = {zrange, xrange, framerange};
range = {zrange, xrange};

% Image refinement and localization parameters
imgRefinementFactor = [10, 10]; % z, x pixel refinement factor
XCThreshold = 0.4;
areaThreshold = 3;

% Load and refine simulated PSF
% load('F:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\PSF sim\PSF.mat')
% % load('D:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\PSF sim\PSF.mat')
% PSFs = PSF(90:110, 58:71); % PSF section, hard code this for now
% refPSF = imresize(PSFs, [size(PSFs, 1) * imgRefinementFactor(1), size(PSFs, 2) * imgRefinementFactor(2)], 'bilinear');

load('F:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\PSF sim\PSF_halfwl.mat')
PSFs = PSF(90:110, 108:120); % PSF section, hard code this for now
refPSF = imresize(PSFs, [size(PSFs, 1) * imgRefinementFactor(1), size(PSFs, 2) * imgRefinementFactor(2)], 'bilinear');

% [~, refPSF_center] = max(abs(refPSF), [], 'all');
% refPSF = refPSF(190:210, 118:138);

allCenters = {};

%% Make video of the filtered IQ data with bubbles
vo = VideoWriter([savepath, 'video_10_80']); % video object

vo.Quality = 100;
vo.FrameRate = 30;
open(vo);

vf = figure;
colormap gray
findfigs

% The resizing to account for the XC padding is not perfect
zOffset = size(refPSF, 1);
if mod(zOffset, 2) ~= 0
    zOffset = zOffset + 1;
end
xOffset = size(refPSF, 2);
if mod(xOffset, 2) ~= 0
    xOffset = xOffset + 1;
end

% test = tc(zOffset/2 : size(tc, 1) - zOffset/2, xOffset/2 : size(tc, 2) - xOffset/2);
% need to do this more rigorously
xCorrection = 3;
zCorrection = -10;

%% Process the data
tic
numFiles = 315;
fileGlobalIndex = 1;

for filenum = 1:numFiles
% for filenum = 1:1
%     load([datapath, 'IQ data\', filename_structure, num2str(filenum), '.mat'])  % load each reconstructed buffer/batch/superframe
    load([datapath, 'IQ data all half wl gain -0.5\', filename_structure, num2str(filenum), '.mat'])  % load each reconstructed buffer/batch/superframe
%     IQr = LA_rollingFrames(IQ);                                                 % rolling method to get more effective frames
    IQr = squeeze(sum(IQ, 3)); % coherent sum across angles
    if filenum == 1
        [zp, xp, nf] = size(IQr);
        range{3} = 1:nf; % set frame range after rolling on the first file
        
    end

    % SVD proc part 1
%     tic
    [PP, EVs, V_sort] = getSVs1D(IQr);
%     disp('SVs decomposed')
%     toc
    % SVD proc part 2
%     tic
    [IQf] = applySVs1D(IQr, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
%     disp('SVD filtered images put together')
%     save([savepath, 'Filtered-Data-', num2str(filenum)], 'IQr', 'PP', 'EVs', 'V_sort', 'IQf', "-v6")

    % Denoise with NLM filter
    IQf_dn = zeros(size(IQf));
    for fr = 1:size(IQf, 3)
        IQf_dn(:, :, fr) = imnlmfilt(abs(IQf(:, :, fr)));
    end
    [centers, refIQs, XC] = localizeBubbles2D(IQf_dn, refPSF, range, imgRefinementFactor, XCThreshold, areaThreshold);

%     [centers, refIQs, XC] = localizeBubbles2D_new(IQf, refPSF, range, imgRefinementFactor, XCThreshold, areaThreshold);

%     save([savepath, 'dataproc-', num2str(filenum)], 'IQf', 'refIQs', 'XC', "-v6")
%     save([savepath, 'centers-', num2str(filenum)], 'centers', "-v6")
    save([savepath, 'dataproc-', num2str(filenum)], 'centers', 'IQf', "-v6")

    allCenters{filenum} = centers;
    disp(strcat("Centroid finding done: file ", num2str(filenum)))

    % Write the frames to the video
    for f = 1:size(refIQs, 3)
        td = abs(refIQs(:, :, f)); % temp, data for 1 frame
        tc = centers(:, :, f); % temp, centers for 1 frame
        
        tcOffset = tc(zOffset/2 + zCorrection: size(tc, 1) - zOffset/2 + zCorrection, xOffset/2 + xCorrection : size(tc, 2) - xOffset/2 + xCorrection);
        
        figure(vf)
%         findfigs %%%%%%%%%%
        imagesc(td)
        hold on
        spy(tcOffset, 'ro') % Plot centers on top
        hold off
        
        % add title
        title(strcat("Frame ", num2str(fileGlobalIndex)))

        axis square
        axis tight

        cp = getframe(vf);     % get the current plane
        rgb = frame2im(cp);      % convert the frame to rgb data
    
        writeVideo(vo, rgb);
        
        fileGlobalIndex = fileGlobalIndex + 1;
    end

    toc
end
save([savepath, 'proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'PSF', 'range', 'imgRefinementFactor', 'XCThreshold', 'areaThreshold')
toc

save([savepath, 'allCenters.mat'], 'allCenters')

close(vo); % Close VideoWriter

%% For LOTUS: turn IQ to coherently summed IQ
% Load parameters and make folder for saving the processed data
datapath = 'F:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\IQ data all half wl gain -0.5\';
load(['F:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\params.mat'])
saveFolderName = 'Coherently summed IQ\';
mkdir([datapath, saveFolderName])

filename_structure = [P.Trans.name, '-IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.PRF), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-'];

numFiles = 315;
for filenum = 1:numFiles
% for filenum = 1:1
%     load([datapath, 'IQ data\', filename_structure, num2str(filenum), '.mat'])  % load each reconstructed buffer/batch/superframe
    load([datapath, filename_structure, num2str(filenum), '.mat'])  % load each reconstructed buffer/batch/superframe
%     IQr = LA_rollingFrames(IQ);                                                 % rolling method to get more effective frames
    IQ = squeeze(sum(IQ, 3)); % coherent sum across angles
    save([datapath, saveFolderName, 'IQcs-', num2str(filenum)], 'IQ', "-v6")

end
%% Test plotting the centroids on top of the filtered IQ data
% figure; imagesc(XC(:, :, 1)); hold on; spy(centers(:, :, 1), 'ro'); hold off

figure;
plotTestInd = 1;
td = abs(refIQs(:, :, plotTestInd)); % test data
tc = centers(:, :, plotTestInd); % test centers

% zOffset = size(XC, 1) - size(td, 1);
zOffset = size(refPSF, 1);
if mod(zOffset, 2) ~= 0
    zOffset = zOffset + 1;
end
% xOffset = size(XC, 2) - size(td, 2);
xOffset = size(refPSF, 2);
if mod(xOffset, 2) ~= 0
    xOffset = xOffset + 1;
end

% test = tc(zOffset/2 : size(tc, 1) - zOffset/2, xOffset/2 : size(tc, 2) - xOffset/2);
% need to do this more rigorously
xCorrection = 3;
zCorrection = -10;
tcOffset = tc(zOffset/2 + zCorrection: size(tc, 1) - zOffset/2 + zCorrection, xOffset/2 + xCorrection : size(tc, 2) - xOffset/2 + xCorrection);

imagesc(td)
% imagesc(abs(XC(:, :, plotTestInd)))
hold on
spy(tcOffset, 'ro') % Plot centers on top
hold off
%% Plot the centroid density map

centerSum = sum(allCenters{1}, 3);
for ci = 2:length(allCenters)
    centerSum = centerSum + sum(allCenters{ci}, 3);
end

%% Calculate bubble count
bubbleCount = zeros(length(allCenters) * size(allCenters{1}, 3), 1); % numFiles/# buffers x # frames per buffer. Count of bubbles in each frame
mbci = 1; % microbubble count index
for ci = 1:length(allCenters)
    bufTemp = allCenters{ci};

    for f = 1:size(bufTemp, 3)
        bubbleCount(mbci) = sum(bufTemp(:, :, f), 'all');
        mbci = mbci + 1;
    end
end
totalCount = sum(centerSum, 'all');

figure; plot(1:mbci-1, bubbleCount)
%%
figure; imagesc(centerSum); colormap turbo

%%
test = centerSum .^ 0.5;
figure; imagesc(test); colormap turbo

%% remove rectangular regions of the test plot
d = drawrectangle
rmvp = round(d.Position);
test(rmvp(2) : rmvp(2) + rmvp(4), rmvp(1) : rmvp(1) + rmvp(3)) = 0;
imagesc(test); colormap turbo
%%
test = centerSum; test(test > 40) = 40;
figure; imagesc(test(117:end, :)); colormap hot
% zpts = [];
% xpts = [];
% 
% for f = 1:size(allCentroids, 1)
%     zpts = [zpts; allCentroids{f}(:, 1)];
%     xpts = [xpts; allCentroids{f}(:, 2)];
% end
% 
% hPixFactor = 10; % increase the pixel count by this factor in each dimension
% figure;
% h = histogram2(zpts, xpts, [zp * hPixFactor, xp * hPixFactor], 'DisplayStyle','tile');
% grid off
% colormap hot
