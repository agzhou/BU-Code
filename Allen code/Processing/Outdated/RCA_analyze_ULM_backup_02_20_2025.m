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
datapath = 'G:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\';
load([datapath, 'params.mat'])

numFiles = 96; % define # of files manually for now

IQfolderName = 'IQ Data - Verasonics Recon\'; % 'IQ data\'
saveFolderName = 'Processed Data 02-21-2025\';
mkdir([datapath, saveFolderName])
savepath = [datapath, saveFolderName];
% savepath = 'F:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\run 1 left eye\FMAS Processed Data\';

filename_structure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];
%% Parameters for processing the data
% Define various processing parameters
% Singular value thresholds
sv_threshold_lower = 20;
sv_threshold_upper = 150;

% % Region of interest
xrange = 1:80;
yrange = 1:80;
zrange = 1:142;

% framerange = 1:200;
% framerange = 1:size(IQf, 3);
% range = {xrange, yrange, zrange, framerange};
range = {xrange, yrange, zrange};

% % Image refinement and localization parameters
imgRefinementFactor = [2, 2, 2]; % z, x pixel refinement factor
XCThreshold = 0.4;
volumeThreshold = 4;
% 
% Load and refine simulated PSF
load('G:\Allen\Data\01-29-2025 AZ001 ULM\RC15gV\PSF sim\PSF.mat', 'PSF')
% figure; imagesc(squeeze(abs(PSF(40, :, :)))')

% refPSF = imresize(PSF, [size(PSF, 1) * imgRefinementFactor(1), size(PSF, 2) * imgRefinementFactor(2)], 'bilinear');
% refPSF = refPSF(190:210, 118:138);

% temporary non-refined PSF...
refPSF = PSF(35:46, 35:46, 95:105);
psfFig = figure; psfV = volshow(abs(refPSF));
psfV.BackgroundColor = [1, 1, 1];
% 
allCentroids = {};

%% Process the data
tic
% for filenum = 1:numFiles
for filenum = 25:25
    load([datapath, IQfolderName, filename_structure, num2str(filenum), '.mat'])  % load each reconstructed buffer/batch/superframe
%     IQr = LA_rollingFrames(IQ);                                                 % rolling method to get more effective frames
    
    IQ = squeeze(IData + 1i .* QData);   % Combine I and Q, which are saved separately. It's easier to save the big reconstructed data with savefast, which doesn't support complex values.
    clear IData QData
    
%     if filenum == 1
%         [zp, xp, nf] = size(IQr);
%         range{end + 1} = 1:nf; % set frame range after rolling on the first file
%         
%     end

    % SVD proc part 1
%     tic
    [PP, EVs, V_sort] = getSVs2D(IQ);
    disp('SVs decomposed')
%     toc
    % SVD proc part 2
%     tic
    [IQf] = applySVs2D(IQ, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
    disp('SVD filtered images put together')
%     toc

%     save([savepath, 'Filtered-Data-', num2str(filenum)], 'IQr', 'PP', 'EVs', 'V_sort', 'IQf', "-v6")
    [centroidCoordinates] = localizeBubbles3D(IQf, refPSF, range, imgRefinementFactor, XCThreshold, volumeThreshold);
%     save([savepath, 'IQf-', num2str(filenum)], 'IQf', "-v6")

% %     save([savepath, 'dataproc-', num2str(filenum)], 'IQf', 'centroidCoordinates', "-v6")
%     save([savepath, 'dataproc-', num2str(filenum)], 'centroidCoordinates', "-v6")

    allCentroids = [allCentroids; centroidCoordinates];
    disp(strcat("Centroid finding done: file ", num2str(filenum)))
%     toc
end
% save([savepath, 'proc_params.mat'], 'sv_threshold_lower', 'sv_threshold_upper', 'PSF', 'range', 'imgRefinementFactor', 'binaryThreshold', 'volumeThreshold')
% save([savepath, 'allCentroids'], 'allCentroids', "-v6")
toc
%% Plot the centroid density map

xpts = [];
ypts = [];
zpts = [];
for f = 1:length(allCentroids)
    xpts = [xpts; allCentroids{f}(:, 1)];
    ypts = [ypts; allCentroids{f}(:, 2)];
    zpts = [zpts; allCentroids{f}(:, 3)];
end

figure; scatter3(xpts, ypts, zpts)

% hPixFactor = 10; % increase the pixel count by this factor in each dimension
% figure;
% h = histogram2(zpts, xpts, [zp * hPixFactor, xp * hPixFactor], 'DisplayStyle','tile');
% grid off
% colormap hot











%% for testing: Make video of the filtered IQ data
vo = VideoWriter([savepath, 'IQf_10_500']);
a = abs(IQf);
a = a ./ max(a); % scale 0-1 so it can make a video

vo.Quality = 100;
vo.FrameRate = 120;
open(vo);
for f = 1:size(a, 3)
    writeVideo(vo, a(:, :, f));
end
close(vo);







%%

% IQ_coherent_sum = squeeze(sum(IQ, 3));

IQ_coherent_sum = squeeze(sum(IQstack, 3));
% I_coherent_sum = abs(IQ_coherent_sum); % intensity
% figure; imagesc(I_coherent_sum(:, :, 100))

%% SVD proc part 1
tic
[PP, EVs, V_sort] = getSVs1D(IQ_coherent_sum);
disp('SVs decomposed')
toc
%% SVD proc part 2
sv_threshold_lower = 50;
sv_threshold_upper = 5000;
tic
[IQ_f_50_5000] = applySVs1D(IQ_coherent_sum, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
disp('SVD filtered images put together')
toc
%% Plot filtered data
abs_IQ_f = abs(IQ_f_50_5000);
figure; imagesc(abs_IQ_f(:, :, 1))
figure; imagesc(abs_IQ_f(:, :, 100))
figure; imagesc(abs_IQ_f(:, :, 500))
figure; imagesc(abs_IQ_f(:, :, 1000))
figure; imagesc(abs_IQ_f(:, :, 2000))

%% Plot square of the singular values
figure; plot(EVs, '-o')
title('Eigenvalues')
xlabel('EV #')
ylabel('EV magnitude')

figure; plot(log10(EVs), '-o')
title('Eigenvalues log plot')
xlabel('EV #')
ylabel('log10(EV magnitude)')

%% Make video
vo = VideoWriter('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processing\video_50_5000');
a = abs(IQ_f_50_5000);
a = a ./ max(a); % scale 0-1 so it can make a video

vo.Quality = 100;
vo.FrameRate = 60;
open(vo);
for f = 1:size(a, 3)
    writeVideo(vo, a(:, :, f));
end
close(vo);

%% Localize bubble function test
IQf = IQ_f_50_5000;
load('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\PSF sim\IQ.mat', 'IQ_coherent_sum')
PSF = IQ_coherent_sum;
zrange = 40:120;
xrange = 1:128;
framerange = 1:size(IQf, 3);
range = {zrange, xrange, framerange};
imgRefinementFactor = [2, 2];
XCThreshold = 0.6;
volumeThreshold = 3;
[centroidCoordinates] = localizeBubbles(IQf, PSF, range, imgRefinementFactor, XCThreshold, volumeThreshold);



%% Bubble tracking attempt
zrange = 40:120;
xrange = 1:128;
% frameRange = 1:10000; 
frameRange = 1:size(IQ_f_50_5000, 3);

IQs = IQ_f_50_5000(zrange, xrange, frameRange); % IQ section
% figure; imagesc(abs(IQs(:, :, 1)))
% figure; imagesc(abs(IQs(:, :, 2)))
% test = img(:, :, 2) - img(:, :, 1);
% figure; imagesc(test)

d = diff(IQs, 1, 3); % take first order difference along the frame dimension, seems to get rid of some background
d(:, :, end+1) = IQs(:, :, end) - IQs(:, :, end-2); % add another value so you get back to orig # frames
% figure; imagesc(abs(d(:, :, 1)))
% figure; imagesc(abs(d(:, :, 2)))

%% image refinement?
rfnZ = 2; % refinement pixel increase factor
rfnX = 2;

refIQs = zeros(size(IQs, 1) * rfnZ, size(IQs, 2) * rfnX, size(IQs, 3)); % refined IQ section
% I = IQs(:, :, 1);

tic
% go through all frames and refine
parfor f = 1:size(IQs, 3)
% for f = 1:1
% parfor f = 1:10
    I_temp = IQs(:, :, f);
    refIQs(:, :, f) = imresize(I_temp/max(I_temp, [], 'all') .* 256 .* 5, [size(I_temp, 1) * rfnZ, size(I_temp, 2) * rfnX], 'bilinear');
end
toc

%% and cross correlation?
% if ~exist('base', 'var', 'PSF')
    load('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\PSF sim\IQ.mat') % Load simulated PSF
%     PSF = IQ_coherent_sum(90:110, 55:75);                             % Crop the PSF
    refPSF = imresize(IQ_coherent_sum, [size(IQ_coherent_sum, 1) * rfnZ, size(IQ_coherent_sum, 2) * rfnX], 'bilinear');
%     PSF = refPSF(198:205, 125:132);
%     PSF = refPSF(160:240, 100:156);
    PSF = refPSF(190:210, 118:138);
% end

% PSFn = PSF/max(PSF, [], 'all') .* 256 .* 5;
% x = normxcorr2(abs(PSF), abs(refIQs));
%%
% XC = zeros(size(refIQs));
XC = normxcorr2(abs(PSF), abs(refIQs(:, :, 1))); % Cross correlate the filtered/refined images and the simulated PSF
parfor f = 2:size(refIQs, 3)
% parfor f = 2:10
    XC(:, :, f) = normxcorr2(abs(PSF), abs(refIQs(:, :, f)));
end

%% Make video
% vo = VideoWriter('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processing\diff_50_5000'); % video object
vo = VideoWriter('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processing\XC_50_5000_2'); % video object

% va = abs(d);
va = XC - min(XC, [], 'all');
va = va ./ max(va, [], 'all');

vo.Quality = 100;
vo.FrameRate = 100;
open(vo);
for f = 1:size(va, 3)
    writeVideo(vo, va(:, :, f));
end
close(vo);

%% binary image conversion of all frames
% a = abs(d);
% a = abs(IQs);
% a = abs(refIQs);
a = XC;

% figure; imagesc(a)
% figure; plot(a(:))

% threshold = 7e3;
% threshold = 800;
threshold = 0.4; % XC threshold, set it manually

tic
mask = a > threshold;
bi = zeros(size(a)); bi(mask) = 1; % binary image with white above the threshold
% binaryImage = figure; imagesc(bi); colormap gray
toc
%% Make video of binary image
vo = VideoWriter('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processing\biXC_50_5000_2'); % video object
va = bi;
va = va ./ max(va);

vo.Quality = 100;
vo.FrameRate = 100;
open(vo);
for f = 1:size(va, 3)
    writeVideo(vo, va(:, :, f));
end
close(vo);
%% Centroid finding
nf = size(bi, 3);        % # frames in the binary image stack
centroids = cell(nf, 1); % initialize
volumeThreshold = 1;

parfor f = 1:nf
    cf = bi(:, :, f); % current frame
    CC = bwconncomp(cf); % connected components (connected regions of 1 in a binary image)

    s = regionprops(CC, cf, 'Area', 'WeightedCentroid'); % Get the weighted centroids and area of each connected component

    % Remove connected regions without enough pixels (probably noise)
    for si = numel(s):-1:1
        if s(si).Area <= volumeThreshold
            s = s(1:si - 1);
        end
    end

    centroidsCurrentFrame = zeros(numel(s), 2); % initialize centroid array. Dimensions: # centroids x 2 (z location, x location)

    % go through the s structure and get the .WeightedCentroid data
    for cn = 1:numel(s)
        centroidsCurrentFrame(cn, :) = s(cn).WeightedCentroid;
    end
    
    centroids{f} = centroidsCurrentFrame; % put back into overall variable

end

%% attempt to graph centroids
% centroidPlotFigure = figure;
% hold on
% for f = 1:nf
%     scatter(centroids{f}(:, 1), centroids{f}(:, 2), '.')
% end
% hold off
% 
zpts = [];
xpts = [];
for f = 1:nf
    zpts = [zpts; centroids{f}(:, 1)];
    xpts = [xpts; centroids{f}(:, 2)];
end

hPixFactor = 10; % increase the pixel count by this factor in each dimension
figure;
h = histogram2(zpts, xpts, [size(a, 1) * hPixFactor, size(a, 2) * hPixFactor], 'DisplayStyle','tile');
grid off
colormap hot

%%


% %% binary image test
% a = abs(d(:, :, 2));
% % figure; imagesc(a)
% % figure; plot(a(:))
% 
% threshold = 1e4;
% 
% mask = a > threshold;
% bi = zeros(size(a)); bi(mask) = 1; % binary image with white above the threshold
% binaryImage = figure; imagesc(bi); colormap gray
% 
% % figure; hold on
% % plot(a(:), 'o')
% % plot(bi(:), '^')
% % hold off
% 
% %%
% CC = bwconncomp(bi)
% s = regionprops(CC, bi, 'WeightedCentroid')
% centroids = zeros(numel(s), 2);
% for cn = 1:numel(s)
%     centroids(cn, :) = s(cn).WeightedCentroid;
% end
% 
% radii = ones(size(centroids, 1), 1);
% hold on
% viscircles(centroids, radii)
% hold off
% % s2 = regionprops(bi,'Area','WeightedCentroid')

%%
taustep = P.SeqControl(1).argument ./ 1e6 .* 1e3 .* P_new.na .* (0:P_new.numSubFrames - 1); % in ms

g1 = g1T_1D(IQ_f);

% pt = [40, 44, 63]; % y, x, z
% pt = [44, 40, 63];
pt = [97, 66];

figure; plot(taustep, abs(squeeze(g1(pt(1), pt(2), :))), '-o')
xlabel('tau (ms)')
ylabel('|g1|')
title(strcat("|g1| at ", num2str(pt), " pixel")) % (z, x)

mag_g1 = abs(g1);

%% Filtered Power Doppler for CBV comparison
IQ_cs_f_sq = IQ_f.^2;
I_f_PowerDoppler = abs(squeeze(sum(IQ_cs_f_sq, 3)));
figure; imagesc(I_f_PowerDoppler)
title('Filtered Power Doppler - xz plane')
xlabel('x pixels')
ylabel('z pixels')
%% Jianbo Power Doppler
[PDI]=sIQ2PDI(IQ_f);

figure; imagesc(PDI(:, :, 1))
title('Power Doppler - positive frequencies')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(PDI(:, :, 2))
title('Power Doppler - negative frequencies')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(PDI(:, :, 3))
title('Power Doppler - all frequencies')
xlabel('x pixels')
ylabel('z pixels')

%% same thing but rescaled
dr_scaling = 0.3; % dynamic range scaling factor

% figure; imagesc(PDI(:, :, 1).^dr_scaling)
% title('Power Doppler - positive frequencies')
% xlabel('x pixels')
% ylabel('z pixels')
%
figure; imagesc(PDI(:, :, 2).^dr_scaling)
title('Power Doppler - negative frequencies')
xlabel('x pixels')
ylabel('z pixels')
%
figure; imagesc(PDI(:, :, 3).^dr_scaling)
title('Power Doppler - all frequencies')
xlabel('x pixels')
ylabel('z pixels')
%%
[CBF, CBV] = g1_to_CBi(g1, taustep, 2, 15, 4); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%% CBFi
figure; imagesc(CBF)
title('CBFi - xz plane')
xlabel('x pixels')
ylabel('z pixels')

%% CBVi
figure; imagesc(CBV)
title('CBVi - xz plane')
xlabel('x pixels')
ylabel('z pixels')
%%
%%
dr_scaling = 0.3; % dynamic range scaling factor

figure; imagesc(CBV.^dr_scaling)

title(strcat("CBVi - xz plane with exponential scaling factor = ", num2str(dr_scaling)))
xlabel('x pixels')
ylabel('z pixels')
colormap gray