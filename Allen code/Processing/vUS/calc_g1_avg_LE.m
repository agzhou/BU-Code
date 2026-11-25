%% Choose processed data directory, params
if ~exist('PDpath', 'var')
    PDpath = uigetdir('D:\Allen\Data\', 'Select the Processed data path');
    PDpath = [PDpath, '\'];
end

% Load processing parameters: fUS_proc_params.mat
load([PDpath, 'fUS_proc_params.mat'])

% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [PDpath, '..\params.mat']);
    load([params_pathname, params_filename])
end

%% Define blocks to use
blocksToUse = 1:numBlocks;
blockFilenames = {};
if isfield(P, 'fsFileName') % If FastSave was used, indices start from 0
    fileOffset = 0;
else
    fileOffset = 1;
end
for bi = 1:length(blocksToUse)
    blockNum = blocksToUse(bi);

    blockStartFile = (blockNum-1)*blockSize + fileOffset;
    blockEndFile = blockStartFile + blockSize - 1;
    blockFilenames{bi} = sprintf('fUSdata-block-%d-%d.mat', blockStartFile, blockEndFile);
end

%%
load([PDpath, blockFilenames{1}], 'g1') % load first file
% voxelTimeseriesGUI(squeeze(g1(:, :, :, :, 2)), abs(g1(:, :, :, 2, 2)), 'Colormap', 'turbo', 'ComplexMode', 'abs')

sliceDim = length(size(g1)); % Slice dimension should be the last
slicesPerBlock = size(g1, sliceDim); % Number of slices per block
fDim = length(size(g1)) - 1; % Frame dimension should be the second to last

% Initialize g1 sum variable
temp = size(g1);
g1_sum = zeros(temp(1:fDim));
for bi = 1:length(blocksToUse)

    disp("Loading block " + num2str(blocksToUse(bi)))
    load([PDpath, blockFilenames{bi}], 'g1')
    g1_sum = g1_sum + sum(g1, sliceDim);
end

g1_avg = g1_sum ./ (slicesPerBlock*length(blocksToUse));
% voxelTimeseriesGUI(g1_avg, abs(g1_avg(:, :, :, 2)), 'Colormap', 'turbo', 'ComplexMode', 'abs')

save([PDpath, 'g1_avg.mat'], 'g1_avg', '-v7.3')
% 
% g1_sum = cell([size(g1)]);
% for j = 1:3
%     g1_sum{j} = zeros(size(g1{j}));
% end
% 
% % for filenum = startFile + 1:endFile
% for filenum = startFile:endFile
%     disp(filenum)
%     tic
%     load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1')
% 
% 
%     for j = 1:3
%         g1_sum{j} = g1_sum{j} + g1{j};
%     end
%     toc
% end
% 
% g1_avg = cell(size(g1_sum));
% for j = 1:3
%     g1_avg{j} = g1_sum{j}./(endFile - startFile + 1);
% end
% 
% save([savepath, 'g1_avg.mat'], 'g1_avg', '-v7.3')
% 
% pixelTimeseriesGUI(g1_avg{3}, abs(g1_avg{3}(:, :, 2)), 'ComplexMode', 'abs', 'Colormap', 'turbo')