figure; imagesc(squeeze(max(PDI(:, :, :, 3), [], 1)).')
figure; imagesc(squeeze(max(abs(g1(:, :, :, 2, 2)), [], 1)).')
voxelTimeseriesGUI(g1(:, :, :, :, 1), abs(g1(:, :, :, 2, 2)), 'ComplexMode', 'abs')

for si = 1:size(PDI, 4)
    % figure; imagesc(squeeze(max(PDI(:, :, :, si), [], 1)).')
    figure; imagesc(squeeze(max(abs(g1(:, :, :, 2, si)), [], 1)).')
end

tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
[U, S, V] = svd(IQblock, 'econ'); % Already sorted in decreasing order
SVs = diag(S);
%     disp('Full SVD done')
toc