%% Test regression of out-of-brain principal components

% Notation: in (functional) and out (noise) data
%% Some basic data resizing

data = PDIallBlocks_reg;
refVol = mean(data, 4); % Reference volume: mean over all frames
ds = size(data); % Data size

%% Create a mask in one direction
% figure; imagesc(abs(squeeze(max(refVol, [], 2))') .^ 0.5); colormap gray
% figure; imagesc(abs(squeeze(max(IQfs, [], 2))'))

% % Draw a ROI on the coronal MIP: outside the volume of interest
% figure; imagesc(abs(squeeze(max(refVol, [], 1))') .^ 0.5); colormap gray
% coronal_out_roi = images.roi.Freehand;
% coronal_out_roi.draw;
% coronal_out_mask = createMask(coronal_out_roi);

% Draw a ROI on the sagittal MIP: OUTSIDE the volume of interest
figure; imagesc(abs(squeeze(max(refVol, [], 2))') .^ 0.5); colormap gray
sagittal_roi_out = images.roi.Freehand;
sagittal_roi_out.draw;
sagittal_out_mask = createMask(sagittal_roi_out);
numVoxelsOut_sagittal = length(find(sagittal_out_mask)) * ds(2); % # of voxels in the 3D-expanded sagittal mask


% % Draw a ROI on the sagittal MIP: INSIDE the volume of interest
% figure; imagesc(abs(squeeze(max(refVol, [], 2))') .^ 0.5)
% sagittal_roi_in = images.roi.Freehand;
% sagittal_roi_in.draw;
% sagittal_in_mask = createMask(sagittal_roi_in);
% numVoxelsIn_sagittal = length(find(sagittal_in_mask)) * ds(2); % # of voxels in the 3D-expanded sagittal mask

%% Mask the PDI inside and outside of the volume of interest
% coronal_out_mask_rep = repmat(permute(~coronal_out_mask, [3, 2, 1]), ds(1), 1, 1, ds(4));
sagittal_out_mask_rep = repmat(permute(sagittal_out_mask, [2, 3, 1]), 1, ds(2), 1, ds(4));

% sagittal_in_mask_rep = repmat(permute(sagittal_in_mask, [2, 3, 1]), 1, ds(2), 1, ds(4));

%% PCA on the voxels outside the volume of interest
data_out = reshape(data(sagittal_out_mask_rep), [numVoxelsOut_sagittal, ds(4)])'; % Should be in the dimensions [time (observations), space (random variables)]
data_out_zm = data_out - mean(data_out, 1); % Zero mean
% Note: we want the features to be time
S_out = data_out*data_out'; % Covariance matrix of out-of-brain voxel timecourses
% S_out = data_out_zm*data_out_zm'; % Covariance matrix of out-of-brain voxel timecourses

[U, S, V] = svd(S_out);
figure; plot(U(:, 1:5))
% [evecs_out, evals_out, W] = eig(S_out); % Eigendecomposition of the covariance matrix of out-of-brain voxel timecourses
% % figure; plot(evecs_out(:, 1:1))
% figure; plot(W(:, 1:1))
%% Choose which PCs to keep and combine
numPCsKept = 5;
% % PC_basis = U(:, 1:numPCsKept);
% % PC_basis = evecs_out(:, 1:numPCsKept);
% % figure; plot(PC_basis)
% SVs = diag(S);
% noise_to_subtract = U(:, 1:numPCsKept)*SVs(1:numPCsKept);
% figure; plot(noise_to_subtract)

A = U(:, 1:numPCsKept);
X = reshape(data, [ds(1)*ds(2)*ds(3), ds(4)])';
W = (A'*A)\(A'*X);

X_denoised = X - A*W;

data_denoised = reshape(X_denoised, ds);

% Mean out-of-brain voxel timecourse subtraction
data_denoised2 = data - repmat(permute(mean(data_out, 2), [2, 3, 4, 1]), ds(1), ds(2), ds(3), 1);
%% Look at data_denoised
coord = [94, 136, 44]; % Vessel in the brain
figure
yyaxis left
plot(squeeze(data(coord(1), coord(2), coord(3), :)))
yyaxis right
plot(squeeze(data_denoised(coord(1), coord(2), coord(3), :)))
hold on
plot(squeeze(data_denoised2(coord(1), coord(2), coord(3), :)))
hold off
legend('Original', 'Denoised with PCR', 'Subtract mean out-of-brain voxel timecourse')
%% Other test
% tempdata = squeeze(data(48, 180, 44, :)); % tempdata = tempdata - mean(tempdata);
% noise_to_subtract = mean(data_out, 2);
% test = tempdata - noise_to_subtract;
% figure; plot(test)

%% old CCA stuff
[U_out, S_out, V_out] = svd(data_out, 'econ'); % Already sorted in decreasing order
% SVs_out = diag(S_out);

data_in = reshape(data(sagittal_in_mask_rep), [numVoxelsIn_sagittal, ds(4)])';
[U_in, S_in, V_in] = svd(data_in, 'econ'); % Already sorted in decreasing order
% SVs_in = diag(S_in);

data_out_whitened = U_out * V_out';
data_in_whitened = U_in * V_in';

data_whitened_cat = cat(2, data_in_whitened, data_out_whitened); % Concatenate the whitened in (functional) and out (noise) data

[U_w, S_w, V_w] = svd(data_whitened_cat, 'econ'); % Already sorted in decreasing order
% SVs_w = diag(S_w);

noise_inds = 1:10; % TESTING
% X_noise = U_out(noise_inds, :)' * U_w(noise_inds, :) * data_whitened_cat';
X_noise = U_out(:, noise_inds) * U_w(:, noise_inds, :)' * data_whitened_cat;