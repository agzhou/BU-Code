%% Test the CCA denoising method

% Notation: in (functional) and out (noise) data
%% Create a mask in one direction

data = PDIallBlocks_reg;
refVol = mean(data, 4); % Reference volume
ds = size(data); % Data size

figure; imagesc(abs(squeeze(max(refVol, [], 2))') .^ 0.5)
% figure; imagesc(abs(squeeze(max(IQfs, [], 2))'))

% Draw a ROI on the coronal MIP: outside the volume of interest
figure; imagesc(abs(squeeze(max(refVol, [], 1))') .^ 0.5)
coronal_out_roi = images.roi.Freehand;
coronal_out_roi.draw;
coronal_out_mask = createMask(coronal_out_roi);

% Draw a ROI on the sagittal MIP: OUTSIDE the volume of interest
figure; imagesc(abs(squeeze(max(refVol, [], 2))') .^ 0.5)
sagittal_roi_out = images.roi.Freehand;
sagittal_roi_out.draw;
sagittal_out_mask = createMask(sagittal_roi_out);
numVoxelsOut_sagittal = length(find(sagittal_out_mask)) * ds(2); % # of voxels in the 3D-expanded sagittal mask


% Draw a ROI on the sagittal MIP: INSIDE the volume of interest
figure; imagesc(abs(squeeze(max(refVol, [], 2))') .^ 0.5)
sagittal_roi_in = images.roi.Freehand;
sagittal_roi_in.draw;
sagittal_in_mask = createMask(sagittal_roi_in);
numVoxelsIn_sagittal = length(find(sagittal_in_mask)) * ds(2); % # of voxels in the 3D-expanded sagittal mask

%% Mask the PDI inside and outside of the volume of interest
% coronal_out_mask_rep = repmat(permute(~coronal_out_mask, [3, 2, 1]), ds(1), 1, 1, ds(4));
sagittal_out_mask_rep = repmat(permute(sagittal_out_mask, [2, 3, 1]), 1, ds(2), 1, ds(4));

sagittal_in_mask_rep = repmat(permute(sagittal_in_mask, [2, 3, 1]), 1, ds(2), 1, ds(4));

%% SVD on the voxels outside and inside the volume of interest
data_out = reshape(data(sagittal_out_mask_rep), [numVoxelsOut_sagittal, ds(4)])'; % Should be in the dimensions [time (observations), space (random variables)]
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