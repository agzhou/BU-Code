%% Test the CCA denoising method


%% Create a mask in one direction

refVol = mean(PDIallBlocks_reg, 4); % Reference volume
ds = size(PDIallBlocks_reg); % Data size
figure; imagesc(abs(squeeze(max(refVol, [], 2))') .^ 0.5)
% figure; imagesc(abs(squeeze(max(IQfs, [], 2))'))

% Draw a ROI on the coronal MIP: outside the volume of interest
figure; imagesc(abs(squeeze(max(refVol, [], 1))') .^ 0.5)
coronal_out_roi = images.roi.Freehand;
coronal_out_roi.draw;
coronal_out_mask = createMask(coronal_out_roi);

% Draw a ROI on the sagittal MIP: outside the volume of interest
figure; imagesc(abs(squeeze(max(refVol, [], 2))') .^ 0.5)
sagittal_roi_out = images.roi.Freehand;
sagittal_roi_out.draw;
sagittal_out_mask = createMask(sagittal_roi_out);
%% Mask the PDI outside of the 
coronal_out_mask_rep = repmat(permute(~coronal_out_mask, [3, 2, 1]), ds(1), 1, 1, ds(4));
sagittal_out_mask_rep = repmat(permute(~sagittal_out_mask, [2, 3, 1]), 1, ds(2), 1, ds(4));