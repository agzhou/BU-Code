%% Manually define a 3D mask with the intersection of three orthogonal masks
rt = 0.1; % r threshold

% Coronal
PDIA_MIP = squeeze(max(PDIA, [], 1)).';
am_rPDIrs_MIP = squeeze(max(am_rPDIrs_TA, [], 1)).';
% overlayActivationMap2D(fliplr(PDIA_MIP), fliplr(am_rPDIrs_MIP), x_mm, z_mm, 'Mask', fliplr(am_rPDIrs_MIP) > rt, 'Alpha', 0.7, 'AmapClim', [min(am_rPDIrs_TA, [], 'all'), max(am_rPDIrs_TA, [], 'all')], 'ColorbarLabel', 'am', 'Title', 'Activation map over power Doppler');
overlayActivationMap2D(PDIA_MIP, am_rPDIrs_MIP, x_mm, z_mm, 'Mask', am_rPDIrs_MIP > rt, 'Alpha', 0.7, 'AmapClim', [min(am_rPDIrs_TA, [], 'all'), max(am_rPDIrs_TA, [], 'all')], 'ColorbarLabel', 'am', 'Title', 'Activation map over power Doppler');
coronal_mask = roipoly; % manually define the ROI
% figure; imagesc(coronal_mask)
coronal_mask_3D = repmat(permute(coronal_mask, [3, 2, 1]), xp, 1, 1);

% Sagittal
PDIA_MIP = squeeze(max(PDIA, [], 2)).';
am_rPDIrs_MIP = squeeze(max(am_rPDIrs_TA, [], 2)).';
% overlayActivationMap2D(fliplr(PDIA_MIP), fliplr(am_rPDIrs_MIP), x_mm, z_mm, 'Mask', fliplr(am_rPDIrs_MIP) > rt, 'Alpha', 0.7, 'AmapClim', [min(am_rPDIrs_TA, [], 'all'), max(am_rPDIrs_TA, [], 'all')], 'ColorbarLabel', 'am', 'Title', 'Activation map over power Doppler');
overlayActivationMap2D(PDIA_MIP, am_rPDIrs_MIP, x_mm, z_mm, 'Mask', am_rPDIrs_MIP > rt, 'Alpha', 0.7, 'AmapClim', [min(am_rPDIrs_TA, [], 'all'), max(am_rPDIrs_TA, [], 'all')], 'ColorbarLabel', 'am', 'Title', 'Activation map over power Doppler');
sagittal_mask = roipoly; % manually define the ROI
% figure; imagesc(sagittal_mask)
sagittal_mask_3D = repmat(permute(sagittal_mask, [2, 3, 1]), 1, yp, 1);

% Axial
PDIA_MIP = squeeze(max(PDIA, [], 3)).';
am_rPDIrs_MIP = squeeze(max(am_rPDIrs_TA, [], 3)).';
% overlayActivationMap2D(fliplr(PDIA_MIP), fliplr(am_rPDIrs_MIP), x_mm, z_mm, 'Mask', fliplr(am_rPDIrs_MIP) > rt, 'Alpha', 0.7, 'AmapClim', [min(am_rPDIrs_TA, [], 'all'), max(am_rPDIrs_TA, [], 'all')], 'ColorbarLabel', 'am', 'Title', 'Activation map over power Doppler');
overlayActivationMap2D(PDIA_MIP, am_rPDIrs_MIP, x_mm, z_mm, 'Mask', am_rPDIrs_MIP > rt, 'Alpha', 0.7, 'AmapClim', [min(am_rPDIrs_TA, [], 'all'), max(am_rPDIrs_TA, [], 'all')], 'ColorbarLabel', 'am', 'Title', 'Activation map over power Doppler');
axial_mask = roipoly; % manually define the ROI
% figure; imagesc(axial_mask)
axial_mask_3D = repmat(permute(axial_mask, [2, 1, 3]), 1, 1, zp);

mask_3D = logical(coronal_mask_3D .* sagittal_mask_3D .* axial_mask_3D);
% volumeViewer(mask_3D)

%% Get timecourse
numPtsUSI = size(PDIrs_TA, fDim);
% Calculate the timecourse from the average within that ROI
roi_rPDIrs_TA = zeros(size(rPDIrs_TA, 3), 1);
% repmat(roi_mask, [1, 1, stim_pattern.trial_duration])
for ti = 1:numPtsUSI
    % for ti = 1
    temp_rPDIrs_TA = rPDIrs_TA(:, :, :, ti);
    roi_rPDIrs_TA(ti) = mean(temp_rPDIrs_TA(mask_3D));
end

% Plot the average timecourse in the ROI
rsT = (0:numPtsUSI-1) ./ rr; % Resampled time vector [s]
figure; plot(rsT, roi_rPDIrs_TA, 'LineWidth', 2); xlabel('Time [s]'); ylabel('rPDI'); title('ROI and trial-averaged rPDI')
hold on
% plot(rsT, stimrs1T)

maxValue = max(roi_rPDIrs_TA) + 0.01; % max value of the stim, for shading
minValue = min(roi_rPDIrs_TA) - 0.01;
yshade = [maxValue, minValue, minValue, maxValue]; % Go top right and clockwise for the shading patch vertices
xshade = [((P.apis.delay_time_ms/1e3 + stim_delay) + P.apis.stim_length_s), ((P.apis.delay_time_ms/1e3 + stim_delay) + P.apis.stim_length_s), (P.apis.delay_time_ms/1e3 + stim_delay), (P.apis.delay_time_ms/1e3 + stim_delay)];
patch(xshade, yshade, 'g', 'FaceAlpha', .1, 'EdgeAlpha', 0) % Plot the shaded region
ylim([minValue, maxValue])
legend("rPDI", "Stim on")
clearvars xshade yshade maxValue minValue