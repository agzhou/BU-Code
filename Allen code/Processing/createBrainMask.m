%%
IQfs = sum(IQ, 4); % Sum IQ over frames
% volumeViewer(abs(IQfs))

%% Use percentiles to automatically find a mask
pctl_cutoff = 90;
IQfs_pctl = prctile(abs(IQfs), pctl_cutoff, 'all');

test = true(size(IQfs));
test(abs(IQfs) >= IQfs_pctl) = false;
volumeViewer(test)
figure; imagesc(squeeze(sum(test, 1))')

%% 
test_mf = medfilt3(test);

%% Use gradients to get edge detection
% [Gmag, Gazimuth, Gelevation] = imgradient3(abs(IQfs));
[Gx, Gy, Gz] = imgradientxyz(abs(IQfs), 'central');
%% Slices
figure; imagesc(abs(squeeze(max(IQfs, [], 1))'))
% figure; imagesc(abs(squeeze(max(IQfs, [], 2))'))

% Draw a ROI on the coronal MIP
coronal_roi = images.roi.Freehand;
coronal_roi.draw;
coronal_mask = createMask(coronal_roi);
% coronal_mask = drawfreehand;
%% Mask the IQ
IQm = IQ; % IQ masked
coronal_mask_rep = repmat(permute(~coronal_mask, [3, 2, 1]), size(IQ, 1), 1, 1, size(IQ, 4));
IQm(coronal_mask_rep) = 0;
figure; imagesc(abs(squeeze(max(IQm(:, :, :, 1), [], 1))'))

volumeViewer(abs(IQm(:, :, :, 1)))

%% Test for surface segmentation
temp_initial_mask = ~coronal_mask_rep(:, :, :, 1);

BW = activecontour(abs(IQfs), temp_initial_mask, 300, 'Chan-Vese', SmoothFactor = 0.2, ContractionBias = -0.0);
%% SVD on the IQ

[xp, yp, zp, nf] = size(IQm);
PP = reshape(IQm, [xp*yp*zp, nf]);
tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
[U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
disp('Full SVD done')
toc

plot_FFT_SVs_function(V, P)
%%
SSM = plotSSM(U);

% Test to look at the individual "weighted images"
k_test = 1; % Which column vector to use
ss_wi = reshape(U(:, k_test) * V(:, k_test)', [xp, yp, zp, nf]); % Subspace weighted image
volumeViewer(abs(ss_wi(:, :, :, 1)))
%     figure; imagesc(abs(mean(test, 4)))

%%
[PP, EVs, V_sort] = getSVs2D(IQm);
disp('SVs decomposed')
[IQfm] = applySVs2D(IQm, PP, EVs, V_sort, 18, 120);
disp('SVD filtered images put together')

%%
volumeViewer(abs(IQfm(:, :, :, 1)))
figure; imagesc(squeeze(abs(max(IQfm(:, :, :, 1), [], 1)))')

PDI = sum(abs(IQfm) .^ 2, 4);
volumeViewer(PDI)
figure; imagesc(squeeze(max(PDI, [], 1))')

%% Get CBVi/CBFsi
g1_tau1_cutoff = 0.3;

numg1pts = 20; % Only calculate the first N points
g1 = g1T(IQfm, numg1pts);

[g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);
[CBFsi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
CBFsi(~g1A_mask) = 0; % Remove noisy points from the CBFspeed index (in theory)

%% Plot CBVi/CBFsi
figure; imagesc(squeeze(max(CBVi(30:50, :, :), [], 1) .^ 0.5)'); colormap hot
figure; imagesc(squeeze(max(CBVi(:, :, :), [], 3) .^ 0.5)'); colormap hot
vcmap = colormap_ULM;
figure; imagesc(squeeze(mean(CBFsi(30:50, :, :), 1))'); colormap(vcmap)

volumeViewer(CBVi)






%% %%%%%%%%%%%%%%%%%%%%%%%%% Same but unmasked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xp, yp, zp, nf] = size(IQ);
PP_um = reshape(IQ, [xp*yp*zp, nf]);
tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
[U_um, S_um, V_um] = svd(PP_um, 'econ'); % Already sorted in decreasing order
disp('Full SVD done')
toc

plot_FFT_SVs_function(V_um, P)

% k_test = 1; % Which column vector to use
% ss_wi = reshape(U_um(:, k_test) * V_um(:, k_test)', [xp, yp, zp, nf]); % Subspace weighted image
% volumeViewer(abs(ss_wi(:, :, :, 1)))
%%
[PP_um, EVs_um, V_sort_um] = getSVs2D(IQ);
disp('SVs decomposed')
[IQf] = applySVs2D(IQ, PP_um, EVs_um, V_sort_um, 5, 120);
disp('SVD filtered images put together')

%%
% volumeViewer(abs(IQf(:, :, :, 1)))
% figure; imagesc(squeeze(abs(max(IQf(:, :, :, 1), [], 1)))')

PDI_um = sum(abs(IQf) .^ 2, 4);
volumeViewer(PDI_um)
figure; imagesc(squeeze(max(PDI_um, [], 1))')

%% Get CBVi/CBFsi
g1_tau1_cutoff = 0.3;

numg1pts = 20; % Only calculate the first N points
g1_um = g1T(IQf, numg1pts);
%%
[g1A_um_mask] = createg1mask(g1_um, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);
[CBFsi_um, CBVi_um] = g1_to_CBi(g1_um, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
CBFsi_um(~g1A_um_mask) = 0; % Remove noisy points from the CBFspeed index (in theory)

%% Plot CBVi/CBFsi
figure; imagesc(squeeze(max(CBVi_um(30:50, :, :), [], 1) .^ 0.5)'); colormap hot
figure; imagesc(squeeze(max(CBVi_um(:, :, :), [], 3) .^ 0.5)'); colormap hot
vcmap = colormap_ULM;
figure; imagesc(squeeze(mean(CBFsi_um(30:50, :, :), 1))'); colormap(vcmap)

volumeViewer(CBVi_um)
%%
% figure; histogram(abs(IQfs), 'Normalization', 'counts')

%% Helper functions

function [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF)

    g1A_T = {};
    
    g1A_T{1} = abs(g1(:, :, :, 2)) > g1_tau1_cutoff; % First treatment: tau1 is above some cutoff (make sure there is some actual blood signal there)
    g1A_T{2} = abs(g1(:, :, :, tau1_index_CBF)) > abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > |g1(tau2)| (noise might have the g1 randomly increase with tau, but it should not happen with a voxel where there is a real blood signal)
    g1A_T{3} = abs(g1(:, :, :, tau1_index_CBF)) > 2 .* abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    % g1A_T{4} = abs(g1(:, :, :, tau1_index_CBF)) - 1 .* abs(g1(:, :, :, tau2_index_CBF)) > tau_difference_cutoff; % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    
    g1A_mask = true(size(g1A_T{1})); % Mask of voxels to keep for the g1 treatments
    for i = 1:length(g1A_T)
        g1A_mask = and(g1A_mask, g1A_T{i});
    end

end