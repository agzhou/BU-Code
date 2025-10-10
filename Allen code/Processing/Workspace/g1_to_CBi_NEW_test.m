n_CBV = 2;

g1_tau1_cutoff = 0.2;
% g1_tau1_cutoff = 0.1;

% g1_tau1_cutoff = 0.0;
% tau_difference_cutoff = 0.2;

% for filenum = startFile:endFile
% % for filenum = [endFile]
% % for filenum = 1
% %     load([savepath, 'g1-', num2str(filenum)], 'g1') % Load the saved g1 mat files
%     load([savepath, 'fUSdata-', num2str(filenum)], 'g1') % Load the saved g1 mat files

    [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);
%     [g1A_mask] = createg1mask(g1Avg, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);

    [CBFsi_new, CBVi_new] = g1_to_CBi_NEW(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV, n_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)

%     CBFsi(~g1A_mask) = -Inf; % Remove noisy points from the CBFspeed index (in theory)
    CBFsi_new(~g1A_mask) = 0; % Remove noisy points from the CBFspeed index (in theory)

%     save([savepath, 'tlfUSdata-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
% %     save([savepath, 'tlfUSdatatest-', num2str(filenum), '.mat'], 'CBFsi', 'CBVi', '-v7.3', '-nocompression');
%     disp("tl-fUS result for file " + num2str(filenum) + " saved" )

% end
% save([savepath, 'tlfUS_proc_params.mat'], 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff', 'g1A_mask');
% % save([savepath, 'tlfUStest_proc_params.mat'], 'tau1_index_CBV', 'tau1_index_CBF', 'tau2_index_CBF', 'g1_tau1_cutoff');
figure; imagesc(squeeze(max(CBVi_new(:, :, :), [], 1) .^ 0.3)'); colormap hot
figure; imagesc(squeeze(max(CBVi_new(:, :, :), [], 3) .^ 0.3)'); colormap hot
vcmap = colormap_ULM;
figure; imagesc(squeeze(mean(CBFsi_new(:, :, :), 1))'); colormap(vcmap)

%% Helper functions

function [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF)

    g1A_T = {};
    
    g1A_T{1} = abs(g1(:, :, :, 2)) > g1_tau1_cutoff; % First treatment: tau1 is above some cutoff (make sure there is some actual blood signal there)
%     g1A_T{2} = abs(g1(:, :, :, tau1_index_CBF)) > abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > |g1(tau2)| (noise might have the g1 randomly increase with tau, but it should not happen with a voxel where there is a real blood signal)
%     g1A_T{3} = abs(g1(:, :, :, tau1_index_CBF)) > 2 .* abs(g1(:, :, :, tau2_index_CBF)); % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
%     % g1A_T{4} = abs(g1(:, :, :, tau1_index_CBF)) - 1 .* abs(g1(:, :, :, tau2_index_CBF)) > tau_difference_cutoff; % Keep the voxels where |g1(tau1)| > 2 * |g1(tau2)| (same as #2, but more severe)
    
    g1A_mask = true(size(g1A_T{1})); % Mask of voxels to keep for the g1 treatments
    for i = 1:length(g1A_T)
        g1A_mask = and(g1A_mask, g1A_T{i});
    end

end

