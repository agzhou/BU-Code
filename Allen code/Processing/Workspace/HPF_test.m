%% Set up the filter
fc = 50; % Cutoff frequency [Hz]
fs = 1000; % Sampling frequency [Hz]
HPF_order = 3; % Butterworth filter order

[b, a] = butter(HPF_order, fc/(fs/2), 'high');

figure
freqz(b, a, [], fs)

subplot(2, 1, 1)
ylim([-100, 20])

%% Apply the filter to IQ data
dim = length(size(IQ)); % Operate on the time dimension
IQ_HPF = filter(b, a, IQ, [], dim);

figure; imagesc(squeeze(max(abs(IQ_HPF(:, :, :, 1)), [], 1))')

%% Redo SVD
IQm = IQ(:, :, 40:end, :);
figure; imagesc(squeeze(max(abs(IQm(:, :, :, 2)), [], 1))')

IQm_HPF = IQ_HPF(:, :, 40:end, :);
figure; imagesc(squeeze(max(abs(IQm_HPF(:, :, :, 2)), [], 1))')

    %%%%%%%%%%%%%% IF USING THE MASK %%%%%%%%%%%%
%     IQm(coronal_mask_rep) = 0; % Apply the brain mask to the IQ: set the non-brain voxels equal to 0
    
    % Determine the optimal SV thresholds with the spatial similarity matrix
    [xp, yp, zp, nf] = size(IQm_HPF);
    PP = reshape(IQm_HPF, [xp*yp*zp, nf]);
    tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
    [U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
    SVs = diag(S);
%     disp('Full SVD done')
    toc

    SSM = plotSSM(U, false);
%     SSM = plotSSM(U, true);
    [~, a_opt, b_opt] = fitSSM(SSM, false); % Get the optimal singular value thresholds
%     [~, a_opt, b_opt] = fitSSM(SSM, true); % Get the optimal singular value thresholds
    

%     [PP, EVs, V_sort] = getSVs2D(IQ);
%     disp('SVs decomposed')
%     [IQf_HPF, noise] = applySVs2D(IQm_HPF, PP, SVs, V, a_opt, b_opt);
    [IQf_HPF, noise] = applySVs2D(IQm_HPF, PP, SVs, V, sv_threshold_lower, sv_threshold_upper);
%     disp('SVD filtered images put together')

%     volumeViewer(abs(IQf_HPF(:, :, :, 1)))
%     figure; imagesc(squeeze(abs(max(IQf_HPF(:, :, :, 1), [], 1)))')
    % clearvars IQ

    % Use the IQf with separated negative and positive frequency components
%     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);
    
    numg1pts = 20; % Only calculate the first N points
%     g1_n = g1T(IQf_separated{1}, numg1pts);
% %     [CBFsi_n, CBVi_n] = g1_to_CBi(g1_n, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%     g1_p = g1T(IQf_separated{2}, numg1pts);
%     [CBFsi_p, CBVi_p] = g1_to_CBi(g1_p, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
    g1 = g1T(IQf_HPF, numg1pts);

%     [PDI] = calcPowerDoppler(IQf_separated);
    PDI = sum(abs(IQf_HPF) .^ 2, 4) ./ size(IQf_HPF, 4);
%     [CDI] = calcColorDoppler(IQf_FT_separated, P);

%     figure; imagesc(squeeze(max(PDI, [], 1))' .^ 0.5); colormap hot
%     figure; imagesc(squeeze(max(PDI ./ noise, [], 1))' .^ 0.5); colormap hot
%     volumeViewer(PDI)

%% tl-fUS
g1_tau1_cutoff = 0.2;
% tau_difference_cutoff = 0.2;

    [g1A_mask] = createg1mask(g1, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);

    [CBFsi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)

    CBFsi(~g1A_mask) = 0; % Remove noisy points from the CBFspeed index (in theory)

figure; imagesc(squeeze(max(CBVi(:, :, :), [], 1) .^ 0.5)'); colormap hot
figure; imagesc(squeeze(max(CBVi(:, :, :), [], 3) .^ 0.5)'); colormap hot
vcmap = colormap_ULM;
figure; imagesc(squeeze(mean(CBFsi(:, :, :), 1))'); colormap(vcmap)

%% tl-fUS
g1_tau1_cutoff = 0.2;
% tau_difference_cutoff = 0.2;

    [g1A_mask2] = createg1mask(g12, g1_tau1_cutoff, tau1_index_CBF, tau2_index_CBF);

    [CBFsi2, CBVi2] = g1_to_CBi(g12, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)

    CBFsi2(~g1A_mask2) = 0; % Remove noisy points from the CBFspeed index (in theory)

figure; imagesc(squeeze(max(CBVi2(:, :, :), [], 1) .^ 0.5)'); colormap hot
figure; imagesc(squeeze(max(CBVi2(:, :, :), [], 3) .^ 0.5)'); colormap hot
vcmap = colormap_ULM;
figure; imagesc(squeeze(mean(CBFsi2(:, :, :), 1))'); colormap(vcmap)

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