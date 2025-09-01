% filenum = 41
% 
% %% Define the mask manually for now
% 
% % % load('E:\Allen BME-BOAS-27 Data Backup\AZ01 fUS\07-21-2025 awake RC15gV manual right whisker stim\coronal_mask_rep_07_24_2025.mat')
% % % load('I:\Ultrasound Data from 04-11-2025 to 05-08-2025\05-06-2025 AZ03 fUS pre-stroke\run 1 all frames stacked\coronal_mask_rep_07_31_2025.mat')
% % load('J:\Ultrasound data from 7-21-2025\08-06-2025 AZ01 RCA fUS\coronal_mask_rep.mat')
% 
% 
% %% Load the IQ data
% tic
% load([IQpath, IQfilenameStructure, num2str(filenum)])
% 
% IQ = single(squeeze(IData + 1i .* QData));
% clearvars IData QData
% 
% % figure; imagesc(squeeze(max(abs(IQ(:, :, :, 2)), [], 1))')
% 
% % Crop the IQ first 
%     zstart = 40;
% %     zstart = 50;
%     zend = size(IQ, 3);
%     IQm = IQ(:, :, zstart:zend, :);
% %     figure; imagesc(squeeze(max(abs(IQm(:, :, :, 2)), [], 1))')
% 
% %%%%%%%%%%%%%% IF USING THE PREDEFINED MASK %%%%%%%%%%%%
% % IQm = IQ;
% % IQm(coronal_mask_rep) = 0; % Apply the brain mask to the IQ: set the non-brain voxels equal to 0
% 
% % Apply the HPF
% % dim = length(size(IQm)); % Operate on the time dimension
% % IQm_HPF = filter(HPF_b, HPF_a, IQm, [], dim);
% 
% % SVD decluttering
% %     [PP, EVs, V_sort] = getSVs2D(IQm);
% [xp, yp, zp, nf] = size(IQm);
% PP = reshape(IQm, [xp*yp*zp, nf]);
% tic
% %     [U, S, V] = svd(PP); % Already sorted in decreasing order
% [U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
% SVs = diag(S);
% %     disp('Full SVD done')
% toc
% disp('SVs decomposed')
% 
% [IQf, noise] = applySVs2D(IQm, PP, SVs, V, sv_threshold_lower, sv_threshold_upper);
% %     [IQf, noise] = applySVs2D(IQm, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
% disp('SVD filtered images put together')
% 
% %     volumeViewer(abs(IQf(:, :, :, 1)))
% %     figure; imagesc(squeeze(abs(max(IQf(:, :, :, 1), [], 1)))')
% % clearvars IQ
% 
% % Use the IQf with separated negative and positive frequency components
% %     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);

%%
numg1pts = 20; % Only calculate the first N points
%     g1_n = g1T(IQf_separated{1}, numg1pts);
% %     [CBFsi_n, CBVi_n] = g1_to_CBi(g1_n, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%     g1_p = g1T(IQf_separated{2}, numg1pts);
%     [CBFsi_p, CBVi_p] = g1_to_CBi(g1_p, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
% 
% g1 = g1T(IQf, numg1pts);
g1 = g1T_test(IQf, numg1pts);