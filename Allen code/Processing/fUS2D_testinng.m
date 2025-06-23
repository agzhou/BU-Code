IQ = squeeze(IData + 1i .* QData);
figure; imagesc(squeeze(abs(IQ(:, :, 1))) .^ 0.5); colormap hot


taustep = 1/P.frameRate;
% tau = taustep:taustep:(P.numFramesPerBuffer * taustep);
tau = 0:taustep:((P.numFramesPerBuffer - 1) * taustep);
tau_ms = tau .* 1000; % Assuming even time spacing between frames

tau1_index_CBF = 2;
tau2_index_CBF = 6;
tau1_index_CBV = 2;

%%
% sv_threshold_lower = 10;
% sv_threshold_upper = 32;

sv_threshold_lower = 90;
% sv_threshold_upper = 200;
sv_threshold_upper = 500;


% clearvars IData QData

[PP, EVs, V_sort] = getSVs1D(IQ);
disp('SVs decomposed')

% plot_FFT_SVs_function(V_sort, P)

[IQf] = applySVs1D(IQ, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
disp('SVD filtered images put together')

figure; imagesc(squeeze(abs(IQf(:, :, 1))) .^ 0.5); colormap hot

numg1pts = 10; % Only calculate the first N points
g1 = g1T(IQf, numg1pts);

[CBFsi, CBVi] = g1_to_CBi(g1, tau_ms, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
figure; imagesc(squeeze(CBVi .^ 0.5)); colormap hot

PDI_test = sum(abs(IQf) .^ 2, 3);
figure; imagesc(squeeze(PDI_test .^ 0.5)); colormap hot