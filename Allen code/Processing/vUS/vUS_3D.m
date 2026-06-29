%% Description:
%   Calculate flow velocity from volumetric fUS data using the vUS method
%   (Tang et al., 2020)
% Inputs:
%   IQ: [x voxels, y voxels, z voxels, frames] complex data matrix
% Outputs:
%   vUS: [x voxels, y voxels, z voxels, 3 (xyz components)] real velocity data matrix

%%
% function [vUS] = vUS_3D(IQ)

%% Set up the High Pass Filter (parameters from the 2020 vUS paper)
HPF.fc = 25; % Cutoff frequency [Hz]
% 25 Hz corresponds to 1 mm/s

HPF.fs = P.frameRate; % Sampling frequency [Hz]
HPF.order = 4; % Butterworth filter order

[HPF.b, HPF.a] = butter(HPF.order, HPF.fc/(HPF.fs/2), 'high');

%% ========= 1. Preprocessing ========= %%

% 1.1 SVD clutter filter
%     [PP, EVs, V_sort] = getSVs2D(IQ);
[xp, yp, zp, nf] = size(IQ);
PP = reshape(IQ, [xp*yp*zp, nf]);
tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
[U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
SVs = diag(S);
%     disp('Full SVD done')
toc
disp('SVs decomposed')

[IQf, noise] = applySVs2D(IQ, PP, SVs, V, sv_threshold_lower, sv_threshold_upper);

% 1.2 High pass filter (apply to the post-SVD clutter filtered data)
HPF.dim = length(size(IQf)); % Operate on the time dimension
IQf_HPF = filter(HPF.b, HPF.a, IQf, [], HPF.dim);

% Testing
% tp = [40, 40, 71]; % Test point
tp = [40, 43, 87]; % Test point
figure; plot(squeeze(abs(IQf(tp(1), tp(2), tp(3), :))))
figure; plot(squeeze(real(IQf(tp(1), tp(2), tp(3), :))))
figure; plot(squeeze(real(IQf_HPF(tp(1), tp(2), tp(3), :))))

%% ========= 2. Directional flow filtering ========= %%

% 2.1 Separate positive and negative frequencies
[IQf_separated, IQf_FT_separated, nFTpts] = separatePosNegFreqs(IQf_HPF); % Outputs are cell arrays in the order of: negative, positive, all frequencies
frameDim = length(size(IQf)); % Get the dimension corresponding to time/frames

% Testing: plot the separated and full Fourier spectrums and reconstructed IQ signals
faxis = linspace(-P.frameRate/2, P.frameRate/2, nFTpts)';
figure; plot(faxis, squeeze(abs(IQf_FT_separated{1}(tp(1), tp(2), tp(3), :))))
figure; plot(faxis, squeeze(abs(IQf_FT_separated{2}(tp(1), tp(2), tp(3), :))))
figure; plot(faxis, squeeze(abs(IQf_FT_separated{3}(tp(1), tp(2), tp(3), :))))
% figure; plot(1:P.numFramesPerBuffer, squeeze(abs(IQf_separated{1}(tp(1), tp(2), tp(3), :))))
% figure; plot(1:P.numFramesPerBuffer, squeeze(abs(IQf_separated{2}(tp(1), tp(2), tp(3), :))))
% figure; plot(1:P.numFramesPerBuffer, squeeze(abs(IQf_separated{3}(tp(1), tp(2), tp(3), :))))

% 2.2 Screen voxels for signal power
IQf_FT_power = {sum(abs(IQf_FT_separated{1}), frameDim), sum(abs(IQf_FT_separated{2}), frameDim), sum(abs(IQf_FT_separated{3}), frameDim)}';
% Get masks for signal quality (Eqs. 13-14 in the vUS paper)
Rneg = IQf_FT_power{1} ./ IQf_FT_power{3};
Rpos = IQf_FT_power{2} ./ IQf_FT_power{3};

% % Apply a median filter to the Rneg and Rpos volumes before screening with a threshold
% % mf_kernel_size = [3, 3, 3];
% % mf_kernel_size = [5, 5, 5];
% mf_kernel_size = [9, 9, 9];
% Rneg_mf = medfilt3(Rneg, mf_kernel_size);
% Rpos_mf = medfilt3(Rpos, mf_kernel_size);

% % Get masks for voxels to keep, according to the Rneg and Rpos volumes
% R_threshold = 0.4;
% Rneg_mask = Rneg > R_threshold;
% Rpos_mask = Rpos > R_threshold;

% Testing/visualization
% figure; imagesc(squeeze(max(sum(abs(IQf_separated{3}), 4), [], 1))')
% volumeViewer(Rneg)
% volumeViewer(Rpos)
figure; imagesc(squeeze(max(Rneg, [], 1))'); colormap gray
figure; imagesc(squeeze(max(Rpos, [], 1))'); colormap gray
% volumeViewer(Rneg_mask)
% volumeViewer(Rpos_mask)
% figure; imagesc(squeeze(max(Rneg_mask, [], 1))'); colormap gray
% figure; imagesc(squeeze(max(Rpos_mask, [], 1))'); colormap gray
% figure; imagesc(squeeze(max(Rneg_mf, [], 1))'); colormap gray
% figure; imagesc(squeeze(max(Rpos_mf, [], 1))'); colormap gray

%% ========= 3. Calculate g1 ========= %%
nTau = ceil(10e-3 *P.frameRate); % # of time lags to consider; empirically set by assuming all g1 for voxels containing actual flow decay within 10 ms
g1neg = g1T(IQf_separated{1}, nTau);
g1pos = g1T(IQf_separated{2}, nTau);

%% ========= 4. Clean data ========= %%

% 4.1 Screen voxels for noisiness, through |g1(tau1)|
g1_tau1_threshold = 0.2;
g1neg_tau1_mask = abs(squeeze(g1neg(:, :, :, 2))) > g1_tau1_threshold; % Use index 2 because index 1 corresponds to tau = 0
g1pos_tau1_mask = abs(squeeze(g1pos(:, :, :, 2))) > g1_tau1_threshold;

% Testing/visualization
figure; plot(squeeze(abs(g1neg(tp(1), tp(2), tp(3), :))), '-o'); title('Negative frequencies')
figure; plot(squeeze(abs(g1pos(tp(1), tp(2), tp(3), :))), '-o'); title('Positive frequencies')
volumeViewer(g1neg_tau1_mask)
volumeViewer(g1pos_tau1_mask)

% 4.2 Apply mask
% ...

%% ========= 5. Fit vUS ========= %%
