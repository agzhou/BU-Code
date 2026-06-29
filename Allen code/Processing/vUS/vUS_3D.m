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
tp = [40, 40, 71]; % Test point
% figure; plot(squeeze(abs(IQf(tp(1), tp(2), tp(3), :))))
% figure; plot(squeeze(real(IQf(tp(1), tp(2), tp(3), :))))
% figure; plot(squeeze(real(IQf_HPF(tp(1), tp(2), tp(3), :))))

%% ========= 2. Directional flow filtering ========= %%
% 2.1 Separate positive and negative frequencies
[IQf_separated, IQf_FT_separated, nFTpts] = separatePosNegFreqs(IQf_HPF);

% Testing
faxis = linspace(-P.frameRate/2, P.frameRate/2, nFTpts)';
% figure; plot(faxis, squeeze(abs(IQf_FT_separated{1}(tp(1), tp(2), tp(3), :))))
% figure; plot(faxis, squeeze(abs(IQf_FT_separated{2}(tp(1), tp(2), tp(3), :))))
% figure; plot(faxis, squeeze(abs(IQf_FT_separated{3}(tp(1), tp(2), tp(3), :))))
figure; plot(1:P.numFramesPerBuffer, squeeze(abs(IQf_separated{1}(tp(1), tp(2), tp(3), :))))
figure; plot(1:P.numFramesPerBuffer, squeeze(abs(IQf_separated{2}(tp(1), tp(2), tp(3), :))))
figure; plot(1:P.numFramesPerBuffer, squeeze(abs(IQf_separated{3}(tp(1), tp(2), tp(3), :))))

% 2.2 Screen voxels for signal power

% 2.3 Screen voxels for noisiness
