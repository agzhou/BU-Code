
%%
% addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Processing\Speckle tracking')
addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Simulations')
%% Load data

%%

IQ = squeeze(IData + 1i .* QData);

%%
%     volumeViewer(squeeze(abs(IQ(:, :, :, 1))))
% figure; imagesc(squeeze(max( abs(IQ(:, :, :, 1)), [], 1) )')
figure; imagesc(squeeze( abs(IQ(40, :, :, 1)) )')
%% Calculate g1T
num_g1_pts = 50;
g1 = sim_g1T(IQ, num_g1_pts);

%%
figure; plot(squeeze(abs(IQ(40, 42, 62, :))))
% figure; plot(squeeze(abs(g1(40, 42, 62, :))))
% figure; plot(squeeze(abs(g1(40, 45, 62, 1:50))))