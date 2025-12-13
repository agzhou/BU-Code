


%% load params and stuff
IQpath = uigetdir('D:\Allen\Data\', 'Select the IQ data path');
IQpath = [IQpath, '\'];

% Load parameters
% if ~exist('P', 'var')
%     load([IQpath, '..\params.mat'])
% end
% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [IQpath, '..\params.mat']);
    load([params_pathname, params_filename])
end

% Load Verasonics reconstruction parameters: datapath\PData.mat
if ~exist('PData', 'var')
    load([IQpath, 'PData.mat'])
end

IQfilenameStructure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

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

%% Calculate time lags
taustep = 1/P.frameRate;
% tau = taustep:taustep:(P.numFramesPerBuffer * taustep);
tau = 0:taustep:((P.numFramesPerBuffer - 1) * taustep);
tau_ms = tau .* 1000; % Assuming even time spacing between frames

%%
% figure; plot(squeeze(abs(IQ(40, 42, 62, :))))
figure; plot(tau_ms(1:num_g1_pts), squeeze(abs(g1(40, 40, 60, :)))); xlabel("Tau [ms]"); ylabel("|g_1|")
% figure; plot(squeeze(abs(g1(40, 45, 62, 1:50))))