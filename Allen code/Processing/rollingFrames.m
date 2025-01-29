%% Create rolling re-sampled frames on delay and summed IQ data (pre-compounding)
% for row column array

% Need to make this into a function
%% load test data
load('F:\Allen\Data\12-16-2024 Phantom\RC15gV\run 2 11 angles -5 to 5 deg\IQ.mat')
load('F:\Allen\Data\12-16-2024 Phantom\RC15gV\run 2 11 angles -5 to 5 deg\params.mat')
%% Initialize variabes and separate R-C and C-R volumes
[xp, yp, zp, nacq, nf] = size(IQ);
na = P.na;

IQ_CR = IQ(:, :, :, 1:na);          % column row volumes
IQ_RC = IQ(:, :, :, na + 1:2*na);   % row column volumes

img_FMAS = zeros(xp, yp, zp, nf);        % initiaize final FMAS image