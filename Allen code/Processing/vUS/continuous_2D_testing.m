%% Description: Run some tests on chunks of continuously acquired 2D fUS data

%% Load data and things
if ~exist('IQpath', 'var')
    IQpath = uigetdir('D:\Allen\Data\', 'Select the IQ data path');
    IQpath = [IQpath, '\'];
end

% Load acquisition parameters: params.mat
if ~exist('P', 'var')
    % Choose and load the params.mat file (from the acquisition)
    [params_filename, params_pathname, ~] = uigetfile('*.mat', 'Select the params file', [IQpath, '..\params.mat']);
    load([params_pathname, params_filename])
end

% load some chunk of continuous superframes
IQfilenameStructure = ['IQ-', num2str(P.maxAngle), '-', num2str(P.na), '-', num2str(P.frameRate), '-', num2str(P.numFramesPerBuffer), '-1-'];

startFile = 1; endFile = 4;
nsfpe = (endFile - startFile + 1)*P.numFramesPerBuffer;% Number of superframes per ensemble
load([IQpath, IQfilenameStructure, num2str(startFile)])
IQle = []; % long ensemble IQ

for fi = startFile:endFile
    load([IQpath, IQfilenameStructure, num2str(fi)])
    IQle = cat(3, IQle, IQ);
end

IQ = IQle; clearvars IQle


faxis = linspace(-P.frameRate/2, P.frameRate/2, nsfpe);
%% Pixelwise IQ average
IQ_raw_pixel_avg = squeeze(mean(IQ, [1, 2]));

% figure; plot(abs(IQ_raw_pixel_avg))
figure; plot(movmean(abs(IQ_raw_pixel_avg), 20)); title("|IQ pixelwise average|")

IQ_raw_pixel_avg_FT = fftshift(fft(IQ_raw_pixel_avg));
figure; plot(faxis, abs(IQ_raw_pixel_avg_FT)); xlabel('Frequency [Hz]'); ylabel('|FFT(raw IQ)|')
%% SVD
%% ========= 1. Preprocessing ========= %%


% Mask the region to actually process
[zpo, xpo, nfo] = size(IQ); % Original sizes
% figure; imagesc(squeeze(abs(IQ(:, :, 1))))
% zrange = 1:100;
% zrange = 1:zpo;
zrange = 16:136;
xrange = 1:xpo;
% xrange = 40:60;
IQm = IQ(zrange, xrange, :);

% 1.1 SVD clutter filter
%     [PP, EVs, V_sort] = getSVs2D(IQ);
[zp, xp, nf] = size(IQm);

CM = reshape(IQm, [zp*xp, nf]); % Covariance matrix
tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
[U, S, V] = svd(CM, 'econ'); % Already sorted in decreasing order
SVs = diag(S);
%     disp('Full SVD done')
toc
disp('SVs decomposed')

SSM = plotSSM(U, true);

% sv_threshold_lower = 30; sv_threshold_upper = size(IQ, 3);


[IQf, noise] = applySVs1D(IQm, CM, SVs, V, sv_threshold_lower, sv_threshold_upper);

% 1.2 High pass filter (apply to the post-SVD clutter filtered data)
HPF.dim = length(size(IQf)); % Operate on the time dimension
IQf_HPF = filter(HPF.b, HPF.a, IQf, [], HPF.dim);

% Testing
% figure; imagesc(squeeze(abs(IQf(:, :, 1))))
PDI = sum(abs(IQf).^2, 3);
figure; imagesc(PDI .^ 0.5)

%% FFT of clutter filtered IQ
IQ_cf_pixel_avg = squeeze(mean(IQf, [1, 2]));
figure; plot(abs(IQ_cf_pixel_avg))

IQ_cf_pixel_avg_FT = fftshift(fft(IQ_cf_pixel_avg));
figure; plot(faxis, abs(IQ_cf_pixel_avg_FT)); xlabel('Frequency [Hz]'); ylabel('|FFT(clutter filtered IQ)|')

%% Calculate g1
nTau = ceil(1000e-3 *P.frameRate); % # of time lags to consider; empirically set by assuming all g1 for voxels containing actual flow decay within 10 ms
tau = (0:nTau - 1)' ./ P.frameRate; % Time lag vector [s]

% g1 = g1T(IQf(53, 99, :));
g1 = g1T(IQf, nTau);

%%
pixelTimeseriesGUI(g1, PDI, 'Colormap', 'hot', 'ComplexMode', 'abs')

%%
testpt_g1 = squeeze(g1(51, 79, :));
figure; plot(tau, abs(testpt_g1))

testpt_g1_FT = fftshift(fft(testpt_g1, nsfpe));
figure; plot(faxis, abs(testpt_g1_FT)); xlabel('Frequency [Hz]'); ylabel('|FFT(clutter filtered IQ)|')

%%
% g1_pixel_avg = squeeze(mean(g1, [1, 2]));
% figure; plot(tau, abs(g1_pixel_avg))

g1_FT = fftshift(fft(g1, nsfpe, 3)); % FFT for every pixel
g1_FT_avg = squeeze(mean(g1_FT, [1, 2]));
figure; plot(faxis, abs(g1_FT_avg)); xlabel('Frequency [Hz]'); ylabel('|FFT(clutter filtered IQ)|')
