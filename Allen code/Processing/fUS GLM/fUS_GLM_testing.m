%% sdfasdf


tbasis(:,1,iConc) = (exp(1)*(tHRF-tau).^2/sigma^2) .* exp( -(tHRF-tau).^2/sigma^2 );
lstNeg = find(tHRF<0);
tbasis(lstNeg,1,iConc) = 0;

if tHRF(1)<tau
    tbasis(1:round((tau-tHRF(1))/dt),1,iConc) = 0;
end

%% Options for HRF shape
% HRFtype = 'boxcar';
HRFtype = 'gamma';

trange = [-2, 30]; % Time range to define HRF over [s]

% Options for the gamma function
% tau = 1;
% sigma = 1;
a = 3;
b = 2;
% test = -2:0.1:15; HRFtest = gampdf(test, a, b); figure; plot(test, HRFtest)

%% Get some basic parameters and stuff
nT = length(ti.t); % # of samples (time points)

dt = median(diff(ti.t)); % Time step; use mean of all the timesteps in case the sampling is not uniform
fq = 1/dt;        % Sampling/measurement frequency
nPre = round(trange(1)/dt);  % Which index for the block average Pre point
nPost = round(trange(2)/dt); % Which index for the block average Post point
tHRF = (nPre * dt:dt:nPost * dt)'; % Time points for the HRF
ntHRF = length(tHRF); % # of points in the HRF

%% Reshape input data: we want it to be of size [# samples, # voxels (vectorized)]

data = reshape(permute(PDIallSF, [4, 1, 2, 3]), nT, size(PDIallSF, 1)*size(PDIallSF, 2)*size(PDIallSF, 3));


%% Define regressors
regressors = {}; % Cell array of each regressor, which are each a vector

% 1. stim/HRF
switch HRFtype
    case 'boxcar'
        stimHRFs = ti.stimAmps; % Technically this might need to be binarized for a true boxcar
    case 'gamma'
        HRF = gampdf(tHRF, a, b); % figure; plot(tHRF, HRF)
        stimHRFs = conv(ti.stimAmps, HRF, 'same'); % figure; plot(ti.t, stimHRFs)

end
regressors{end + 1} = stimHRFs;

% 2. Drift
driftOrder = 2;
for do = 0:driftOrder % Start at 0 (DC)
    drift = (t - t(1)) .^ do;
    regressors{end + 1} = drift ./ drift(end); % As in Homer3: Rescale the polynomial to make its last value = 1
end

nR = length(regressors); % # of regressors

%% testing: plot regressors
figure; hold on
for ind = 1:nR
    plot(regressors{ind})
end

%% Design matrix

X = zeros(nT, nR); % Initialize design matrix
% Store regressors
for ind = 1:nR
    X(:, ind) = regressors{ind};
end

%% Solve for the GLM weights
% weights = inv(X'*X)* (X'*data);
weights = (X'*X)\(X'*data);

%% 
HRFonly = regressors{1}*weights(1, :);
HRFonlyRS = permute(reshape(HRFonly, nT, size(PDIallSF, 1), size(PDIallSF, 2), size(PDIallSF, 3)), [2, 3, 4, 1]);

mapHRF = reshape(weights(1, :), size(PDIallSF, 1), size(PDIallSF, 2), size(PDIallSF, 3));

%% test
figure; imagesc(squeeze(max(mapHRF, [], 1))')
figure; imagesc(squeeze(max(HRFonlyRS(:, :, :, 40), [], 1))')
figure; imagesc(squeeze(max(mean(HRFonlyRS, 4), [], 1))')
figure; plot(squeeze(-mean(HRFonlyRS, [1, 2, 3])))