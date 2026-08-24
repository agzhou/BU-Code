%% Description:
%   Calculate flow velocity from planar fUS data using the vUS method
%   (Tang et al., 2020)
% Inputs:
%   IQ: [z voxels, x voxels, frames] complex data matrix
% Outputs:
%   vUS: [z voxels, x voxels, 2 (xz components)] real velocity data matrix

%%
% function [vUS] = vUS_2D(IQ)
%% Only use if needed: convert Jianbo/Bingxue's acquisition parameters to something I can use
% First: manually load an IQ file, like: load('E:\PROJ_tlfUS\IQdata\0806_2021_BL3_vUS_run1(good)\IQ-10-5-5000-1000-1-BL3-1.mat')
if ~exist('P_old', 'var')
    P_old = P; clearvars P
end
[P] = oldP2P(P_old);

%% Add the Speckle tracking folder to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
AllenProcessingCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing\");
addpath(genpath(AllenProcessingCodePath))

%% Set up the High Pass Filter (parameters from the 2020 vUS paper)
HPF.fc = 25; % Cutoff frequency [Hz]
% 25 Hz corresponds to 1 mm/s

HPF.fs = P.frameRate; % Sampling frequency [Hz]
HPF.order = 4; % Butterworth filter order

[HPF.b, HPF.a] = butter(HPF.order, HPF.fc/(HPF.fs/2), 'high');

%% Define some parameters
% sigma = [113, 999999, 151].*1e-6; % 1/e PSF values [m] for the L22-14v probe at 15.625 MHz and 17 angles from -10 to 10 deg. The y component is set to some arbitrary positive number but it won't really be used. (G:\My Drive\Data\PSF Simulations\L22-14v PSF sim - 17 angles from -10 to 10 deg)
% sigma = [113, 151].*1e-6; % 1/e PSF values (x, z) [m] for the L22-14v probe at 15.625 MHz and 17 angles from -10 to 10 deg. The y component is set to some arbitrary positive number but it won't really be used. (G:\My Drive\Data\PSF Simulations\L22-14v PSF sim - 17 angles from -10 to 10 deg)
% sigma = [41.6564, 52.3236].*1e-6; % Intensity-based 1/e PSF values (x, z) [m] for the L22-14v probe at 15.625 MHz and 17 angles from -10 to 10 deg. The y component is set to some arbitrary positive number but it won't really be used. (G:\My Drive\Data\PSF Simulations\L22-14v PSF sim - 17 angles from -10 to 10 deg)
sigma = [58.9110, 73.9967].*1e-6; % Field-based 1/e PSF values (x, z) [m] for the L22-14v probe at 15.625 MHz and 17 angles from -10 to 10 deg. The y component is set to some arbitrary positive number but it won't really be used. (G:\My Drive\Data\PSF Simulations\L22-14v PSF sim - 17 angles from -10 to 10 deg)

fDim = 3; % Dimension of the data corresponding to frequency (or time)
zDim = 1; % Dimension of the data corresponding to z (axial direction)
xDim = 2; % Dimension of the data corresponding to x (lateral direction)

%% ========= 1. Preprocessing ========= %%
% IQ = squeeze(complex(IData, QData));
% clearvars IData QData

% sv_threshold_lower = 20; sv_threshold_upper = size(IQ, 3);

% Mask the region to actually process
[zpo, xpo, nfo] = size(IQ); % Original sizes
% figure; imagesc(squeeze(abs(IQ(:, :, 1))))
% zrange = 1:100;
zrange = 1:zpo;
% zrange = 20:140;
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

[IQf, noise] = applySVs1D(IQm, CM, SVs, V, sv_threshold_lower, sv_threshold_upper);

% 1.2 High pass filter (apply to the post-SVD clutter filtered data)
HPF.dim = length(size(IQf)); % Operate on the time dimension
IQf_HPF = filter(HPF.b, HPF.a, IQf, [], HPF.dim);

% Testing
% figure; imagesc(squeeze(abs(IQf(:, :, 1))))
temp = sum(abs(IQf).^2, 3);
figure; imagesc(temp .^ 0.5)
% tp = [128, 39]; % Test point
% tp = [87, 26]; % Test point
% figure; plot(squeeze(abs(IQf_HPF(tp(1), tp(2), :))))
% figure; plot(squeeze(real(IQf(tp(1), tp(2), :))))
% figure; plot(squeeze(real(IQf_HPF(tp(1), tp(2), :))))

%% ========= 2. Directional flow filtering ========= %%

% 2.1 Separate positive and negative frequencies
[IQf_separated, IQf_FT_separated, nFTpts] = separatePosNegFreqs(IQf_HPF); % Outputs are cell arrays in the order of: negative, positive, all frequencies
% frameDim = length(size(IQf)); % Get the dimension corresponding to time/frames

ctp = 1:length(IQf_FT_separated); % Indices of which frequency Components To Process (typically [1, 2, 3]: negative, positive, all)
ctp_labels = {"Down flows", "Up flows", "All flows"};

% Testing: plot the separated and full Fourier spectrums and reconstructed IQ signals
faxis = linspace(-P.frameRate/2, P.frameRate/2, nFTpts)';
% figure; plot(faxis, squeeze(abs(IQf_FT_separated{1}(tp(1), tp(2), :))))
% figure; plot(faxis, squeeze(abs(IQf_FT_separated{2}(tp(1), tp(2), :))))
% figure; plot(faxis, squeeze(abs(IQf_FT_separated{3}(tp(1), tp(2), :))))
% figure; plot(1:P.numFramesPerBuffer, squeeze(abs(IQf_separated{1}(tp(1), tp(2), tp(3), :))))
% figure; plot(1:P.numFramesPerBuffer, squeeze(abs(IQf_separated{2}(tp(1), tp(2), tp(3), :))))
% figure; plot(1:P.numFramesPerBuffer, squeeze(abs(IQf_separated{3}(tp(1), tp(2), tp(3), :))))

% 2.2 Mask out all frequencies outside of +/- 1100 Hz, and other bands
% considered 'system noise' -- see Jianbo's sysNoiseRemove.m
freqMask = abs(faxis) > 1100; % [Hz]
% *********** add the system noise stuff later ***********

% Create IQf_FT_separated_masked: the Fourier-transformed filtered IQ data,
% with some frequencies masked out
IQf_FT_separated_masked = cell(size(IQf_FT_separated));
IQf_separated_masked = cell(size(IQf_separated));
for j = ctp
    IQf_FT_separated_masked{j} = IQf_FT_separated{j};
    IQf_FT_separated_masked{j}(:, :, freqMask) = 0;

    IQf_separated_masked{j} = ifft(ifftshift(IQf_FT_separated_masked{j}, fDim), nFTpts, fDim);
end
% figure; plot(faxis, squeeze(abs(IQf_FT_separated{3}(tp(1), tp(2), :))), '-', 'LineWidth', 2)
% hold on
% plot(faxis, squeeze(abs(IQf_FT_separated_masked{3}(tp(1), tp(2), :))), '--', 'LineWidth', 1)
% hold off

%% ========= 3. Calculate g1 at the first couple time lags, for cleaning criteria ========= %%

g1_dirty = g1T(IQf_separated{3}, 4); % Use all frequencies

%% 3.5 Calculate g1
% startTau = 1; % Index for the first tau point (tau1) for subsequent analysis. Changed this from 2 to 1 on 7/8/26 because I changed the g1T.m function to output g1 starting from tau = tau1 instead of tau = 0.
startTau = 2; % Index for the first tau point (tau1) for subsequent analysis.

nTau = ceil(20e-3 *P.frameRate); % # of time lags to consider; empirically set by assuming all g1 for voxels containing actual flow decay within 10 ms
% nTau = ceil(100e-3 *P.frameRate); % # of time lags to consider; empirically set by assuming all g1 for voxels containing actual flow decay within 10 ms
tau = (0:nTau - 1)' ./ P.frameRate; % Time lag vector [s]

% g1neg = g1T(IQf_separated{1}, nTau + startTau - 1); % Add the startTau-1 because the values start at startTau, but we still want nTau points total
% g1pos = g1T(IQf_separated{2}, nTau + startTau - 1); % Add the startTau-1 because the values start at startTau, but we still want nTau points total

% Store g1 for each frequency component in a cell array
g1 = cell(size(IQf_separated));

for j = ctp
    % g1{j} = g1T(IQf_separated{j}, nTau); % Use the base filtered IQ
    g1{j} = g1T(IQf_separated_masked{j}, nTau); % Use the filtered IQ with system noise removed
end

% Testing
figure; plot(tau, squeeze(abs(g1{1}(tp(1), tp(2), :))), '-o')
figure; plot(tau, squeeze(abs(g1{3}(tp(1), tp(2), :))), '-o')

figure; imagesc(squeeze(mean(abs(IQf_separated{1}), fDim))); title('Down flow')
figure; imagesc(squeeze(mean(abs(IQf_separated{2}), fDim))); title('Up flow')
figure; imagesc(squeeze(mean(abs(IQf_separated{3}), fDim))); title('All flow')

%% Create a struct for all the relevant processing parameters
dimensionality = 2; % 2D data
frameRate = P.frameRate;
wl = P.wl;
k0 = 2*pi/wl;
PP = createStruct(zp, xp, nf, nTau, xDim, zDim, fDim, dimensionality, faxis, freqMask, frameRate, wl, k0); % Processing Parameters ======> adjust as needed

%% ========= 4. Clean data ========= %%

% % 4.1 Screen voxels for noisiness, through |g1(tau1)|
% g1_tau1_threshold = 0.2;
% 
% g1_tau1_mask = cell(size(IQf_separated)); % Cell array of masks using the g1(tau1) threshold
% for j = ctp
%     g1_tau1_mask{j} = abs(squeeze(g1{j}(:, :, startTau))) > g1_tau1_threshold; % Use index 2 because index 1 corresponds to tau = 0
% end
% 
% % Testing/visualization
% figure; plot(squeeze(abs(g1{1}(tp(1), tp(2), :))), '-o'); title('Negative frequencies')
% figure; plot(squeeze(abs(g1{2}(tp(1), tp(2), :))), '-o'); title('Positive frequencies')
% % volumeViewer(g1_tau1_mask{1})
% % volumeViewer(g1_tau1_mask{2})
% 
% % 4.2 Apply mask
% % ...


% 4.1 Frequency-based SNR: whole frequency spectrum
[fbSNR, fbSNR_mask] = spectralSNR(IQf_FT_separated_masked{3}, IQf_FT_separated{3}, PP, 'full');

% 4.2 g1-based SNR: whole frequency spectrum
% % gR = mean(real(g1_dirty(:, :, 2:3)), fDim); % Average real(g1) over the first two time lags
% gR = mean(abs(g1_dirty(:, :, 2:3)), fDim); % Average abs(g1) over the first two time lags
% gR_pixel_avg = mean(gR, [zDim, xDim]);
% gR_std = std(gR, 0, [1, 2]);
% gR_mask = gR > max( (gR_pixel_avg - 0.4*gR_std), 0.08 );
% g1_express_mask = abs(g1_dirty(:, :, 2)) > 0.4; % Immediate pass mask according to |g1(tau1)| > threshold
[g1SNR, g1SNR_mask, g1SNR_express_mask] = g1BasedSNR(g1_dirty, PP, 'full');

% 4.3 Create the overall whole-frequency-spectrum mask (of pixels to keep)
overall_mask = and( or(fbSNR_mask, g1SNR_express_mask), g1SNR_mask);
figure; imagesc(overall_mask); title('Overall whole-frequency-spectrum mask')

%% ========= 5. Fit vUS ========= %%
% Initial guesses for parameters; separate fitting for negative and positive frequencies (down and up flows)

% CHANGE THIS LATER, WHEN I ACTUALLY IMPLEMENT VOXEL SCREENING!!!!!!!!!!!!!!!!!!
% num_voxels = size(g1neg, 1)*size(g1neg, 2)*size(g1neg, 3);
ps = size(IQf); ps = ps(1:end-1); % Plane size [voxels]
num_voxels = size(IQf, 1)*size(IQf, 2);

% Create new variables for experimental g1, with spatial dimensions stacked
g1_exp = cell(size(IQf_separated)); % Cell array of experimental g1 data with spatial dimensions vectorized/stacked
for j = ctp
    g1_exp{j} = reshape(g1{j}, num_voxels, nTau);
end

t1i = 2; % Index for tau1 --> 2 for my code, because it calculates g1 starting at tau = 0

% ---- Loop through directional components and go through the fitting process ---- %
% for j = ctp
% for j = 1:2 % Fit only negative and positive frequencies (down and up flows)
for j = 2
    % ---- Create masks for this direction's signal ---- %
    [fbSNR_j, fbSNR_mask_j, fbSNR_express_mask_j] = spectralSNR(IQf_FT_separated_masked{j}, IQf_FT_separated{j}, PP, 'half');
    [g1SNR_j, g1SNR_mask_j, g1SNR_express_mask] = g1BasedSNR(g1{j}, PP, 'half');
    [pnSNR_j, pnSNR_mask_j] = pnSpectralSNR(IQf_FT_separated_masked, PP, j);
    overall_mask_j = and( and( and(or(pnSNR_mask_j, fbSNR_express_mask_j), or(fbSNR_mask_j, g1SNR_express_mask)), g1SNR_mask_j), overall_mask);
    overall_mask_stacked_j = stackData(overall_mask_j, PP);
    % figure; imagesc(overall_mask_j)

    % ---- Adjust the g1 for this direction's signal ---- %
    % Create masks for potentially bad pixels
    g1adj_mask1_j = abs( abs(g1{j}(:, :, t1i)) - abs(g1{j}(:, :, t1i+1)) ) > 2.*abs( abs(g1{j}(:, :, t1i+1)) - abs(g1{j}(:, :, t1i+2)) ); % Flag a pixel if the |g1| drop from tau1 → tau2 is more than double the drop from tau2 → tau3 (is there some extra noise decorrelation in that first interval)
    g1adj_mask2_j = and( and( abs(g1{j}(:, :, t1i)) > 0.55, abs(g1{j}(:, :, t1i + 1)) < 0.25 ), abs(g1{j}(:, :, t1i + 1)) < abs(g1{j}(:, :, t1i + 2)) ); % (|g1(tau1)| > 0.55) AND (|g1(tau2)| < 0.25) AND (|g1(tau2)| < |g1(tau3)|) --> Flag a pixel if the |g1| at tau1 is high, low at tau2, and then goes back up at tau3 (which would be strange)
    % g1adj_mask_j = or(g1adj_mask1_j, g1adj_mask2_j);
    g1adj_mask_j = g1adj_mask2_j; % TESTING
    % Adjust g1(tau1) for these potentially bad pixels
    g1tau1_temp_j = (1 - g1adj_mask_j).*squeeze(g1{j}(:, :, t1i)) + (g1adj_mask_j).*( g1{j}(:, :, t1i + 1) + complex( abs(real(g1{j}(:, :, t1i) - g1{j}(:, :, t1i + 1))), imag(g1{j}(:, :, t1i + 1) - g1{j}(:, :, t1i + 2)) ) );
    % figure; imagesc(abs(g1tau1_temp_j)); clim([0, 1]); colorbar
    tp = [93, 174]; % test point
    % figure; plot(squeeze(real(g1{j}(tp(1), tp(2), :))), '-o')
    % figure; plot(squeeze(abs(g1{j}(tp(1), tp(2), :))), '-o')
    % temp = repmat(g1adj_mask_j, [1, 1, nTau]);
    g1adj_j = g1{j}; g1adj_j(:, :, t1i) = g1tau1_temp_j;
    % figure; plot(squeeze(abs(g1adj_j(tp(1), tp(2), :))), '-o')
    g1adj_stacked_j = stackData(g1adj_j, PP);

    % ---- Find initial guesses for fit parameters, for this direction's signal ---- %
    % Static (DC) component -- complex valued
    RotCtr_j = FindCOR( g1adj_stacked_j(:, round(nTau/2):end) ); % [nz*nx, nTau] -- for each pixel, take its last ½ of values (where the complex g1 spiral has theoretically started to slow down and look like a circle) and fit a circle to those points using FindCOR.m. The resulting center point of the fit is theoretically the center/end point of the complex g1 spiral, which represents the steady-state value. Take the real component of that output and use either this value if positive, or 0
    DCR0_v1_j = max(real(RotCtr_j), 0); % Version 1 of the DC component's Real component
    DCR0_v2_j = max(mean( real(g1adj_stacked_j(:, floor(end*2/3):end)), 2 ), 0); % Version 2 of the DC component's Real component -- for each pixel, take the temporal mean of the real components of its last ⅓ of values (where there is theoretically a steady-state/plateau), and use either this value if positive, or 0 otherwise.
    % DCR0_j = min(DCR0_v1_j, DCR0_v2_j); % Take the minimum of the two guesses above [nz*nx, nTau]
    % DCR0_j = complex(min(real(DCR0_v1_j), real(DCR0_v2_j)), min(imag(DCR0_v1_j), imag(DCR0_v2_j))); % Take the minimum of the two guesses above [nz*nx, nTau]
    DCR0_j = DCR0_v1_j; % Testing

    % Absolute "error" due to the noise decorrelation at tau1
    tau1_decorr_drop_j = 1 - abs(g1adj_stacked_j(:, t1i)); % How much does |g1(tau1)| drop from 1 (This is not used as its own explicit parameter)
    
    % Dynamic (noise decorrelation) component -- 'F' in the vUS paper
    FR_j = max(min(1 - abs(DCR0_j) - tau1_decorr_drop_j, 1), 0); % F Real component, clamped to [0, 1]

    % Axial component of the blood flow's group velocity -- v_zgp
    [Vz0, tau_V] = findVzPhaseDiff(stackData(g1{j}, PP), PP, DCR0_j); % v_zgp [m/s]
    % [Vz0, tau_V] = findVzPhaseDiff(g1adj_stacked_j, PP); % v_zgp [m/s]

    % Mesh method for finding v_xgp0, p0
    [v_zgp0, v_xgp0, p0, DC0, F0, R20] = InitvUS2DParamsWithMesh(g1adj_stacked_j, Vz0, DCR0_j, FR_j, PP, sigma, tau);

    % Testing: vessel angle for v_xgp0
    spacing = [1, 1];
    sigmas = [1:0.5:10];
    tau = 0.5;
    brightondark = true;
    vesselnessThreshold = 0.01;
    minBranchLengthPix = 5;
    minSegLengthPix = 9;
    PDI_j = sum(abs(IQf_separated_masked{j}).^2, fDim); % PDI for this direction
    gamma = 0.5;
    figure; imagesc(PDI_j .^ gamma)
    PDI_US_j = imresize(PDI_j, 5, "bilinear");
    figure; imagesc(PDI_US_j .^ gamma)
    % PDIN = PDI ./ max(PDIN, [], 'all'); % Normalized PDI [0, 1]
    % [vessels, angleMap, vesselness, vesselMask, segLabel, skel] = vesselAngle2D(abs(unstackData(Vz03, PP)) .* 1, sigmas, spacing, tau, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix);
    [vessels, angleMap, vesselness, vesselMask, segLabel, skel] = vesselAngle2D(abs(PDI_US_j .^ gamma), sigmas, spacing, tau, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix);
    figure; imagesc(vesselness)
    figure; imagesc(skel)
    figure; imagesc(angleMap); colormap hsv

    PDI3_US = imresize(PDI3, 5, "bilinear");
    figure; imagesc(PDI3_US .^ gamma)
    % PDIN = PDI ./ max(PDIN, [], 'all'); % Normalized PDI [0, 1]
    % [vessels, angleMap, vesselness, vesselMask, segLabel, skel] = vesselAngle2D(abs(unstackData(Vz03, PP)) .* 1, sigmas, spacing, tau, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix);
    [vessels, angleMap, vesselness, vesselMask, segLabel, skel] = vesselAngle2D(abs(PDI3_US .^ gamma), sigmas, spacing, tau, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix);
    figure; imagesc(vesselness)
    figure; imagesc(skel)
    figure; imagesc(angleMap); colormap hsv

    % ---- Fit this direction's signal ---- %
    % anon_fun = @(x) g1vUS2D_Jac(x, tau, sigma, PP.k0);
    anon_fun = @(x) g1vUS2D_vec_split(x, tau, sigma, PP.k0, useF, useDC);
    % opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true);
    opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', false);

    % v_xgp = zeros(PP.zp, PP.xp);
    % v_zgp = zeros(PP.zp, PP.xp);
    % p = zeros(PP.zp, PP.xp);
    % F = zeros(PP.zp, PP.xp);
    % DC = zeros(PP.zp, PP.xp);
    v_xgp_stacked = zeros(PP.zp*PP.xp, 1);
    v_zgp_stacked = zeros(PP.zp*PP.xp, 1);
    p_stacked = zeros(PP.zp*PP.xp, 1);
    F_stacked = zeros(PP.zp*PP.xp, 1);
    DC_stacked = zeros(PP.zp*PP.xp, 1);

    tic
    % for vi = 1:num_voxels % voxel index
    % for vi = 1:300
    for vi = ind
        % [zi, xi] = 
        if overall_mask_stacked_j(vi)
            x0 = [v_xgp0(vi), v_zgp0(vi), p0(vi), F0(vi), DC0(vi)];
            % **** NEED TO FIX lb AND ub --> VELOCITIES CAN BE NEGATIVE ****
            % lb = [x0(1)*0.5, x0(2)*0.5, max(p0(vi) - 0.2, 0), max(F0(vi) - 0.2, 0), max(DC0(vi) - 0.2, 0)]; % TESTING
            % ub = [x0(1)*1.5, x0(2)*1.5, min(p0(vi) + 0.2, 1), min(F0(vi) + 0.2, 1), min(DC0(vi) + 0.2, 1)]; % TESTING
            lb = [x0(1)*(1 - 0.25), x0(2)*(1 - 0.25), max(p0(vi) - 0.2, 0), max(F0(vi) - 0.2, 0), max(DC0(vi) - 0.2, 0)]; % TESTING
            ub = [x0(1)*(1 + 0.25), x0(2)*(1 + 0.25), min(p0(vi) + 0.2, 1), min(F0(vi) + 0.2, 1), min(DC0(vi) + 0.2, 1)]; % TESTING
            x = lsqnonlin(anon_fun, x0, lb, ub, opts); % x = [v_xgp, v_zgp, p, F, DC]
            v_xgp_stacked(vi) = x(1);
            v_zgp_stacked(vi) = x(2);
            p_stacked(vi) = x(3);
            F_stacked(vi) = x(4);
            DC_stacked(vi) = x(5);
        end
    end
    toc

    v_xgp = unstackData(v_xgp_stacked, PP);
    v_zgp = unstackData(v_zgp_stacked, PP);
    p = unstackData(p_stacked, PP);
    F = unstackData(F_stacked, PP);
    DC = unstackData(DC_stacked, PP);

    test = g1vUS2D_Jac(x, tau, sigma, k0);
    figure; plot(tau, abs(g1adj_stacked_j(vi, :)), tau, abs(test))
end








%%
% % Create structs that store parameters (including initial guesses) for the vUS fitting
% 
% % General stuff
% vf_gen = struct(); % vUS fitting struct
% vf_gen.k0 = 2*pi/P.wl; % Angular wavenumber [rad/m]
% % p_all.v_xgp_range = [1, 30]./1e3; % Min and max values [m/s] for v_xgp to use in the mesh initial guessing
% vf_gen.v_xgp_range = [1, 10]./1e3; % Min and max values [m/s] for v_xgp to use in the mesh initial guessing
% vf_gen.v_xgp_step = 1e-3; % Increment for the v_xgp grid [m/s]
% vf_gen.v_xgp_grid = vf_gen.v_xgp_range(1):vf_gen.v_xgp_step:vf_gen.v_xgp_range(2);
% 
% % p_all.v_ygp_range = [1, 30]./1e3; % Min and max values [m/s] for v_ygp to use in the mesh initial guessing
% vf_gen.v_ygp_range = [1, 10]./1e3; % Min and max values [m/s] for v_ygp to use in the mesh initial guessing
% vf_gen.v_ygp_step = 1e-3; % Increment for the v_ygp grid [m/s]
% vf_gen.v_ygp_grid = vf_gen.v_ygp_range(1):vf_gen.v_ygp_step:vf_gen.v_ygp_range(2);
% 
% vf_gen.p_range = [1, 0]; % Min and max values [unitless] for p to use in the mesh initial guessing
% vf_gen.p_step = -0.1; % Increment for the p grid
% vf_gen.p_grid = vf_gen.p_range(1):vf_gen.p_step:vf_gen.p_range(2);
% 
% vf_gen.meshgrid = meshgrid(vf_gen.v_xgp_grid, vf_gen.v_ygp_grid, vf_gen.p_grid); % Create mesh for the guessing of initial values for v_xgp0, v_ygp0, and p0
% 
% vf = struct(); % vUS fitting struct; one component per frequency component with the following fields each
% for j = ctp
%     vf(j).F0 = reshape( abs(squeeze(g1{j}(:, :, startTau))), num_voxels, 1); % Initial guess for F
%     vf(j).tau_V = findFirstLocalMin(g1_exp{j}, nTau, 'smooth') ./ P.frameRate; % Time lag [s] at which g1 reaches its first minimum, per voxel. Here, I'm reshaping g1 to pass in a matrix where voxels are stacked.
%     vf(j).v_zgp0 = squeeze(P.wl./(4.*vf(j).tau_V)); % Initial guess for v_zgp (Eq. 16)
% end
% 
% % Rmesh = 
% 
% %% Fit voxels individually
% tic
% 
% % useF = false; % Use the F parameter or not
% useF = true;
% 
% useDC = true; % Use a DC offset for fitting or not
% % useDC = false;
% 
% % tp = [40, 43, 87];
% % tp = [40, 39, 87];
% % tp = [10, 70, 142];
% 
% % fit_roi = {tp(1), tp(2), tp(3)}; % Define a spatial region to fit within
% % k = [2, 5, 10];
% % fit_roi = {tp(1) - k(1) : tp(1) + k(1), tp(2) - k(2) : tp(2) + k(2), tp(3) - k(3):tp(3) + k(3)}; % Define a spatial region to fit within
% % fit_roi = {50:70, 130:170};
% fit_roi = {1:ps(1), 1:ps(2)}; % Full volume
% 
% % Fitting options
% options = optimoptions('lsqcurvefit', 'Display', 'off');
% if useF
%     if useDC
%         % ****** NEED TO CHANGE THE BELOW TO ONLY FIT X AND Z ****** %
%         % [vx, vz, p, F, DC offset]
%         lb = [0, 0, 0, 0, 0];             % Lower bounds for parameters [SI units]
%         ub = [50e-3, 50e-3, 1, 1, 1]; % Upper bounds for parameters [SI units]
%         % ub = [100e-3, 100e-3, 1, 1, 1]; % Upper bounds for parameters [SI units]
%     else
%         % [vx, vz, p, F]
%         lb = [0, 0, 0, 0];             % Lower bounds for parameters [SI units]
%         ub = [50e-3, 50e-3, 1, 1]; % Upper bounds for parameters [SI units]
%         % ub = [100e-3, 100e-3, 1, 1]; % Upper bounds for parameters [SI units]
%     end
% else
%     % [vx, vz, p]
%     lb = [0, 0, 0];             % Lower bounds for parameters [SI units]
%     ub = [50e-3, 50e-3, 1]; % Upper bounds for parameters [SI units]
%     % ub = [100e-3, 100e-3, 1]; % Upper bounds for parameters [SI units]
% end
% % fit_roi = {11, 81};
% fit_roi = {1:ps(1), 1:ps(2)}; % Full volume
% 
% % Fitting options
% options = optimoptions('lsqcurvefit', 'Display', 'off');
% 
% % Create variables to store vUS fitting results
% vUS = cell(size(IQf_separated)); % vUS results for each frequency component
% p = cell(size(IQf_separated)); % p results for each frequency component
% if useF
%     F = cell(size(IQf_separated)); % F results for each frequency component
% end
% if useDC
%     DC = cell(size(IQf_separated)); % DC component results for each frequency component
% end
% 
% for j = ctp
%     vUS{j} = zeros([ps, 2]); % [z pix, x pix, x or z velocity component]
%     p{j} = zeros(ps); % [z pix, x pix]
%     if useF
%         F{j} = zeros(ps);
%     end
%     if useDC
%         DC{j} = zeros(ps);
%     end
% end
% 
% 
% tic
% for zi = fit_roi{1}
%     for xi = fit_roi{2}
%         ind = sub2ind(ps, zi, xi);
%         for j = ctp % For each pixel, go through and fit each 
% 
%             % Only fit if the voxel meets some criterion (after screening). For testing, don't do this.
%             % if 1
%             % if vf_all.F0(ind) > 0.2
%             if overall_mask(zi, xi)
% 
%                 % % Fitting the complex data all-in-one
%                 % % x0_neg = [p_neg.F0, p_neg.p0, p_neg.v_xgp0, p_neg.v_ygp0, p_neg.v_zgp0]; % ICs: [F0, p0, v_xgp0, v_ygp0, v_zgp0]
%                 % % x0_neg = [p_neg.p0(ind), p_neg.v_xgp0(ind), p_neg.v_ygp0(ind), p_neg.v_zgp0(ind)]; % ICs: [F0, p0, v_xgp0, v_ygp0, v_zgp0]
%                 % x0_neg = [1, 10e-3, 10e-3, vf_neg.v_zgp0(ind)]; % ICs: [p0, v_xgp0, v_ygp0, v_zgp0]
%                 % % lb = [0, 0, 0, 0];             % Lower bounds for parameters [SI units]
%                 % % ub = [1, 50e-3, 50e-3, 50e-3]; % Upper bounds for parameters [SI units]
%                 % f_temp = @(x, tau) g1vUS3D_vec(x, tau, sigma, vf_gen.k0, useF); % Use "anonymous function" to pass in the g1 vUS model function to the fitting
%                 % ydata = squeeze(g1neg_exp(ind, :));
%                 % x_neg = lsqcurvefit(f_temp, x0_neg, squeeze(tau(1:nTau)), ydata);
%                 % % x_neg = lsqcurvefit(f_temp, x0_neg, squeeze(tau(1:nTau)), squeeze(g1neg_exp(ind, :)), lb, ub);
% 
%                 % Get the data and initial conditions for this voxel
%                 if useF
%                     ydata = squeeze(g1_exp{j}(ind, 2:end)); ydata = ydata(:); % Experimental data to fit to
% 
%                     if useDC
%                         x0 = [10e-3, vf(j).v_zgp0(ind), 1, vf(j).F0(ind), 0]; % ICs: [v_xgp0, v_zgp0, p0, F0, DC]
%                     else
%                         x0 = [10e-3, vf(j).v_zgp0(ind), 1, vf(j).F0(ind)]; % ICs: [v_xgp0, v_zgp0, p0, F0]
%                     end
%                 else
%                     ydata = squeeze(g1_exp{j}(ind, :)); ydata = ydata(:); % make sure it's a column vector
%                     x0 = [10e-3, vf(j).v_zgp0(ind), 1]; % ICs: [v_xgp0, v_zgp0, p0]
% 
%                 end
% 
%                 % % Splitting the real and complex components
%                 % ydata_neg = ydata_neg(:);
%                 % ydata_neg_split = [real(ydata_neg), imag(ydata_neg)];
%                 % 
%                 % ydata_all = ydata_all(:);
%                 % ydata_all_split = [real(ydata_all), imag(ydata_all)];
% 
%                 % % Magnitude only
%                 % ydata_neg = abs(ydata_neg(:));
%                 % ydata_neg_split = [real(ydata_neg), zeros(size(ydata_neg))];
%                 % 
%                 % ydata_all = abs(ydata_all(:));
%                 % ydata_all_split = [real(ydata_all), zeros(size(ydata_all))];
% 
%                 % Perform the fitting
%                 % f_temp = @(x, tau) g1vUS3D_vec_split(x, tau, sigma, vf_gen.k0, useF); % Use "anonymous function" to pass in the g1 vUS model function to the fitting
%                 % % x_neg = lsqcurvefit(f_temp, x0_neg, squeeze(tau(1:nTau)), ydata_split);
%                 % % x_all = lsqcurvefit(f_temp, x0_all, squeeze(tau(1:nTau)), ydata_all_split, lb, ub, options);
%                 % if useF
%                 %     x_neg = lsqcurvefit(f_temp, x0_neg, squeeze(tau(2:nTau)), ydata_neg_split, lb, ub, options);
%                 %     F_neg(xi, yi, zi) = x_neg(5);
%                 % 
%                 %     x_all = lsqcurvefit(f_temp, x0_all, squeeze(tau(2:nTau)), ydata_all_split, lb, ub, options);
%                 %     F_all(xi, yi, zi) = x_all(5);
%                 % else
%                 %     x_neg = lsqcurvefit(f_temp, x0_neg, squeeze(tau(1:nTau)), ydata_neg_split, lb, ub, options);
%                 % 
%                 %     x_all = lsqcurvefit(f_temp, x0_all, squeeze(tau(1:nTau)), ydata_all_split, lb, ub, options);
%                 % end
%                 % 
%                 % % Store more results
%                 % vUS_neg(xi, yi, zi, :) = x_neg(2:4);
%                 % p_neg(xi, yi, zi) = x_neg(1);
%                 % 
%                 % vUS_all(xi, yi, zi, :) = x_all(2:4);
%                 % p_all(xi, yi, zi) = x_all(1);
% 
% 
%                 % ======== Fit only |g1| ========
%                 f_temp = @(x, tau) g1vUS2D_mag_vec(x, tau, sigma, vf_gen.k0, useF, useDC); % Use "anonymous function" to pass in the g1 vUS model function to the fitting
%                 if useF
%                     x = lsqcurvefit(f_temp, x0, squeeze(tau(2:nTau)), abs(ydata), lb, ub, options);
%                     F{j}(zi, xi) = x(4);
% 
%                     if useDC
%                         DC{j}(zi, xi) = x(5);
%                     end
% 
%                 else
%                     x = lsqcurvefit(f_temp, x0, squeeze(tau(1:nTau)), abs(ydata), lb, ub, options);
%                 end
% 
%                 % Store more results
%                 vUS{j}(zi, xi, :) = x(1:2);
%                 p{j}(zi, xi) = x(3);
% 
%             else % If we don't fit that pixel due to noisiness, set the velocities and other parameters to 0
%                 vUS{j}(zi, xi, :) = [0, 0];
%                 p{j}(zi, xi) = 0;
%                 if useF
%                     F{j}(zi, xi) = 0;
%                 end
%                 if useDC
%                     DC{j}(zi, xi) = 0;
%                 end
% 
%             end
%         end
%     end
% end
% toc
% 
% % % testg1 = g1vUS3D_vec(x_neg, tau, sigma, vf_gen.k0, useF);
% % % figure; plot(ydata_neg, '-x', 'LineWidth', 2); hold on; plot(testg1, '-o', 'LineWidth', 1); hold off
% % % figure; plot(abs(ydata_neg), '-x', 'LineWidth', 2); hold on; plot(abs(testg1), ':', 'LineWidth', 1); hold off
% % % testspeed = sqrt(sum(x_neg(2:4).^2))
% % 
% % testg1 = g1vUS2D_vec(x_all, tau(2:end), sigma, vf_gen.k0, useF, useDC);
% % figure; plot(ydata_all); hold on; plot(testg1); hold off
% % figure; plot(abs(ydata_all)); hold on; plot(abs(testg1)); hold off
% % testspeed = sqrt(sum(x_all(1:2).^2))
% % 
% % testg1 = g1vUS2D_vec(x, tau(2:end), sigma, vf_gen.k0, useF, useDC);
% % figure; plot(ydata); hold on; plot(testg1); hold off
% % figure; plot(abs(ydata)); hold on; plot(abs(testg1)); hold off
% % testspeed = sqrt(sum(x(1:2).^2))
% % 
% % testg1 = g1vUS2D_vec(x_pos, tau(2:end), sigma, vf_gen.k0, useF, useDC);
% % figure; plot(ydata_pos); hold on; plot(testg1); hold off
% % figure; plot(abs(ydata_pos)); hold on; plot(abs(testg1)); hold off
% % testspeed = sqrt(sum(x_pos(1:2).^2))
% 
% % testg1 = g1vUS3D_vec(x_all, tau(1:end), sigma, vf_gen.k0, useF);
% % figure; plot(ydata_all); hold on; plot(testg1); hold off
% % figure; plot(abs(ydata_all)); hold on; plot(abs(testg1)); hold off
% % testspeed = sqrt(sum(x_all(2:4).^2))

%% Testing: visualize vUS results
vUS_speed = cell(size(vUS));
for j = ctp
    vUS_speed{j} = squeeze(sqrt(sum(vUS{j}(fit_roi{1}, fit_roi{2}, :).^2, 3)));
end

for j = ctp
    figure; imagesc(vUS_speed{j}); title(ctp_labels{j}); colorbar
end

% figure; imagesc(squeeze(max(PDI, [], 1))')
PDI = squeeze(mean(abs(IQf_HPF).^2, 3));
figure; imagesc(squeeze(PDI(fit_roi{1}, fit_roi{2})) .^ 0.5); title('PDI')

%% Plot the vUS fit at a test point
% testpt = [43, 11];
% testpt = [60, 133];
% testpt = [11, 81];
testpt = [43, 65];
for j = ctp
    if useF
        if useDC
            plotg1pt(testpt(1), testpt(2), useF, useDC, tau, sigma, vf_gen.k0, g1{j}, vUS{j}, p{j}, F{j}, DC{j})
        else
            plotg1pt(testpt(1), testpt(2), useF, useDC, tau, sigma, vf_gen.k0, g1{j}, vUS{j}, p{j}, F{j})
        end
    else
        plotg1pt(testpt(1), testpt(2), useF, useDC, tau, sigma, vf_gen.k0, g1{j}, vUS{j}, p{j})
    end
end

findfigs
%% Overlay up and down flows

% Load Jianbo's colormaps
[VzCmap, VzCmapDn, VzCmapUp, pdiCmapUp, PhtmCmap] = Colormaps_fUS;
vUSMapFig = figure;

% Plot with two linked axes (one for up Z, other for down Z)
% % hold on
% vCrange = [-maxSpeedExpectedMMPerS, maxSpeedExpectedMMPerS];
figure(vUSMapFig)
h1 = axes;
imagesc(h1, vUS_speed{2})
% alpha(h1, double(abs(zvUpMap) > 1))
% alpha(h1, 1)
alpha(h1, vUS_speed_pos ./ max(vUS_speed{2}, [], 'all'));
colormap(vUSMapFig, flipud(VzCmapUp))
% caxis(vCrange);
axis tight
colorbar
hold on

% figure
h2 = axes;
imagesc(h2, vUS_speed{1})
% alpha(h2, double(abs(vUS_speed_neg) > 1))
% alpha(h2, 0.5)
alpha(h2, vUS_speed_neg ./ max(vUS_speed{1}, [], 'all'));
colormap(vUSMapFig, VzCmapDn)
% caxis(vCrange);
axis tight
cb = colorbar;
axis off
linkaxes([h1, h2]);
% ylabel(cb, 'Speed (m/s)') % label colorbar

% clim([])
findfigs

%% Helper functions

% Plot the vUS fit at one point (z, x) against the experimental data. Do
% up, down, all flow directions separately
function plotg1pt(z, x, useF, useDC, tau, sigma, k0, g1exp, vUS, p, varargin)
    if nargin > 10
        F = varargin{1};
        if nargin > 11
            DC = varargin{2};
            X = paramsToX(z, x, useF, useDC, vUS, p, F, DC);
        else
            X = paramsToX(z, x, useF, useDC, vUS, p, F);
        end
    end
    
    testg1 = g1vUS2D_vec(X, tau(2:end), sigma, k0, useF, useDC);
    % figure; plot(squeeze(g1exp(z, x, :))); hold on; plot(testg1); hold off
    figure; plot(tau, squeeze(abs(g1exp(z, x, :))), '-x', 'LineWidth', 2); hold on; plot(tau(2:end), abs(testg1), '-o', 'LineWidth', 2); hold off; ylabel('|g1|'); xlabel('Time lag [s]'); legend('Data', 'Fit')
    % testspeed = sqrt(sum(X(1:2).^2))
    
end

% Convert fitted params into a vector X
%   Optional inputs: F, DC
function X = paramsToX(z, x, useF, useDC, vUS, p, varargin)
    if nargin > 6
        F = varargin{1};
        if nargin > 7
            DC = varargin{2};
        end
    end

    if useF
        if useDC
            % error('Have not added this in the code yet')
            X = [squeeze(vUS(z, x, :)); p(z, x); F(z, x); DC(z, x)];
        else
            X = [squeeze(vUS(z, x, :)); p(z, x); F(z, x)];
        end
    else
        X = [squeeze(vUS(z, x, :)); p(z, x)];
    end
end