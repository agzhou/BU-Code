%% Description:
%   Calculate flow velocity from planar fUS data using the vUS method
%   (Tang et al., 2020)
% Inputs:
%   IQ: [z voxels, x voxels, frames] complex data matrix
% Outputs:
%   vUS: [z voxels, x voxels, 2 (xz components)] real velocity data matrix

%%
% function [vUS] = vUS_2D(IQ)

%% Add the Speckle tracking folder to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
AllenSpeckleTrackingCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing\Speckle tracking");
addpath(AllenSpeckleTrackingCodePath)

%% Set up the High Pass Filter (parameters from the 2020 vUS paper)
HPF.fc = 25; % Cutoff frequency [Hz]
% 25 Hz corresponds to 1 mm/s

HPF.fs = P.frameRate; % Sampling frequency [Hz]
HPF.order = 4; % Butterworth filter order

[HPF.b, HPF.a] = butter(HPF.order, HPF.fc/(HPF.fs/2), 'high');

%% ========= 1. Preprocessing ========= %%

% 1.1 SVD clutter filter
%     [PP, EVs, V_sort] = getSVs2D(IQ);
[zp, xp, nf] = size(IQ);
PP = reshape(IQ, [zp*xp, nf]);
tic
%     [U, S, V] = svd(PP); % Already sorted in decreasing order
[U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
SVs = diag(S);
%     disp('Full SVD done')
toc
disp('SVs decomposed')

[IQf, noise] = applySVs1D(IQ, PP, SVs, V, sv_threshold_lower, sv_threshold_upper);

% 1.2 High pass filter (apply to the post-SVD clutter filtered data)
HPF.dim = length(size(IQf)); % Operate on the time dimension
IQf_HPF = filter(HPF.b, HPF.a, IQf, [], HPF.dim);

% Testing
% figure; imagesc(squeeze(abs(IQf(:, :, 1))))
temp = sum(abs(IQf).^2, 3);
figure; imagesc(temp)
% tp = [128, 39]; % Test point
tp = [87, 26]; % Test point
figure; plot(squeeze(abs(IQf_HPF(tp(1), tp(2), :))))
% figure; plot(squeeze(real(IQf(tp(1), tp(2), :))))
% figure; plot(squeeze(real(IQf_HPF(tp(1), tp(2), :))))

%% ========= 2. Directional flow filtering ========= %%

% 2.1 Separate positive and negative frequencies
[IQf_separated, IQf_FT_separated, nFTpts] = separatePosNegFreqs(IQf_HPF); % Outputs are cell arrays in the order of: negative, positive, all frequencies
frameDim = length(size(IQf)); % Get the dimension corresponding to time/frames

% Testing: plot the separated and full Fourier spectrums and reconstructed IQ signals
faxis = linspace(-P.frameRate/2, P.frameRate/2, nFTpts)';
figure; plot(faxis, squeeze(abs(IQf_FT_separated{1}(tp(1), tp(2), :))))
figure; plot(faxis, squeeze(abs(IQf_FT_separated{2}(tp(1), tp(2), :))))
figure; plot(faxis, squeeze(abs(IQf_FT_separated{3}(tp(1), tp(2), :))))
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


% startTau = 1; % Index for the first tau point (tau1) for subsequent analysis. Changed this from 2 to 1 on 7/8/26 because I changed the g1T.m function to output g1 starting from tau = tau1 instead of tau = 0.
startTau = 2; % Index for the first tau point (tau1) for subsequent analysis.

nTau = ceil(20e-3 *P.frameRate); % # of time lags to consider; empirically set by assuming all g1 for voxels containing actual flow decay within 10 ms
tau = (0:nTau - 1)' ./ P.frameRate; % Time lag vector [s]

% g1neg = g1T(IQf_separated{1}, nTau + startTau - 1); % Add the startTau-1 because the values start at startTau, but we still want nTau points total
% g1pos = g1T(IQf_separated{2}, nTau + startTau - 1); % Add the startTau-1 because the values start at startTau, but we still want nTau points total

% Store g1 for each frequency component in a cell array
g1 = cell(size(IQf_separated));
ctp = 1:length(g1); % Indices of which frequency Components To Process (typically [1, 2, 3]: negative, positive, all)
ctp_labels = {"Down flows", "Up flows", "All flows"};

for j = ctp
    g1{j} = g1T(IQf_separated{j}, nTau);
end

% Testing
figure; plot(squeeze(abs(g1{1}(tp(1), tp(2), :))), '-o')
figure; plot(squeeze(abs(g1{3}(tp(1), tp(2), :))), '-o')

figure; imagesc(squeeze(mean(abs(IQf_separated{1}), frameDim))); title('Down flow')
figure; imagesc(squeeze(mean(abs(IQf_separated{2}), frameDim))); title('Up flow')
figure; imagesc(squeeze(mean(abs(IQf_separated{3}), frameDim))); title('All flow')

%% ========= 4. Clean data ========= %%

% 4.1 Screen voxels for noisiness, through |g1(tau1)|
g1_tau1_threshold = 0.2;

g1_tau1_mask = cell(size(IQf_separated)); % Cell array of masks using the g1(tau1) threshold
for j = ctp
    g1_tau1_mask{j} = abs(squeeze(g1{j}(:, :, startTau))) > g1_tau1_threshold; % Use index 2 because index 1 corresponds to tau = 0
end

% Testing/visualization
figure; plot(squeeze(abs(g1{1}(tp(1), tp(2), :))), '-o'); title('Negative frequencies')
figure; plot(squeeze(abs(g1{2}(tp(1), tp(2), :))), '-o'); title('Positive frequencies')
% volumeViewer(g1_tau1_mask{1})
% volumeViewer(g1_tau1_mask{2})

% 4.2 Apply mask
% ...

%% ========= 5. Fit vUS ========= %%
% Initial guesses for parameters; separate fitting for negative and positive frequencies (down and up flows)

% CHANGE THIS LATER, WHEN I ACTUALLY IMPLEMENT VOXEL SCREENING!!!!!!!!!!!!!!!!!!
% num_voxels = size(g1neg, 1)*size(g1neg, 2)*size(g1neg, 3);
ps = size(IQf); ps = ps(1:end-1); % Plane size [voxels]
num_voxels = size(IQf, 1)*size(IQf, 2);

g1_exp = cell(size(IQf_separated)); % Cell array of experimental g1 data with spatial dimensions vectorized/stacked
for j = ctp
    g1_exp{j} = reshape(g1{j}, num_voxels, nTau);
end

% TESTING
% sigma = [113, 999999, 151].*1e-6; % 1/e PSF values [m] for the L22-14v probe at 15.625 MHz and 17 angles from -10 to 10 deg. The y component is set to some arbitrary positive number but it won't really be used. (G:\My Drive\Data\PSF Simulations\L22-14v PSF sim - 17 angles from -10 to 10 deg)
% sigma = [113, 151].*1e-6; % 1/e PSF values (x, z) [m] for the L22-14v probe at 15.625 MHz and 17 angles from -10 to 10 deg. The y component is set to some arbitrary positive number but it won't really be used. (G:\My Drive\Data\PSF Simulations\L22-14v PSF sim - 17 angles from -10 to 10 deg)
% sigma = [41.6564, 52.3236].*1e-6; % Intensity-based 1/e PSF values (x, z) [m] for the L22-14v probe at 15.625 MHz and 17 angles from -10 to 10 deg. The y component is set to some arbitrary positive number but it won't really be used. (G:\My Drive\Data\PSF Simulations\L22-14v PSF sim - 17 angles from -10 to 10 deg)
sigma = [58.9110, 73.9967].*1e-6; % Field-based 1/e PSF values (x, z) [m] for the L22-14v probe at 15.625 MHz and 17 angles from -10 to 10 deg. The y component is set to some arbitrary positive number but it won't really be used. (G:\My Drive\Data\PSF Simulations\L22-14v PSF sim - 17 angles from -10 to 10 deg)

% Create structs that store parameters (including initial guesses) for the vUS fitting

% General stuff
vf_gen = struct(); % vUS fitting struct
vf_gen.k0 = 2*pi/P.wl; % Angular wavenumber [rad/m]
% p_all.v_xgp_range = [1, 30]./1e3; % Min and max values [m/s] for v_xgp to use in the mesh initial guessing
vf_gen.v_xgp_range = [1, 10]./1e3; % Min and max values [m/s] for v_xgp to use in the mesh initial guessing
vf_gen.v_xgp_step = 1e-3; % Increment for the v_xgp grid [m/s]
vf_gen.v_xgp_grid = vf_gen.v_xgp_range(1):vf_gen.v_xgp_step:vf_gen.v_xgp_range(2);

% p_all.v_ygp_range = [1, 30]./1e3; % Min and max values [m/s] for v_ygp to use in the mesh initial guessing
vf_gen.v_ygp_range = [1, 10]./1e3; % Min and max values [m/s] for v_ygp to use in the mesh initial guessing
vf_gen.v_ygp_step = 1e-3; % Increment for the v_ygp grid [m/s]
vf_gen.v_ygp_grid = vf_gen.v_ygp_range(1):vf_gen.v_ygp_step:vf_gen.v_ygp_range(2);

vf_gen.p_range = [1, 0]; % Min and max values [unitless] for p to use in the mesh initial guessing
vf_gen.p_step = -0.1; % Increment for the p grid
vf_gen.p_grid = vf_gen.p_range(1):vf_gen.p_step:vf_gen.p_range(2);

vf_gen.meshgrid = meshgrid(vf_gen.v_xgp_grid, vf_gen.v_ygp_grid, vf_gen.p_grid); % Create mesh for the guessing of initial values for v_xgp0, v_ygp0, and p0

vf = struct(); % vUS fitting struct; one component per frequency component with the following fields each
for j = ctp
    vf(j).F0 = reshape( abs(squeeze(g1{j}(:, :, startTau))), num_voxels, 1); % Initial guess for F
    vf(j).tau_V = findFirstLocalMin(g1_exp{j}, nTau, 'smooth') ./ P.frameRate; % Time lag [s] at which g1 reaches its first minimum, per voxel. Here, I'm reshaping g1 to pass in a matrix where voxels are stacked.
    vf(j).v_zgp0 = squeeze(P.wl./(4.*vf(j).tau_V)); % Initial guess for v_zgp (Eq. 16)
end

% Rmesh = 

%% Fit voxels individually
tic

% useF = false; % Use the F parameter or not
useF = true;

useDC = true; % Use a DC offset for fitting or not
% useDC = false;

% tp = [40, 43, 87];
% tp = [40, 39, 87];
% tp = [10, 70, 142];

% fit_roi = {tp(1), tp(2), tp(3)}; % Define a spatial region to fit within
% k = [2, 5, 10];
% fit_roi = {tp(1) - k(1) : tp(1) + k(1), tp(2) - k(2) : tp(2) + k(2), tp(3) - k(3):tp(3) + k(3)}; % Define a spatial region to fit within
% fit_roi = {50:70, 130:170};
fit_roi = {1:ps(1), 1:ps(2)}; % Full volume

% Fitting options
options = optimoptions('lsqcurvefit', 'Display', 'off');
if useF
    if useDC
        % ****** NEED TO CHANGE THE BELOW TO ONLY FIT X AND Z ****** %
        % [vx, vz, p, F, DC offset]
        lb = [0, 0, 0, 0, 0];             % Lower bounds for parameters [SI units]
        ub = [50e-3, 50e-3, 1, 1, 1]; % Upper bounds for parameters [SI units]
        % ub = [100e-3, 100e-3, 1, 1, 1]; % Upper bounds for parameters [SI units]
    else
        % [vx, vz, p, F]
        lb = [0, 0, 0, 0];             % Lower bounds for parameters [SI units]
        ub = [50e-3, 50e-3, 1, 1]; % Upper bounds for parameters [SI units]
        % ub = [100e-3, 100e-3, 1, 1]; % Upper bounds for parameters [SI units]
    end
else
    % [vx, vz, p]
    lb = [0, 0, 0];             % Lower bounds for parameters [SI units]
    ub = [50e-3, 50e-3, 1]; % Upper bounds for parameters [SI units]
    % ub = [100e-3, 100e-3, 1]; % Upper bounds for parameters [SI units]
end
                
% Create variables to store vUS fitting results
vUS = cell(size(IQf_separated)); % vUS results for each frequency component
p = cell(size(IQf_separated)); % p results for each frequency component
if useF
    F = cell(size(IQf_separated)); % F results for each frequency component
end
if useDC
    DC = cell(size(IQf_separated)); % DC component results for each frequency component
end

for j = ctp
    vUS{j} = zeros([ps, 2]); % [z pix, x pix, x or z velocity component]
    p{j} = zeros(ps); % [z pix, x pix]
    if useF
        F{j} = zeros(ps);
    end
    if useDC
        DC{j} = zeros(ps);
    end
end


tic
for zi = fit_roi{1}
    for xi = fit_roi{2}
        ind = sub2ind(ps, zi, xi);
        for j = ctp % For each pixel, go through and fit each 

            % Only fit if the voxel meets some criterion (after screening). For testing, don't do this.
            if 1
            % if vf_all.F0(ind) > 0.2
    
                % % Fitting the complex data all-in-one
                % % x0_neg = [p_neg.F0, p_neg.p0, p_neg.v_xgp0, p_neg.v_ygp0, p_neg.v_zgp0]; % ICs: [F0, p0, v_xgp0, v_ygp0, v_zgp0]
                % % x0_neg = [p_neg.p0(ind), p_neg.v_xgp0(ind), p_neg.v_ygp0(ind), p_neg.v_zgp0(ind)]; % ICs: [F0, p0, v_xgp0, v_ygp0, v_zgp0]
                % x0_neg = [1, 10e-3, 10e-3, vf_neg.v_zgp0(ind)]; % ICs: [p0, v_xgp0, v_ygp0, v_zgp0]
                % % lb = [0, 0, 0, 0];             % Lower bounds for parameters [SI units]
                % % ub = [1, 50e-3, 50e-3, 50e-3]; % Upper bounds for parameters [SI units]
                % f_temp = @(x, tau) g1vUS3D_vec(x, tau, sigma, vf_gen.k0, useF); % Use "anonymous function" to pass in the g1 vUS model function to the fitting
                % ydata = squeeze(g1neg_exp(ind, :));
                % x_neg = lsqcurvefit(f_temp, x0_neg, squeeze(tau(1:nTau)), ydata);
                % % x_neg = lsqcurvefit(f_temp, x0_neg, squeeze(tau(1:nTau)), squeeze(g1neg_exp(ind, :)), lb, ub);
                
                % Get the data and initial conditions for this voxel
                if useF
                    ydata = squeeze(g1_exp{j}(ind, 2:end)); ydata = ydata(:); % Experimental data to fit to

                    if useDC
                        x0 = [10e-3, vf(j).v_zgp0(ind), 1, vf(j).F0(ind), 0]; % ICs: [v_xgp0, v_zgp0, p0, F0, DC]
                    else
                        x0 = [10e-3, vf(j).v_zgp0(ind), 1, vf(j).F0(ind)]; % ICs: [v_xgp0, v_zgp0, p0, F0]
                    end
                else
                    ydata = squeeze(g1_exp{j}(ind, :)); ydata = ydata(:); % make sure it's a column vector
                    x0 = [10e-3, vf(j).v_zgp0(ind), 1]; % ICs: [v_xgp0, v_zgp0, p0]
    
                end
    
                % % Splitting the real and complex components
                % ydata_neg = ydata_neg(:);
                % ydata_neg_split = [real(ydata_neg), imag(ydata_neg)];
                % 
                % ydata_all = ydata_all(:);
                % ydata_all_split = [real(ydata_all), imag(ydata_all)];
    
                % % Magnitude only
                % ydata_neg = abs(ydata_neg(:));
                % ydata_neg_split = [real(ydata_neg), zeros(size(ydata_neg))];
                % 
                % ydata_all = abs(ydata_all(:));
                % ydata_all_split = [real(ydata_all), zeros(size(ydata_all))];
    
                % Perform the fitting
                % f_temp = @(x, tau) g1vUS3D_vec_split(x, tau, sigma, vf_gen.k0, useF); % Use "anonymous function" to pass in the g1 vUS model function to the fitting
                % % x_neg = lsqcurvefit(f_temp, x0_neg, squeeze(tau(1:nTau)), ydata_split);
                % % x_all = lsqcurvefit(f_temp, x0_all, squeeze(tau(1:nTau)), ydata_all_split, lb, ub, options);
                % if useF
                %     x_neg = lsqcurvefit(f_temp, x0_neg, squeeze(tau(2:nTau)), ydata_neg_split, lb, ub, options);
                %     F_neg(xi, yi, zi) = x_neg(5);
                % 
                %     x_all = lsqcurvefit(f_temp, x0_all, squeeze(tau(2:nTau)), ydata_all_split, lb, ub, options);
                %     F_all(xi, yi, zi) = x_all(5);
                % else
                %     x_neg = lsqcurvefit(f_temp, x0_neg, squeeze(tau(1:nTau)), ydata_neg_split, lb, ub, options);
                % 
                %     x_all = lsqcurvefit(f_temp, x0_all, squeeze(tau(1:nTau)), ydata_all_split, lb, ub, options);
                % end
                % 
                % % Store more results
                % vUS_neg(xi, yi, zi, :) = x_neg(2:4);
                % p_neg(xi, yi, zi) = x_neg(1);
                % 
                % vUS_all(xi, yi, zi, :) = x_all(2:4);
                % p_all(xi, yi, zi) = x_all(1);
    
    
                % ======== Fit only |g1| ========
                f_temp = @(x, tau) g1vUS2D_mag_vec(x, tau, sigma, vf_gen.k0, useF, useDC); % Use "anonymous function" to pass in the g1 vUS model function to the fitting
                if useF
                    x = lsqcurvefit(f_temp, x0, squeeze(tau(2:nTau)), abs(ydata), lb, ub, options);
                    F{j}(zi, xi) = x(4);

                    if useDC
                        DC{j}(zi, xi) = x(5);
                    end
                    
                else
                    x = lsqcurvefit(f_temp, x0, squeeze(tau(1:nTau)), abs(ydata), lb, ub, options);
                end
    
                % Store more results
                vUS{j}(zi, xi, :) = x(1:2);
                p{j}(zi, xi) = x(3);

            else % If we don't fit that pixel due to noisiness, set the velocities and other parameters to 0
                vUS{j}(zi, xi, :) = [0, 0];
                p{j}(zi, xi) = 0;
                if useF
                    F{j}(zi, xi) = 0;
                end
                if useDC
                    DC{j}(zi, xi) = 0;
                end
                
            end
        end
    end
end
toc

% % testg1 = g1vUS3D_vec(x_neg, tau, sigma, vf_gen.k0, useF);
% % figure; plot(ydata_neg, '-x', 'LineWidth', 2); hold on; plot(testg1, '-o', 'LineWidth', 1); hold off
% % figure; plot(abs(ydata_neg), '-x', 'LineWidth', 2); hold on; plot(abs(testg1), ':', 'LineWidth', 1); hold off
% % testspeed = sqrt(sum(x_neg(2:4).^2))
% 
% testg1 = g1vUS2D_vec(x_all, tau(2:end), sigma, vf_gen.k0, useF, useDC);
% figure; plot(ydata_all); hold on; plot(testg1); hold off
% figure; plot(abs(ydata_all)); hold on; plot(abs(testg1)); hold off
% testspeed = sqrt(sum(x_all(1:2).^2))
% 
% testg1 = g1vUS2D_vec(x, tau(2:end), sigma, vf_gen.k0, useF, useDC);
% figure; plot(ydata); hold on; plot(testg1); hold off
% figure; plot(abs(ydata)); hold on; plot(abs(testg1)); hold off
% testspeed = sqrt(sum(x(1:2).^2))
% 
% testg1 = g1vUS2D_vec(x_pos, tau(2:end), sigma, vf_gen.k0, useF, useDC);
% figure; plot(ydata_pos); hold on; plot(testg1); hold off
% figure; plot(abs(ydata_pos)); hold on; plot(abs(testg1)); hold off
% testspeed = sqrt(sum(x_pos(1:2).^2))

% testg1 = g1vUS3D_vec(x_all, tau(1:end), sigma, vf_gen.k0, useF);
% figure; plot(ydata_all); hold on; plot(testg1); hold off
% figure; plot(abs(ydata_all)); hold on; plot(abs(testg1)); hold off
% testspeed = sqrt(sum(x_all(2:4).^2))

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
testpt = [60, 133];
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