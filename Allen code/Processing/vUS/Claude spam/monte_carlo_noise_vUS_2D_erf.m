%% Description:
%   Monte Carlo test of how accurately the new (Poiseuille / uniform-
%   velocity-distribution) g1 model's parameters -- especially v_xgp vs.
%   v_zgp -- can be recovered from noisy g1(tau) data, across a sweep of
%   SNR regimes. Extends fitting_test_vUS_2D_erf.m's "Sensitivity to SNR
%   analysis" section (single SNR, single ground truth) to a full
%   SNR x ground-truth-case sweep with bias/std/RMSE summary statistics.
%
%   Three ground-truth cases are tested at each SNR level to show whether
%   recoverability depends on the vx/vz balance -- expected from theory
%   (see project discussion) since v_xgp only enters the model through
%   the amplitude-decay envelope (no phase information), while v_zgp
%   carries phase information too:
%     - balanced:      v_xgp and v_zgp comparable
%     - vx-dominant:   mostly transverse flow
%     - vz-dominant:   mostly axial flow
%
%   Uses the analytic-Jacobian residual (vUS_2D_erf_residJac.m) for
%   speed, given the large number of total fits (cases x SNR levels x
%   repeats). Bounds/x0 are FIXED and wide (not proportional to a
%   per-voxel initial guess), unlike vUS_2D_newmodel.m's per-voxel mesh-
%   based bounds -- this isolates the model/optimization's inherent
%   recoverability from the bounds-collapse and other implementation
%   issues flagged separately, so it answers "can the true velocities be
%   recovered at all, given enough optimizer freedom" cleanly.

%% Add path to the erfz code -- error function with complex inputs
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
addpath(genpath(ErrorFunctionCodePath))

outDir = fileparts(mfilename('fullpath'));

%% Fixed simulation parameters (matching fitting_test_vUS_2D_erf.m)
tau = (0:1/5000:20e-3).'; % s
tau_mask = 2:length(tau); % skip tau=0 (NaN in the model)
c0 = 1540; % m/s
fc = 15.625e6; % Hz
lambda0 = c0/fc; % m
k0 = 2*pi/lambda0; % rad/m
sigma = [58.9110, 73.9967] .* 1e-6; % [sigma_x, sigma_z], L22-14v PSF [m]

F_true = 0.85;   % dynamic fraction ground truth
DC_true = 0.05;  % static offset ground truth

%% Ground-truth cases: [v_xgp, v_zgp] in m/s
cases = { ...
    [ 8e-3,  8e-3], 'balanced'; ...
    [15e-3,  2e-3], 'vx-dominant'; ...
    [ 2e-3, 15e-3], 'vz-dominant'};

SNR_dB_list = [30, 25, 20, 15, 10, 5];
numSamples = 100; % Monte Carlo repeats per (case, SNR)

lb = [0, -50e-3, 0, 0];
ub = [50e-3, 50e-3, 1, 1];
x0 = [0.1e-3, 0.1e-3, 1, 0]; % deliberately generic/uninformed initial guess

opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true);

nCases = size(cases, 1);
nSNR = numel(SNR_dB_list);

% Results: [nCases, nSNR, numSamples, 4] fitted parameters
xFit = zeros(nCases, nSNR, numSamples, 4);

fprintf('Running %d total fits (%d cases x %d SNR levels x %d repeats)...\n', ...
    nCases*nSNR*numSamples, nCases, nSNR, numSamples);
ticAll = tic;
for ci = 1:nCases
    v_xgp_true = cases{ci,1}(1); v_zgp_true = cases{ci,1}(2);
    g1_true = vUS_2D_erf(tau, k0, sigma, v_xgp_true, v_zgp_true, F_true, DC_true);
    g1_true_masked = g1_true(tau_mask);

    for si = 1:nSNR
        SNR_dB = SNR_dB_list(si);
        for ri = 1:numSamples
            g1_noisy = awgn(g1_true_masked, SNR_dB);
            g1_noisy_split = [real(g1_noisy), imag(g1_noisy)];

            fun = @(x) vUS_2D_erf_residJac(x, tau(tau_mask), k0, sigma, ...
                g1_noisy_split(:,1), g1_noisy_split(:,2));

            xFit(ci, si, ri, :) = lsqnonlin(fun, x0, lb, ub, opts);
        end
    end
    fprintf('  case %d/%d (%s) done, %.1f s elapsed\n', ci, nCases, cases{ci,2}, toc(ticAll));
end
fprintf('Total time: %.1f s\n', toc(ticAll));

%% Bias / std / RMSE summary, per case x SNR x parameter
paramNames = {'v_xgp', 'v_zgp', 'F', 'DC'};
trueVals = zeros(nCases, 4);
for ci = 1:nCases
    trueVals(ci, :) = [cases{ci,1}(1), cases{ci,1}(2), F_true, DC_true];
end

bias = zeros(nCases, nSNR, 4);
sd   = zeros(nCases, nSNR, 4);
rmse = zeros(nCases, nSNR, 4);
for ci = 1:nCases
    for si = 1:nSNR
        samples = squeeze(xFit(ci, si, :, :)); % [numSamples, 4]
        bias(ci, si, :) = mean(samples, 1) - trueVals(ci, :);
        sd(ci, si, :) = std(samples, 0, 1);
        rmse(ci, si, :) = sqrt(bias(ci, si, :).^2 + sd(ci, si, :).^2);
    end
end

%% Print summary table (v_xgp and v_zgp, in mm/s)
fprintf('\n%-14s %6s   %14s %14s %14s   %14s %14s %14s\n', ...
    'case', 'SNR', 'vxgp bias(mm/s)', 'vxgp sd(mm/s)', 'vxgp rmse(mm/s)', ...
    'vzgp bias(mm/s)', 'vzgp sd(mm/s)', 'vzgp rmse(mm/s)');
for ci = 1:nCases
    for si = 1:nSNR
        fprintf('%-14s %5ddB   %16.3f %14.3f %14.3f   %14.3f %14.3f %14.3f\n', ...
            cases{ci,2}, SNR_dB_list(si), ...
            bias(ci,si,1)*1e3, sd(ci,si,1)*1e3, rmse(ci,si,1)*1e3, ...
            bias(ci,si,2)*1e3, sd(ci,si,2)*1e3, rmse(ci,si,2)*1e3);
    end
end

%% Save numeric results for later inspection
save(fullfile(outDir, 'monte_carlo_noise_vUS_2D_erf_results.mat'), ...
    'cases', 'SNR_dB_list', 'numSamples', 'xFit', 'trueVals', 'bias', 'sd', 'rmse', ...
    'tau', 'sigma', 'k0', 'F_true', 'DC_true', 'paramNames');

%% Plot 1: RMSE of v_xgp and v_zgp vs SNR, one line per case
fig1 = figure('Visible', 'off', 'Position', [100 100 1000 420]);
colors = lines(nCases);
subplot(1,2,1)
hold on
for ci = 1:nCases
    plot(SNR_dB_list, squeeze(rmse(ci,:,1))*1e3, '-o', 'LineWidth', 2, 'Color', colors(ci,:))
end
hold off
set(gca, 'XDir', 'reverse')
xlabel('SNR [dB]'); ylabel('v_{xgp} RMSE [mm/s]')
title('v_{xgp} recovery error vs. noise')
legend(cases(:,2), 'Location', 'northwest')
grid on

subplot(1,2,2)
hold on
for ci = 1:nCases
    plot(SNR_dB_list, squeeze(rmse(ci,:,2))*1e3, '-o', 'LineWidth', 2, 'Color', colors(ci,:))
end
hold off
set(gca, 'XDir', 'reverse')
xlabel('SNR [dB]'); ylabel('v_{zgp} RMSE [mm/s]')
title('v_{zgp} recovery error vs. noise')
legend(cases(:,2), 'Location', 'northwest')
grid on

sgtitle('Recovery accuracy vs. SNR: v_{xgp} (no phase info) vs. v_{zgp} (has phase info)')
exportgraphics(fig1, fullfile(outDir, 'monte_carlo_rmse_vs_snr.png'), 'Resolution', 150)

%% Plot 2: bias vs SNR (signed), to show systematic bias trends
fig2 = figure('Visible', 'off', 'Position', [100 100 1000 420]);
subplot(1,2,1)
hold on
for ci = 1:nCases
    plot(SNR_dB_list, squeeze(bias(ci,:,1))*1e3, '-o', 'LineWidth', 2, 'Color', colors(ci,:))
end
yline(0, 'k--')
hold off
set(gca, 'XDir', 'reverse')
xlabel('SNR [dB]'); ylabel('v_{xgp} bias [mm/s]')
title('v_{xgp} bias vs. noise')
legend(cases(:,2), 'Location', 'northwest')
grid on

subplot(1,2,2)
hold on
for ci = 1:nCases
    plot(SNR_dB_list, squeeze(bias(ci,:,2))*1e3, '-o', 'LineWidth', 2, 'Color', colors(ci,:))
end
yline(0, 'k--')
hold off
set(gca, 'XDir', 'reverse')
xlabel('SNR [dB]'); ylabel('v_{zgp} bias [mm/s]')
title('v_{zgp} bias vs. noise')
legend(cases(:,2), 'Location', 'northwest')
grid on

sgtitle('Systematic bias vs. SNR (v_{xgp} is bounded >= 0, so noise pushes it up, not down)')
exportgraphics(fig2, fullfile(outDir, 'monte_carlo_bias_vs_snr.png'), 'Resolution', 150)

fprintf('\nSaved:\n  %s\n  %s\n  %s\n', ...
    fullfile(outDir, 'monte_carlo_noise_vUS_2D_erf_results.mat'), ...
    fullfile(outDir, 'monte_carlo_rmse_vs_snr.png'), ...
    fullfile(outDir, 'monte_carlo_bias_vs_snr.png'));
