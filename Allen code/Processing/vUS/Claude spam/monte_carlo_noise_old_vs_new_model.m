%% Description:
%   Same Monte Carlo design as monte_carlo_noise_vUS_2D_erf.m (SNR sweep
%   x 3 vx/vz ground-truth cases x 100 repeats), but fitting v_xgp with
%   BOTH models to compare recoverability:
%     - NEW model (Poiseuille/uniform-velocity-distribution): vUS_2D_erf,
%       3 free params [v_xgp, v_zgp, F] (+ DC)
%     - OLD model (Gaussian velocity distribution, Tang et al. 2020 Eq.
%       15): g1vUS2D, 4 free params [v_xgp, v_zgp, p, F] (+ DC), where p
%       is the sigma_vz = p*v_zgp linear factor
%
%   Each model is tested on data generated from ITS OWN formula (fair
%   like-for-like test of each model's intrinsic v_xgp identifiability,
%   not a test of which model matches reality better) and fit with its
%   own analytic-Jacobian residual function -- vUS_2D_erf_residJac.m for
%   the new model, g1vUS2D_residJac.m for the old model (both already
%   validated separately). Ground-truth v_xgp/v_zgp pairs, SNR levels,
%   and repeat count are identical between the two models for a direct,
%   paired comparison.
%
%   The OLD model's g1vUS2D_residJac.m header already documents a known
%   v_xgp/p identifiability problem (a full degenerate ridge, not just a
%   sign ambiguity) from a noise-free test -- this simulation quantifies
%   how that plays out under realistic noise, alongside the NEW model's
%   sign-ambiguity-only limitation.

%% Add path to the erfz code -- error function with complex inputs
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
addpath(genpath(ErrorFunctionCodePath))

outDir = fileparts(mfilename('fullpath'));

%% Fixed simulation parameters (matching monte_carlo_noise_vUS_2D_erf.m)
tau = (0:1/5000:20e-3).'; % s
tau_mask = 2:length(tau);
c0 = 1540; fc = 15.625e6; lambda0 = c0/fc; k0 = 2*pi/lambda0;
sigma = [58.9110, 73.9967] .* 1e-6; % [sigma_x, sigma_z] [m]

F_true = 0.85; DC_true = 0.05;
p_true = 0.3; % old model's sigma_vz = p*v_zgp linear factor, ground truth

cases = { ...
    [ 8e-3,  8e-3], 'balanced'; ...
    [15e-3,  2e-3], 'vx-dominant'; ...
    [ 2e-3, 15e-3], 'vz-dominant'};

SNR_dB_list = [30, 25, 20, 15, 10, 5];
numSamples = 100;

% NEW model: x = [v_xgp, v_zgp, F, DC]
lb_new = [0, -50e-3, 0, 0];
ub_new = [50e-3, 50e-3, 1, 1];
x0_new = [0.1e-3, 0.1e-3, 1, 0];

% OLD model: x = [v_xgp, v_zgp, p, F, DC]
lb_old = [0, -50e-3, 0, 0, 0];
ub_old = [50e-3, 50e-3, 1, 1, 1];
x0_old = [0.1e-3, 0.1e-3, 0.5, 1, 0]; % p0=0.5, a neutral (non-informative) midpoint guess

opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true);

nCases = size(cases, 1);
nSNR = numel(SNR_dB_list);

xFitNew = zeros(nCases, nSNR, numSamples, 4);
xFitOld = zeros(nCases, nSNR, numSamples, 5);

fprintf('Running %d fits per model (%d cases x %d SNR levels x %d repeats)...\n', ...
    nCases*nSNR*numSamples, nCases, nSNR, numSamples);
ticAll = tic;
for ci = 1:nCases
    v_xgp_true = cases{ci,1}(1); v_zgp_true = cases{ci,1}(2);

    g1_true_new = vUS_2D_erf(tau, k0, sigma, v_xgp_true, v_zgp_true, F_true, DC_true);
    g1_true_new_masked = g1_true_new(tau_mask);

    g1_true_old_split = g1vUS2D_vec_split([v_xgp_true, v_zgp_true, p_true, F_true, DC_true], ...
        tau(tau_mask), sigma, k0, true, true);
    g1_true_old_masked = complex(g1_true_old_split(:,1), g1_true_old_split(:,2));

    for si = 1:nSNR
        SNR_dB = SNR_dB_list(si);
        for ri = 1:numSamples
            % NEW model fit
            g1_noisy_new = awgn(g1_true_new_masked, SNR_dB);
            g1_noisy_new_split = [real(g1_noisy_new), imag(g1_noisy_new)];
            funNew = @(x) vUS_2D_erf_residJac(x, tau(tau_mask), k0, sigma, ...
                g1_noisy_new_split(:,1), g1_noisy_new_split(:,2));
            xFitNew(ci, si, ri, :) = lsqnonlin(funNew, x0_new, lb_new, ub_new, opts);

            % OLD model fit
            g1_noisy_old = awgn(g1_true_old_masked, SNR_dB);
            g1_noisy_old_split = [real(g1_noisy_old), imag(g1_noisy_old)];
            funOld = @(x) g1vUS2D_residJac(x, tau(tau_mask), sigma, k0, ...
                g1_noisy_old_split(:,1), g1_noisy_old_split(:,2));
            xFitOld(ci, si, ri, :) = lsqnonlin(funOld, x0_old, lb_old, ub_old, opts);
        end
    end
    fprintf('  case %d/%d (%s) done, %.1f s elapsed\n', ci, nCases, cases{ci,2}, toc(ticAll));
end
fprintf('Total time: %.1f s\n', toc(ticAll));

%% Bias / std / RMSE for v_xgp (column 1 in both models), per case x SNR
trueVals = zeros(nCases, 2); % [v_xgp, v_zgp] true, same for both models
for ci = 1:nCases
    trueVals(ci, :) = cases{ci,1};
end

biasNew = zeros(nCases, nSNR); sdNew = zeros(nCases, nSNR); rmseNew = zeros(nCases, nSNR);
biasOld = zeros(nCases, nSNR); sdOld = zeros(nCases, nSNR); rmseOld = zeros(nCases, nSNR);
for ci = 1:nCases
    for si = 1:nSNR
        sNew = squeeze(xFitNew(ci, si, :, 1)); % v_xgp samples, new model
        sOld = squeeze(xFitOld(ci, si, :, 1)); % v_xgp samples, old model

        biasNew(ci,si) = mean(sNew) - trueVals(ci,1);
        sdNew(ci,si) = std(sNew);
        rmseNew(ci,si) = sqrt(biasNew(ci,si)^2 + sdNew(ci,si)^2);

        biasOld(ci,si) = mean(sOld) - trueVals(ci,1);
        sdOld(ci,si) = std(sOld);
        rmseOld(ci,si) = sqrt(biasOld(ci,si)^2 + sdOld(ci,si)^2);
    end
end

%% Print summary table
fprintf('\n%-14s %6s   %12s %12s %12s   %12s %12s\n', ...
    'case', 'SNR', 'new rmse', 'old rmse', 'old/new', 'new bias', 'old bias');
for ci = 1:nCases
    for si = 1:nSNR
        fprintf('%-14s %5ddB   %10.3f  %10.3f  %10.2fx   %10.3f  %10.3f   [mm/s]\n', ...
            cases{ci,2}, SNR_dB_list(si), ...
            rmseNew(ci,si)*1e3, rmseOld(ci,si)*1e3, rmseOld(ci,si)/rmseNew(ci,si), ...
            biasNew(ci,si)*1e3, biasOld(ci,si)*1e3);
    end
end

%% Save results
save(fullfile(outDir, 'monte_carlo_old_vs_new_results.mat'), ...
    'cases', 'SNR_dB_list', 'numSamples', 'xFitNew', 'xFitOld', 'trueVals', ...
    'biasNew', 'sdNew', 'rmseNew', 'biasOld', 'sdOld', 'rmseOld', ...
    'p_true', 'F_true', 'DC_true');

%% Plot: v_xgp RMSE (top row) and bias (bottom row) vs SNR, old vs new, per case
fig = figure('Visible', 'off', 'Position', [100 100 1400 700]);
for ci = 1:nCases
    subplot(2, nCases, ci)
    plot(SNR_dB_list, rmseNew(ci,:)*1e3, '-o', 'LineWidth', 2, 'DisplayName', 'new model')
    hold on
    plot(SNR_dB_list, rmseOld(ci,:)*1e3, '-s', 'LineWidth', 2, 'DisplayName', 'old model')
    hold off
    set(gca, 'XDir', 'reverse'); grid on
    xlabel('SNR [dB]'); ylabel('v_{xgp} RMSE [mm/s]')
    title(sprintf('%s (v_{xgp}=%.0f, v_{zgp}=%.0f mm/s)', cases{ci,2}, cases{ci,1}(1)*1e3, cases{ci,1}(2)*1e3))
    legend('Location', 'northwest')

    subplot(2, nCases, nCases+ci)
    plot(SNR_dB_list, biasNew(ci,:)*1e3, '-o', 'LineWidth', 2, 'DisplayName', 'new model')
    hold on
    plot(SNR_dB_list, biasOld(ci,:)*1e3, '-s', 'LineWidth', 2, 'DisplayName', 'old model')
    yline(0, 'k--', 'HandleVisibility', 'off')
    hold off
    set(gca, 'XDir', 'reverse'); grid on
    xlabel('SNR [dB]'); ylabel('v_{xgp} bias [mm/s]')
    legend('Location', 'northwest')
end
sgtitle('v_{xgp} recoverability: new (uniform-distribution, 3 free params) vs. old (Gaussian-distribution, 4 free params incl. p) model')
exportgraphics(fig, fullfile(outDir, 'monte_carlo_old_vs_new_vxgp.png'), 'Resolution', 150)

fprintf('\nSaved:\n  %s\n  %s\n', ...
    fullfile(outDir, 'monte_carlo_old_vs_new_results.mat'), ...
    fullfile(outDir, 'monte_carlo_old_vs_new_vxgp.png'));
