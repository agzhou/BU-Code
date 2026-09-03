%% Description:
%   Monte Carlo test of the acquisition trade-off: fewer frames with more
%   compounded angles each, vs. more frames with fewer compounded angles
%   each, at a FIXED total acquisition time T -- which is more beneficial
%   for v_xgp/v_zgp recoverability?
%
%   Unlike the noise-regime Monte Carlo scripts (which add noise directly
%   to the clean g1(tau) curve), this simulates realistic per-frame data:
%   synthesizeCorrelatedG1Process.m generates a complex time series with
%   the target autocorrelation, additive per-frame complex Gaussian noise
%   is added (amplitude SNR improving as sqrt(N) with angle count, the
%   standard coherent-compounding result), and g1(tau) is ESTIMATED from
%   that noisy series via g1T.m -- the same estimator used on real data.
%   This is necessary to capture the mechanisms actually in question:
%     - frame rate = PRF/N sets tau sampling density within the decay
%       window (fewer angles -> higher frame rate -> finer tau sampling)
%     - total frames = T*frameRate sets how many frame-PAIRS g1T.m
%       averages at each lag (more frames -> lower estimator variance)
%     - per-frame SNR improves ~sqrt(N) with angle count (more angles ->
%       better per-frame image quality)
%     - v_zgp phase-aliasing risk (Nyquist: 2*k0*v_zgp/pi < frameRate)
%       grows as frame rate drops
%
%   PRF = 30 kHz (this system's max single-angle plane-wave rate, per
%   Tang et al. 2020). T = 200 ms (their actual acquisition window).
%   Per-frame SNR is anchored at SNR0 dB at the paper's actual N=5 choice
%   -- this anchor point is a free assumption (not independently
%   measured), clearly labeled as such; the qualitative N-dependence
%   trends are the actual result of interest, not the absolute numbers.
%
%   Two ground-truth cases: "moderate flow" (well within Nyquist at every
%   N tested) and "fast flow" (v_zgp chosen to alias at the largest N /
%   lowest frame rate), to separately show the tau-resolution/estimator-
%   noise trade-off and the hard aliasing failure mode.

%% Add path to the erfz code
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
addpath(genpath(ErrorFunctionCodePath))
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Speckle tracking')); % g1T.m

outDir = fileparts(mfilename('fullpath'));

%% Fixed system/acquisition parameters
c0 = 1540; fc = 15.625e6; lambda0 = c0/fc; k0 = 2*pi/lambda0;
sigma = [58.9110, 73.9967] .* 1e-6; % [sigma_x, sigma_z] [m] -- held FIXED across N
                                      % (isolating the temporal-sampling/
                                      % estimator-noise/per-frame-SNR
                                      % mechanisms from any possible
                                      % angle-count effect on lateral
                                      % resolution, which is a separate,
                                      % not-independently-characterized
                                      % effect -- see caveats)
PRF = 30000; % [Hz] max single-angle plane-wave rate
T_acq = 0.200; % [s] total acquisition time, fixed across all N

F_true = 0.85; DC_true = 0.05;
SNR0_dB = 8; % per-frame SNR at the reference angle count N_ref
N_ref = 5;   % Tang et al. 2020's actual choice

N_list = [1, 3, 5, 9, 17]; % angle counts to sweep
numSamples = 60; % Monte Carlo repeats per (case, N) -- fewer than the
                  % other MC scripts since per-repeat cost is higher here
                  % (full time-series synthesis, not just a clean-curve
                  % noise draw)

cases = { ...
    8e-3, 8e-3, 'moderate flow'; ...
    8e-3, 45e-3, 'fast flow (aliases at low frame rate)'};

lb = [0, -100e-3, 0, 0];
ub = [50e-3, 100e-3, 1, 1];
x0 = [0.1e-3, 0.1e-3, 1, 0];
opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true);

nCases = size(cases, 1);
nN = numel(N_list);

frameRate_list = PRF ./ N_list;
nFrames_list = round(T_acq .* frameRate_list);
SNR_dB_list = SNR0_dB + 5*log10(N_list ./ N_ref); % sqrt(N) amplitude-SNR gain

fprintf('N       frameRate[Hz]  nFrames  SNR[dB]  Nyquist_vzmax[mm/s]\n');
for i = 1:nN
    nyquist_vzmax = frameRate_list(i) * pi / (2*k0) * 1e3;
    fprintf('%-7d %-14.0f %-8d %-8.2f %-.1f\n', N_list(i), frameRate_list(i), nFrames_list(i), SNR_dB_list(i), nyquist_vzmax);
end

xFit = zeros(nCases, nN, numSamples, 4);

fprintf('\nRunning %d total fits (%d cases x %d N-values x %d repeats)...\n', ...
    nCases*nN*numSamples, nCases, nN, numSamples);
ticAll = tic;
for ci = 1:nCases
    v_xgp_true = cases{ci,1}; v_zgp_true = cases{ci,2};
    g1fun = @(tau) vUS_2D_erf(max(tau,eps), k0, sigma, v_xgp_true, v_zgp_true, F_true, DC_true) .* (tau~=0) + (tau==0);

    for ni = 1:nN
        frameRate = frameRate_list(ni);
        nFrames = nFrames_list(ni);
        noiseVar = 10^(-SNR_dB_list(ni)/10); % E[|s|^2] ~= 1, so this is the noise power directly

        tau_full = (0:nFrames-1)'/frameRate;
        nTau = min(round(0.02*frameRate), nFrames-1); % fit within a 20 ms decay window, same convention as vUS_2D_newmodel.m
        tau_fit = tau_full(2:nTau+1); % skip tau=0

        for ri = 1:numSamples
            s = synthesizeCorrelatedG1Process(g1fun, frameRate, nFrames);
            s = s + sqrt(noiseVar/2)*(randn(nFrames,1) + 1i*randn(nFrames,1)); % additive per-frame complex noise

            g1_est = g1T(reshape(s, 1, []), nTau+1);
            g1_est = g1_est(2:nTau+1); % drop tau=0
            g1_est_split = [real(g1_est(:)), imag(g1_est(:))];

            fun = @(x) vUS_2D_erf_residJac(x, tau_fit, k0, sigma, g1_est_split(:,1), g1_est_split(:,2));
            xFit(ci, ni, ri, :) = lsqnonlin(fun, x0, lb, ub, opts);
        end
    end
    fprintf('  case %d/%d (%s) done, %.1f s elapsed\n', ci, nCases, cases{ci,3}, toc(ticAll));
end
fprintf('Total time: %.1f s\n', toc(ticAll));

%% Bias / std / RMSE for v_xgp and v_zgp
trueVals = cell2mat(cases(:,1:2));
bias = zeros(nCases, nN, 2); sd = zeros(nCases, nN, 2); rmse = zeros(nCases, nN, 2);
for ci = 1:nCases
    for ni = 1:nN
        for pi_ = 1:2
            samples = squeeze(xFit(ci, ni, :, pi_));
            bias(ci,ni,pi_) = mean(samples) - trueVals(ci,pi_);
            sd(ci,ni,pi_) = std(samples);
            rmse(ci,ni,pi_) = sqrt(bias(ci,ni,pi_)^2 + sd(ci,ni,pi_)^2);
        end
    end
end

%% Print summary
fprintf('\n%-30s %4s   %14s %14s   %14s %14s\n', 'case', 'N', 'vxgp rmse(mm/s)', 'vxgp bias(mm/s)', 'vzgp rmse(mm/s)', 'vzgp bias(mm/s)');
for ci = 1:nCases
    for ni = 1:nN
        fprintf('%-30s %4d   %14.3f %14.3f   %14.3f %14.3f\n', cases{ci,3}, N_list(ni), ...
            rmse(ci,ni,1)*1e3, bias(ci,ni,1)*1e3, rmse(ci,ni,2)*1e3, bias(ci,ni,2)*1e3);
    end
end

save(fullfile(outDir, 'monte_carlo_angles_vs_frames_results.mat'), ...
    'cases', 'N_list', 'frameRate_list', 'nFrames_list', 'SNR_dB_list', 'numSamples', ...
    'xFit', 'trueVals', 'bias', 'sd', 'rmse', 'PRF', 'T_acq');

%% Plot: v_xgp and v_zgp RMSE vs. N (angle count), one line per case
fig = figure('Visible', 'off', 'Position', [100 100 1300 700]);
colors = lines(nCases);
for ci = 1:nCases
    subplot(2,2,1); hold on
    plot(N_list, squeeze(rmse(ci,:,1))*1e3, '-o', 'LineWidth', 2, 'Color', colors(ci,:), 'DisplayName', cases{ci,3})
    subplot(2,2,2); hold on
    plot(N_list, squeeze(rmse(ci,:,2))*1e3, '-o', 'LineWidth', 2, 'Color', colors(ci,:), 'DisplayName', cases{ci,3})
    subplot(2,2,3); hold on
    plot(N_list, squeeze(bias(ci,:,1))*1e3, '-o', 'LineWidth', 2, 'Color', colors(ci,:), 'DisplayName', cases{ci,3})
    subplot(2,2,4); hold on
    plot(N_list, squeeze(bias(ci,:,2))*1e3, '-o', 'LineWidth', 2, 'Color', colors(ci,:), 'DisplayName', cases{ci,3})
end
subplot(2,2,1); hold off; grid on; xlabel('angles per frame (N)'); ylabel('v_{xgp} RMSE [mm/s]'); title('v_{xgp} RMSE vs. compounding'); legend('Location','northwest')
subplot(2,2,2); hold off; grid on; xlabel('angles per frame (N)'); ylabel('v_{zgp} RMSE [mm/s]'); title('v_{zgp} RMSE vs. compounding'); legend('Location','northwest')
subplot(2,2,3); hold off; grid on; xlabel('angles per frame (N)'); ylabel('v_{xgp} bias [mm/s]'); title('v_{xgp} bias vs. compounding'); yline(0,'k--','HandleVisibility','off'); legend('Location','best')
subplot(2,2,4); hold off; grid on; xlabel('angles per frame (N)'); ylabel('v_{zgp} bias [mm/s]'); title('v_{zgp} bias vs. compounding'); yline(0,'k--','HandleVisibility','off'); legend('Location','best')
sgtitle(sprintf('Fixed total acquisition time (%.0f ms): fewer angles/more frames (left) vs. more angles/fewer frames (right)', T_acq*1e3))
exportgraphics(fig, fullfile(outDir, 'monte_carlo_angles_vs_frames.png'), 'Resolution', 150)

fprintf('\nSaved:\n  %s\n  %s\n', fullfile(outDir,'monte_carlo_angles_vs_frames_results.mat'), fullfile(outDir,'monte_carlo_angles_vs_frames.png'));
