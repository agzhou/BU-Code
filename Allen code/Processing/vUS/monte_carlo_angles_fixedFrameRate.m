%% Description:
%   Variant of monte_carlo_angles_vs_frames.m that decouples frame rate
%   from angle count N: frame rate is held CONSTANT at 5000 Hz across all
%   N, with the system PRF assumed able to go as high as needed
%   (PRF = frameRate*N, reported for reference only, not a constraint).
%
%   This isolates the per-frame-SNR (~sqrt(N) compounding gain) channel
%   from the temporal-sampling/frame-count/aliasing channel that dominated
%   monte_carlo_angles_vs_frames.m's results -- with frame rate, tau
%   sampling density, total frame count, and v_zgp Nyquist margin all now
%   IDENTICAL across N, any remaining trend in v_xgp/v_zgp recoverability
%   vs. N must come from per-frame SNR alone. Answers: "if temporal
%   sampling were free (unlimited PRF), would more compounding just
%   straightforwardly help?"
%
%   Same synthesis/estimation/fitting methodology as
%   monte_carlo_angles_vs_frames.m (see that file and
%   synthesizeCorrelatedG1Process.m for the validated methodology and its
%   caveats -- SNR0 anchor is an assumption, sigma_x held fixed).

%% Add path to the erfz code
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
addpath(genpath(ErrorFunctionCodePath))
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Speckle tracking')); % g1T.m

outDir = fileparts(mfilename('fullpath'));

%% Fixed system/acquisition parameters
c0 = 1540; fc = 15.625e6; lambda0 = c0/fc; k0 = 2*pi/lambda0;
sigma = [58.9110, 73.9967] .* 1e-6; % held fixed across N, as before

frameRate = 5000; % [Hz] CONSTANT across all N -- the change from before
T_acq = 0.200; % [s] total acquisition time, fixed
nFrames = round(T_acq * frameRate); % now also constant across N

F_true = 0.85; DC_true = 0.05;
SNR0_dB = 8; % per-frame SNR at the reference angle count N_ref (same assumption as before)
N_ref = 5;

N_list = [1, 3, 5, 9, 17];
numSamples = 60;

cases = { ...
    8e-3, 8e-3, 'moderate flow'; ...
    8e-3, 45e-3, 'fast flow'};

lb = [0, -100e-3, 0, 0];
ub = [50e-3, 100e-3, 1, 1];
x0 = [0.1e-3, 0.1e-3, 1, 0];
opts = optimoptions('lsqnonlin', 'Display', 'off', 'SpecifyObjectiveGradient', true);

nCases = size(cases, 1);
nN = numel(N_list);

SNR_dB_list = SNR0_dB + 5*log10(N_list ./ N_ref); % sqrt(N) amplitude-SNR gain, same model as before
PRF_implied = frameRate .* N_list; % reported for reference only -- not a constraint per this run's assumption

nyquist_vzmax_mm_s = frameRate * pi / (2*k0) * 1e3;
fprintf('Frame rate fixed at %d Hz for all N -> Nyquist v_zgp limit = %.1f mm/s for every condition\n', frameRate, nyquist_vzmax_mm_s);
fprintf('N       PRF_implied[Hz]  SNR[dB]\n');
for i = 1:nN
    fprintf('%-7d %-16.0f %-.2f\n', N_list(i), PRF_implied(i), SNR_dB_list(i));
end

xFit = zeros(nCases, nN, numSamples, 4);

fprintf('\nRunning %d total fits (%d cases x %d N-values x %d repeats)...\n', ...
    nCases*nN*numSamples, nCases, nN, numSamples);
ticAll = tic;
for ci = 1:nCases
    v_xgp_true = cases{ci,1}; v_zgp_true = cases{ci,2};
    g1fun = @(tau) vUS_2D_erf(max(tau,eps), k0, sigma, v_xgp_true, v_zgp_true, F_true, DC_true) .* (tau~=0) + (tau==0);

    tau_full = (0:nFrames-1)'/frameRate;
    nTau = min(round(0.02*frameRate), nFrames-1); % same 20 ms decay window convention, now constant across N
    tau_fit = tau_full(2:nTau+1);

    for ni = 1:nN
        noiseVar = 10^(-SNR_dB_list(ni)/10);

        for ri = 1:numSamples
            s = synthesizeCorrelatedG1Process(g1fun, frameRate, nFrames);
            s = s + sqrt(noiseVar/2)*(randn(nFrames,1) + 1i*randn(nFrames,1));

            g1_est = g1T(reshape(s, 1, []), nTau+1);
            g1_est = g1_est(2:nTau+1);
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
fprintf('\n%-16s %4s   %14s %14s   %14s %14s\n', 'case', 'N', 'vxgp rmse(mm/s)', 'vxgp bias(mm/s)', 'vzgp rmse(mm/s)', 'vzgp bias(mm/s)');
for ci = 1:nCases
    for ni = 1:nN
        fprintf('%-16s %4d   %14.3f %14.3f   %14.3f %14.3f\n', cases{ci,3}, N_list(ni), ...
            rmse(ci,ni,1)*1e3, bias(ci,ni,1)*1e3, rmse(ci,ni,2)*1e3, bias(ci,ni,2)*1e3);
    end
end

save(fullfile(outDir, 'monte_carlo_angles_fixedFrameRate_results.mat'), ...
    'cases', 'N_list', 'frameRate', 'nFrames', 'SNR_dB_list', 'PRF_implied', 'numSamples', ...
    'xFit', 'trueVals', 'bias', 'sd', 'rmse', 'T_acq');

%% Plot
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
subplot(2,2,1); hold off; grid on; xlabel('angles per frame (N)'); ylabel('v_{xgp} RMSE [mm/s]'); title('v_{xgp} RMSE vs. compounding'); legend('Location','best')
subplot(2,2,2); hold off; grid on; xlabel('angles per frame (N)'); ylabel('v_{zgp} RMSE [mm/s]'); title('v_{zgp} RMSE vs. compounding'); legend('Location','best')
subplot(2,2,3); hold off; grid on; xlabel('angles per frame (N)'); ylabel('v_{xgp} bias [mm/s]'); title('v_{xgp} bias vs. compounding'); yline(0,'k--','HandleVisibility','off'); legend('Location','best')
subplot(2,2,4); hold off; grid on; xlabel('angles per frame (N)'); ylabel('v_{zgp} bias [mm/s]'); title('v_{zgp} bias vs. compounding'); yline(0,'k--','HandleVisibility','off'); legend('Location','best')
sgtitle(sprintf('Frame rate FIXED at %d Hz for all N (PRF assumed unconstrained) -- isolates the sqrt(N) SNR-only effect', frameRate))
exportgraphics(fig, fullfile(outDir, 'monte_carlo_angles_fixedFrameRate.png'), 'Resolution', 150)

fprintf('\nSaved:\n  %s\n  %s\n', fullfile(outDir,'monte_carlo_angles_fixedFrameRate_results.mat'), fullfile(outDir,'monte_carlo_angles_fixedFrameRate.png'));
