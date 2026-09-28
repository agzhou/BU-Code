%% Description:
%   Builds the general (any-vessel-orientation) RCA g1 model into the fitting pipeline and
%   compares it against the two models it generalizes:
%     - vUS_3D_quad_vec.m         (single-Gaussian PSF, isotropic sigma)
%     - vUS_3D_quad_RCA_vec.m     (two-step RCA model: PSF assumed much wider than the vessel)
%     - vUS_3D_quad_RCA_vessel_vec.m (joint model: PSF weighting and flow profile combined)
%
%   Part 1: forward-model comparison. |g1(tau)| for all three models, at matched physical flow
%   parameters, across vessel radii from well below to well above the PSF width -- the "how much
%   it matters" table from RCA_g1s_derivation.md, now as curves.
%
%   Part 2: fitting demonstration. Simulates noisy g1(tau) data from the joint model at a fixed,
%   known vessel radius (200 um, comparable to sigma_wide) and known flow parameters, then fits the
%   SAME data three ways (single-Gaussian, two-step RCA, joint RCA+vessel with R fixed at its true
%   value) via lsqnonlin. Reports the recovered v_tgp, v_zgp against the ground truth for each.

clearvars
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
addpath(genpath(fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\Processing\")))

sigma = [57.6, 286.5, 50.1] * 1e-6;              % [sigma_narrow, sigma_wide, sigma_axial], fit_RCA_PSF_v2.m
sigma_iso = [130.1474, 130.1474, 44.7070] * 1e-6; % Per-axis single-Gaussian fit, vUS_3D_newmodel.m
k0 = 2*pi / (1540/13.6e6);
[s, w] = gaussLegendre01(48);
nth = 32;
frameRate = 2500; % [Hz]
tau = (0:40)'/frameRate;

%% ======================= Part 1: forward-model comparison ======================= %%
v_tgp = 8e-3; v_zgp = 5e-3; phi = deg2rad(20); F = 0.85; DC = 0.05; k = 2.5; a = 0.7;
x_iso6 = [v_tgp, v_zgp, F, DC, k, a];
x_twostep7 = [x_iso6, phi];

R_list = [30, 100, 200, 500, 1000] * 1e-6; % Vessel radius, well below to well above sigma_wide

g_iso = vUS_3D_quad_vec(x_iso6, tau, k0, sigma_iso, s, w); % Doesn't depend on R
g_twostep = vUS_3D_quad_RCA_vec(x_twostep7, tau, k0, sigma, s, w); % Doesn't depend on R either

figure('Name', 'RCA g1 model comparison', 'Position', [100, 100, 1400, 700]);
tiledlayout(2, 3, 'TileSpacing', 'compact');
for i = 1:numel(R_list)
    R = R_list(i);
    g_joint = vUS_3D_quad_RCA_vessel_vec([x_twostep7, R], tau, k0, sigma, s, w, nth);

    nexttile;
    plot(tau*1e3, abs(g_iso), ':', 'LineWidth', 1.5); hold on
    plot(tau*1e3, abs(g_twostep), '--', 'LineWidth', 1.5);
    plot(tau*1e3, abs(g_joint), '-', 'LineWidth', 2); hold off
    xlabel('\tau [ms]'); ylabel('|g_1(\tau)|'); ylim([0, 1]); grid on
    title(sprintf('R = %.0f \\mum (%.2f\\times\\sigma_{wide})', R*1e6, R/sigma(2)))
    if i == 1, legend('single-Gaussian', 'two-step RCA', 'joint RCA+vessel', 'Location', 'northeast'); end
end
nexttile; axis off
text(0, 0.9, sprintf(['Fixed physical parameters:\nv_{tgp} = %.0f mm/s, v_{zgp} = %.0f mm/s\n', ...
    '\\phi = %.0f deg, F = %.2f, DC = %.2f\nk = %.1f, a = %.1f\n\n', ...
    '\\sigma_{narrow} = %.1f, \\sigma_{wide} = %.1f, \\sigma_{axial} = %.1f \\mum\n\n', ...
    'Single-Gaussian and two-step RCA are\nR-independent by construction; only the\njoint model changes with R.'], ...
    v_tgp*1e3, v_zgp*1e3, rad2deg(phi), F, DC, k, a, sigma(1)*1e6, sigma(2)*1e6, sigma(3)*1e6), ...
    'VerticalAlignment', 'top', 'FontSize', 10)

fprintf('=== Part 1: forward comparison ===\n')
fprintf('%8s %10s %12s %12s %12s\n', 'R [um]', 'R/sw', '|g1| single', '|g1| 2-step', '|g1| joint')
tau_report = 16; % ms
[~, ir] = min(abs(tau*1e3 - tau_report));
for R = R_list
    g_joint_r = vUS_3D_quad_RCA_vessel_vec([x_twostep7, R], tau, k0, sigma, s, w, nth);
    fprintf('%8.0f %10.2f %12.4f %12.4f %12.4f  (at tau = %.1f ms)\n', ...
        R*1e6, R/sigma(2), abs(g_iso(ir)), abs(g_twostep(ir)), abs(g_joint_r(ir)), tau(ir)*1e3)
end

%% ======================= Part 2: fitting demonstration ======================= %%
fprintf('\n=== Part 2: fit noisy synthetic data (true R = 200 um) with all three models ===\n')
rng(7)
R_true = 200e-6;
x_true8 = [x_twostep7, R_true];
g_true = vUS_3D_quad_RCA_vessel_vec(x_true8, tau, k0, sigma, s, w, nth);

noise_sigma = 0.02; % Per-component (real, imag) additive Gaussian noise
g_noisy = g_true + noise_sigma*(randn(size(g_true)) + 1i*randn(size(g_true)));
yReal = real(g_noisy); yImag = imag(g_noisy);

opts = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);

% v_tgp is poorly identified from g1(tau) alone (a known issue -- see vUS_3D_newmodel.m's comments
% on v_tgp initialization) and the objective has multiple local minima in it, so each model gets a
% small multi-start over v_tgp0 and phi0 and keeps its best (lowest resnorm) result. v_zgp0 is set
% to its true value throughout, standing in for the independent, well-conditioned phase-based
% estimate (findVzPhaseDiff.m) the real pipeline uses to initialize it.
vtgp0_list = [2, 5, 10, 20] * 1e-3;
phi0_list = deg2rad([5, 20, 35]);

% ---- (a) Single-Gaussian isotropic model (6 params, no phi, no R) ----
lb6 = [0, -50e-3, 0, 0, 2, 0]; ub6 = [250e-3, 50e-3, 1, 1, 3, 1];
resid6 = @(x) reshape(vUS_3D_quad_vec(x, tau, k0, sigma_iso, s, w), [], 1) - complex(yReal, yImag);
resid6_split = @(x) [real(resid6(x)); imag(resid6(x))];
[x_fit6, resn6] = multiStartFit(resid6_split, vtgp0_list, [], v_zgp, [0.5, 0.1, 2.5, 0.8], lb6, ub6, opts);

% ---- (b) Two-step RCA model (7 params, phi in [0, 45 deg] by symmetry, no R) ----
lb7 = [lb6, 0]; ub7 = [ub6, pi/4];
resid7 = @(x) reshape(vUS_3D_quad_RCA_vec(x, tau, k0, sigma, s, w), [], 1) - complex(yReal, yImag);
resid7_split = @(x) [real(resid7(x)); imag(resid7(x))];
[x_fit7, resn7] = multiStartFit(resid7_split, vtgp0_list, phi0_list, v_zgp, [0.5, 0.1, 2.5, 0.8], lb7, ub7, opts);

% ---- (c) Joint RCA+vessel model (7 params, R fixed at its TRUE value) ----
resid_joint = @(x) vUS_3D_quad_RCA_vessel_residual(x, R_true, tau, k0, sigma, s, w, nth, yReal, yImag);
[x_fit_joint, resn_joint] = multiStartFit(resid_joint, vtgp0_list, phi0_list, v_zgp, [0.5, 0.1, 2.5, 0.8], lb7, ub7, opts);

fprintf('(resnorm at best start: single-Gaussian %.4f, two-step %.4f, joint %.4f)\n', resn6, resn7, resn_joint)

g_fit6 = vUS_3D_quad_vec(x_fit6, tau, k0, sigma_iso, s, w);
g_fit7 = vUS_3D_quad_RCA_vec(x_fit7, tau, k0, sigma, s, w);
g_fit_joint = vUS_3D_quad_RCA_vessel_vec([x_fit_joint, R_true], tau, k0, sigma, s, w, nth);

rmse = @(g) sqrt(mean(abs(g - g_noisy).^2));

fprintf('\nTrue:               v_tgp = %6.2f mm/s, v_zgp = %6.2f mm/s, phi = %5.1f deg\n', v_tgp*1e3, v_zgp*1e3, rad2deg(phi))
fprintf('Single-Gaussian fit: v_tgp = %6.2f mm/s, v_zgp = %6.2f mm/s  (%+5.1f%% / %+5.1f%%), RMSE = %.4f\n', ...
    x_fit6(1)*1e3, x_fit6(2)*1e3, 100*(x_fit6(1)/v_tgp-1), 100*(x_fit6(2)/v_zgp-1), rmse(g_fit6))
fprintf('Two-step RCA fit:    v_tgp = %6.2f mm/s, v_zgp = %6.2f mm/s, phi = %5.1f deg (%+5.1f%% / %+5.1f%%), RMSE = %.4f\n', ...
    x_fit7(1)*1e3, x_fit7(2)*1e3, rad2deg(x_fit7(7)), 100*(x_fit7(1)/v_tgp-1), 100*(x_fit7(2)/v_zgp-1), rmse(g_fit7))
fprintf('Joint RCA+vessel fit:v_tgp = %6.2f mm/s, v_zgp = %6.2f mm/s, phi = %5.1f deg (%+5.1f%% / %+5.1f%%), RMSE = %.4f\n', ...
    x_fit_joint(1)*1e3, x_fit_joint(2)*1e3, rad2deg(x_fit_joint(7)), 100*(x_fit_joint(1)/v_tgp-1), 100*(x_fit_joint(2)/v_zgp-1), rmse(g_fit_joint))

figure('Name', 'RCA g1 fit comparison'); tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; plot(tau*1e3, abs(g_noisy), '.', 'Color', [0.6 0.6 0.6]); hold on
plot(tau*1e3, abs(g_true), 'k-', 'LineWidth', 1);
plot(tau*1e3, abs(g_fit6), ':', 'LineWidth', 1.5);
plot(tau*1e3, abs(g_fit7), '--', 'LineWidth', 1.5);
plot(tau*1e3, abs(g_fit_joint), '-', 'LineWidth', 1.5); hold off
xlabel('\tau [ms]'); ylabel('|g_1(\tau)|'); grid on
legend('noisy data', 'true (noiseless)', 'single-Gaussian fit', 'two-step RCA fit', 'joint RCA+vessel fit')
title('Magnitude')

nexttile; plot(real(g_noisy), imag(g_noisy), '.', 'Color', [0.6 0.6 0.6]); hold on
plot(real(g_true), imag(g_true), 'k-', 'LineWidth', 1);
plot(real(g_fit6), imag(g_fit6), ':', 'LineWidth', 1.5);
plot(real(g_fit7), imag(g_fit7), '--', 'LineWidth', 1.5);
plot(real(g_fit_joint), imag(g_fit_joint), '-', 'LineWidth', 1.5); hold off
axis equal; xlim([-0.2, 1]); ylim([-0.6, 0.6]); grid on
xlabel('Re(g_1)'); ylabel('Im(g_1)'); title('Complex plane')

%% Helper: multi-start lsqnonlin over v_tgp0 (and phi0, if given), keeping the best (lowest resnorm)
% fit. x0_rest is [F0, DC0, k0, a0] (shared across starts); v_zgp0 is fixed across starts, standing
% in for an independent, well-conditioned estimate (see the comment above).
function [x_best, resn_best] = multiStartFit(residFun, vtgp0_list, phi0_list, vzgp0, x0_rest, lb, ub, opts)
    if isempty(phi0_list), phi0_list = 0; end % 6-param model: no phi, loop runs once with a dummy value
    nParam = numel(lb);
    resn_best = Inf; x_best = [];
    for vtgp0 = vtgp0_list
        for phi0 = phi0_list
            if nParam == 6
                x0 = [vtgp0, vzgp0, x0_rest];
            else
                x0 = [vtgp0, vzgp0, x0_rest, phi0];
            end
            x0 = min(max(x0, lb), ub); % Clip into bounds
            [x, resn] = lsqnonlin(residFun, x0, lb, ub, opts);
            if resn < resn_best
                resn_best = resn; x_best = x;
            end
        end
    end
end
