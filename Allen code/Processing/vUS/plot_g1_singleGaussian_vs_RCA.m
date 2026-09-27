%% Description:
%   Compare the g1 model for the single-Gaussian PSF (vUS_3D_quad_vec.m) against
%   the RCA-specific PSF (vUS_3D_quad_RCA_vec.m), for the SAME group velocity
%   components (v_xgp, v_ygp, v_zgp), and for several flow profiles (k, a). The
%   single-Gaussian model only knows the transverse speed v_tgp = sqrt(v_xgp^2 +
%   v_ygp^2); the RCA model also uses the transverse flow direction
%   phi = atan2(v_ygp, v_xgp).
%
%   Figure 1: |g1| vs. time lag; one row per velocity case, one column per flow profile.
%   Figure 2: complex g1 (real/imaginary parts, and the complex-plane trajectory) for
%       the velocity cases with axial flow; one column per flow profile.
%
%   Note: at a = 0 the profile is plug flow (f == 1), so k has no effect and
%   (k, a) = (2, 0) and (3, 0) are the same curve; they share one column.

clearvars

%% Add paths
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
AllenCodeDir = join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\";
addpath(genpath(fullfile(AllenCodeDir + "Processing\")))

%% Parameters
k0 = 2*pi / (1540/13.6e6);                          % [rad/m]
sigma_single = [130.1474, 130.1474, 44.7070]*1e-6;  % [sigma_x, sigma_y, sigma_z], field sigmas [m], from vUS_3D_newmodel.m (RC15gV, 5 x 2 angles from -6 to 6 deg)
sigma_RCA = [57.6, 286.5, 50.1]*1e-6;               % [sigma_narrow, sigma_wide, sigma_axial], field sigmas [m], from fit_RCA_PSF_v2.m (same data set)

F = 1; DC = 0;                                      % Purely dynamic signal

% Flow profiles: [k, a]. (k = 2, a = 1) is the "resulting model" on slides 24-25.
profiles = [2, 0; ...   % a = 0: plug flow (any k)
            2, 1; ...
            3, 1];
profile_labels = ["a = 0 (plug flow, any k)", "k = 2, a = 1", "k = 3, a = 1"];

tau = linspace(0, 20e-3, 401)';                     % Time lags [s]
[s, w] = gaussLegendre01(48);

% Group velocity cases: [v_xgp, v_ygp, v_zgp] [m/s]. Same transverse speed (15 mm/s) along an array axis and
% at 45 deg to it, each without and with axial flow. (For v_z >~ 8 mm/s the axial dephasing dominates g1 and
% the two PSF models become indistinguishable.)
vt_case = 15; vz_case = 3; % [mm/s]
vcases = [vt_case,          0,                 0; ...
          vt_case/sqrt(2),  vt_case/sqrt(2),   0; ...
          vt_case,          0,                 vz_case; ...
          vt_case/sqrt(2),  vt_case/sqrt(2),   vz_case] * 1e-3;

%% Evaluate both models for each velocity case and flow profile
nc = size(vcases, 1); np = size(profiles, 1);
g1_single = zeros(numel(tau), nc, np);
g1_RCA = zeros(numel(tau), nc, np);
case_labels = strings(1, nc);
for ci = 1:nc
    v_xgp = vcases(ci, 1); v_ygp = vcases(ci, 2); v_zgp = vcases(ci, 3);
    v_tgp = hypot(v_xgp, v_ygp);
    phi = atan2(v_ygp, v_xgp); % Transverse flow direction, measured from the x axis
    case_labels(ci) = sprintf('(%.1f, %.1f, %.1f) mm/s', v_xgp*1e3, v_ygp*1e3, v_zgp*1e3);

    for pi_ = 1:np
        k = profiles(pi_, 1); a = profiles(pi_, 2);
        g1_single(:, ci, pi_) = vUS_3D_quad_vec([v_tgp, v_zgp, F, DC, k, a], tau, k0, sigma_single, s, w);
        g1_RCA(:, ci, pi_) = vUS_3D_quad_RCA_vec([v_tgp, v_zgp, F, DC, k, a, phi], tau, k0, sigma_RCA, s, w);
    end
end

model_title = sprintf('Single Gaussian: \\sigma = [%.0f, %.0f, %.0f] \\mum. RCA: \\sigma_{narrow}, \\sigma_{wide}, \\sigma_{axial} = [%.1f, %.1f, %.1f] \\mum. Rows: (v_x, v_y, v_z)', ...
    sigma_single*1e6, sigma_RCA*1e6);

cSingle = [0, 0.447, 0.741];
cRCA = [0.85, 0.325, 0.098];
tms = tau*1e3; % [ms]

%% Figure 1: |g1|
figure('Name', '|g1|: single Gaussian vs. RCA PSF', 'Position', [50, 50, 1300, 1000]);
tl = tiledlayout(nc, np, 'TileSpacing', 'compact', 'Padding', 'compact');
for ci = 1:nc
    for pi_ = 1:np
        nexttile
        plot(tms, abs(g1_single(:, ci, pi_)), 'Color', cSingle, 'LineWidth', 2); hold on
        plot(tms, abs(g1_RCA(:, ci, pi_)), 'Color', cRCA, 'LineWidth', 2); hold off
        ylim([0, 1.05]); grid on
        if ci == 1, title(profile_labels(pi_)); end
        if ci == nc, xlabel('\tau [ms]'); end
        if pi_ == 1, ylabel({case_labels(ci), '|g_1|'}); end
        if ci == 1 && pi_ == 1, legend('Single-Gaussian PSF', 'RCA-specific PSF'); end
    end
end
title(tl, {'|g_1| for the same velocity components', model_title}, 'FontSize', 10);

%% Figure 2: complex g1, for the cases with axial flow
axial_cases = find(vcases(:, 3) ~= 0);
na = numel(axial_cases);
figure('Name', 'g1: single Gaussian vs. RCA PSF (axial flow)', 'Position', [50, 50, 1300, 1500]);
tl = tiledlayout(2*na, np, 'TileSpacing', 'compact', 'Padding', 'compact');
for ai = 1:na
    ci = axial_cases(ai);
    for pi_ = 1:np
        nexttile((2*ai - 2)*np + pi_) % Real and imaginary parts vs. tau
        plot(tms, real(g1_single(:, ci, pi_)), '-', 'Color', cSingle, 'LineWidth', 1.8); hold on
        plot(tms, imag(g1_single(:, ci, pi_)), '--', 'Color', cSingle, 'LineWidth', 1.8)
        plot(tms, real(g1_RCA(:, ci, pi_)), '-', 'Color', cRCA, 'LineWidth', 1.8)
        plot(tms, imag(g1_RCA(:, ci, pi_)), '--', 'Color', cRCA, 'LineWidth', 1.8); hold off
        ylim([-1.05, 1.05]); grid on
        if ai == 1, title(profile_labels(pi_)); end
        xlabel('\tau [ms]')
        if pi_ == 1, ylabel({case_labels(ci), 'g_1'}); end
        if ai == 1 && pi_ == 1, legend('Single, Re', 'Single, Im', 'RCA, Re', 'RCA, Im', 'Location', 'northeast'); end

        nexttile((2*ai - 1)*np + pi_) % Complex plane
        plot(real(g1_single(:, ci, pi_)), imag(g1_single(:, ci, pi_)), 'Color', cSingle, 'LineWidth', 1.8); hold on
        plot(real(g1_RCA(:, ci, pi_)), imag(g1_RCA(:, ci, pi_)), 'Color', cRCA, 'LineWidth', 1.8)
        plot(1, 0, 'ko', 'MarkerFaceColor', 'k', 'HandleVisibility', 'off'); hold off % tau = 0
        axis equal; xlim([-1.05, 1.05]); ylim([-1.05, 1.05]); grid on
        xlabel('Re(g_1)'); if pi_ == 1, ylabel('Im(g_1)'); end
        if ai == 1 && pi_ == 1, legend('Single-Gaussian PSF', 'RCA-specific PSF', 'Location', 'southwest'); end
    end
end
title(tl, {'Complex g_1 for the same velocity components (complex plane: dot = \tau = 0)', model_title}, 'FontSize', 10);
