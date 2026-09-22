%% Add path to the erfz code -- error function with complex inputs
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
addpath(genpath(ErrorFunctionCodePath))

%% Parameters for testing
v_xgp = 0.01; % m/s
v_ygp = 0.01; % m/s
v_zgp = 0.02; % m/s
tau = 0:1/5000:20e-3; % s
c0 = 1540; % m/s
fc = 15.625e6; % Hz
lambda0 = c0./fc; % m
k0 = 2*pi/lambda0; % m^-1
sigma = [429, 429, 126].*1e-6; % 1/e PSF values [m] for the RC15gV probe at 13.6 MHz and 11 x 2 angles from -5 to 5 deg (G:\My Drive\Data\RC15gV PSF sim - 11 angles from -5 to 5 deg)
% sigma = [100, 100, 100].*1e-6;

%% Calculate g1 with the new model vs. old
[g1_erf] = vUS_3D_erf(tau, k0, sigma, v_xgp, v_ygp, v_zgp);

% v_gp = v_max/2; % According to the uniform distribution, the mean velocity would be half v_max
% v_min = 0;
% variance_v = (v_max - v_min)^2 / 12;
% sigma_v = sqrt(variance_v);
% F = 1; [g1_old] = g1vUS1D(F, v_gp, sigma_v, sigma, k0, tau);

%% Compare results
% |g1|
figure
plot(tau.*1e3, abs(g1_erf), '-o', 'LineWidth', 2)
hold on
% plot(tau.*1e3, abs(g1_old), '-x', 'LineWidth', 2)
hold off
xlabel('Time lag [ms]')
ylabel('|g1|')
legend('Uniform velocity probability distribution', 'Gaussian velocity probability distribution')

% Full complex g1
figure
plot(g1_erf, '-o', 'LineWidth', 2)
hold on
% plot(g1_old, '-x', 'LineWidth', 2)
hold off
axis square
% legend('Uniform velocity probability distribution', 'Gaussian velocity probability distribution')