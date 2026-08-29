%% Add path to the erfz code -- error function with complex inputs
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
addpath(genpath(ErrorFunctionCodePath))

%% Parameters for testing
v_xgp = 0.01; % m/s
v_zgp = 0.02; % m/s
tau = 0:1/5000:20e-3; % s
c0 = 1540; % m/s
fc = 15.625e6; % Hz
lambda0 = c0./fc; % m
k0 = 2*pi/lambda0; % m^-1
sigma = [58.9110, 73.9967].*1e-6; % Field-based 1/e PSF values (x, z) [m] for the L22-14v probe at 15.625 MHz and 17 angles from -10 to 10 deg. The y component is set to some arbitrary positive number but it won't really be used. (G:\My Drive\Data\PSF Simulations\L22-14v PSF sim - 17 angles from -10 to 10 deg)

%% Calculate g1 with the new model vs. old
[g1_erf] = vUS_2D_erf(tau, k0, sigma, v_xgp, v_zgp);

v_zmax = v_zgp*2; % According to the uniform distribution, the mean velocity would be half v_max
v_min = 0;
variance_v = (v_max - v_min)^2 / 12;
sigma_v = sqrt(variance_v);
F = 1; [g1_old] = g1vUS2D(F, v_xgp, v_zgp, p, sigma, k0, tau);

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
legend('Uniform velocity probability distribution', 'Gaussian velocity probability distribution')