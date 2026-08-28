%% Add path to the erfz code -- error function with complex inputs
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
addpath(genpath(ErrorFunctionCodePath))

%% Parameters for testing
v_max = 0.03; % m/s
tau = 0:1/5000:20e-3; % s
c0 = 1540; % m/s
fc = 15.625e6; % Hz
lambda0 = c0./fc; % m
k0 = 2*pi/lambda0; % m^-1
sigma = 75e-6; % m

%% Calculate g1 with the new model vs. old
[g1_erf] = vUS_1D_erf(tau, k0, sigma, v_max);

v_gp = v_max/2; % According to the uniform distribution, the mean velocity would be half v_max
v_min = 0;
variance_v = (v_max - v_min)^2 / 12;
sigma_v = sqrt(variance_v);
F = 1; [g1_old] = g1vUS1D(F, v_gp, sigma_v, sigma, k0, tau);

%% Compare results
% |g1|
figure
plot(tau.*1e3, abs(g1_erf), '-o', 'LineWidth', 2)
hold on
plot(tau.*1e3, abs(g1_old), '-x', 'LineWidth', 2)
hold off
xlabel('Time lag [ms]')
ylabel('|g1|')
legend('Uniform velocity probability distribution', 'Gaussian velocity probability distribution')

% Full complex g1
figure
plot(g1_erf, '-o', 'LineWidth', 2)
hold on
plot(g1_old, '-x', 'LineWidth', 2)
hold off
axis square
legend('Uniform velocity probability distribution', 'Gaussian velocity probability distribution')