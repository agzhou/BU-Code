%% Description: test for how effectively my new g1 model can be fit

%% Add path to the erfz code -- error function with complex inputs
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
% AllenVerasonicsCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "Allen code"))), '\') + "\Verasonics");
ErrorFunctionCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "BU-Code"))), '\') + "\Allen Code\ErrorFunction\");
addpath(genpath(ErrorFunctionCodePath))

%% Parameters for testing
v_xgp = 0.005; % m/s
v_zgp = 0.01; % m/s
tau = (0:1/5000:20e-3).'; % s
c0 = 1540; % m/s
fc = 15.625e6; % Hz
lambda0 = c0./fc; % m
k0 = 2*pi/lambda0; % m^-1
sigma = [58.9110, 73.9967].*1e-6; % Field-based 1/e PSF values (x, z) [m] for the L22-14v probe at 15.625 MHz and 17 angles from -10 to 10 deg. The y component is set to some arbitrary positive number but it won't really be used. (G:\My Drive\Data\PSF Simulations\L22-14v PSF sim - 17 angles from -10 to 10 deg)

%% Calculate g1 with the new model vs. old
[g1_erf] = vUS_2D_erf(tau, k0, sigma, v_xgp, v_zgp);

% v_zmax = v_zgp*2; % According to the uniform distribution, the mean velocity would be half v_max
% v_min = 0;
variance_f = (2 - 0)^2 / 12;
sigma_Vz = v_zgp * sqrt(variance_f);
F = 1; [g1_old] = g1vUS2D_sigmaVz(F, v_xgp, v_zgp, sigma_Vz, sigma, k0, tau);

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

%% Fitting test - can we recover the ground truth parameters
tau_mask = 2:length(tau); % Don't fit at tau = 0, which results in a NaN
% g1_exp = g1_erf(tau_mask);
% SNR_dB = 30;
SNR_dB = 50;
g1_exp = awgn(g1_erf(tau_mask), SNR_dB);
g1_exp_split = [real(g1_exp), imag(g1_exp)]; % Can add noise if desired
fun_new = @(x) vUS_2D_erf_vec_split(x, tau(tau_mask), k0, sigma) - g1_exp_split;

lb = [0, -50e-3, 0];
ub = [50e-3, 50e-3, 1];
x0 = [0.0001, 0.0001, 1]; % v_xgp0, v_zgp0, F guesses [m/s]

x = lsqnonlin(fun_new, x0, lb, ub);

%% Plot fitting results
% |g1|
g1_fitted_erf = vUS_2D_erf_vec(x, tau, k0, sigma);
figure
plot(tau(tau_mask).*1e3, abs(g1_exp), '-o', 'LineWidth', 2)
hold on
plot(tau.*1e3, abs(g1_fitted_erf), '-x', 'LineWidth', 2)
hold off
xlabel('Time lag [ms]')
ylabel('|g1|')
legend("Noisy model data: v_{xgp} = "+num2str(v_xgp*1e3)+" mm/s, v_{zgp} = "+num2str(v_zgp*1e3) + " mm/s", "Fit: v_{xgp} = "+num2str(x(1)*1e3)+"mm/s, v_{zgp} = "+num2str(x(2)*1e3) + " mm/s")
title("Signal with Gaussian white noise added; SNR: " + num2str(SNR_dB) + " dB")