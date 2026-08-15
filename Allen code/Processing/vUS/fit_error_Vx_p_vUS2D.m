%% Parameters
sigma = [60, Inf, 75] .* 1e-6; % sigma x, y, z [m]
c0 = 1540; % Speed of sound [m/s]
fc = 15.625e6; % Frequency [Hz]
lambda0 = c0/fc; % Wavelength [m]
k0 = 2*pi/lambda0; % Wavenumber [m^-1]

%% Create grids
v_xgp_bounds = [0, 30].*1e-3; % Bounds for v_xgp [m/s]
v_zgp_bounds = [0, 30].*1e-3; % Bounds for v_zgp [m/s]
vstep = 1*1e-3; % Interval for the velocity meshes/grids [m/s]
v_xgp_vec = v_xgp_bounds(1):vstep:v_xgp_bounds(2);
v_zgp_vec = v_zgp_bounds(1):vstep:v_zgp_bounds(2);
[Vx, Vz] = meshgrid(v_xgp_vec, v_zgp_vec);

%% Create function to calculate v_xgp error, if we don't fit the p term separately
error_v_xgp_func = @(v_xgp, v_zgp, p) sqrt(v_xgp.^2 + (4*pi*p*sigma(1)/lambda0 .*v_zgp).^2) - v_xgp; % v_xgp and v_zgp can be matrices, p is a single value

%% Calculate errors
p = 1;
error_v_xgp = error_v_xgp_func(Vx, Vz, p);
% figure; imagesc(error_v_xgp)
relative_error = error_v_xgp ./ Vx; % Relative error: error/v_xgp_actual

%% Plots
fs = 12;
% Plot absolute error
figure; surf(Vx.*1e3, Vz.*1e3, error_v_xgp.*1e3); xlabel('True v_{xgp} [mm/s]'); ylabel('True v_{zgp} [mm/s]'); zlabel('Absolute error in fitted v_{xgp} [mm/s]'); title("p = " + num2str(p))
fontsize(fs, 'points')

% Plot relative error
% figure; surf(Vx.*1e3, Vz.*1e3, relative_error .* 100); xlabel('True v_xgp [mm/s]'); ylabel('True v_{zgp} [mm/s]'); zlabel('Relative error in fitted v_xgp [%]')
figure; surf(Vx.*1e3, Vz.*1e3, log10(relative_error)); xlabel('True v_{xgp} [mm/s]'); ylabel('True v_{zgp} [mm/s]'); zlabel('log10(relative error in fitted v_{xgp})'); title("p = " + num2str(p))
fontsize(fs, 'points')
