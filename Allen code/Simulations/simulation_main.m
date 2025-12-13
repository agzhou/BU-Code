clearvars
%% Define Simulation Parameter struct
% SP.endDepthMM = 1; % End depth [mm]
% SP.startDepthMM = 0; % Start depth [mm]
SP.c = 1540; % Speed of sound [m/s]
SP.f = 13.8889 * 1e6; % Ultrasound frequency [Hz]
SP.wl = SP.c / SP.f; % Wavelength [m]
SP.frameRate = 2500; % Frame rate [Hz]
% vesselX = 100e-6;    % x dimension
% vesselY = 100e-6;    % y dimension
% vesselZ = endDepthMM/1e3;  % z dimension
SP.scatterReflectivity = 1.0;
SP.sigma = [300e-6, 300e-6, 150e-6]; %%%% PSF testing %%%%

SP.snr = 50; % Choose the SNR for the data vs. Gaussian white noise (5 is what Bingxue used)

% SP.vesselDiam = 50e-6; % Vessel diameter [m]
SP.vesselDiam = 100e-6; % Vessel diameter [m]

% SP.vesselLength = (SP.endDepthMM - SP.startDepthMM)/1e3;  % Vessel length [m]
SP.vesselLength = 8 * 1e-3;  % Vessel length [m]

% Define the center of the vessel  
SP.xstart = 0;
SP.ystart = 0;
SP.zstart = 0;

% Get points in a cylindrical "vessel"
% cyl_vessel = genRandomPts3D_cyl(vesselDiam, vesselLength, startDepthMM/1e3, xstart, ystart, zstart);
[cyl_vessel, SP] = genRandomPts3D_cyl(SP);
% plotPoints(cyl_vessel, SP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define a rotation matrix for final manipulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SP.dim = 3;

SP.flow_v_mm_s = 30;
% SP.flow_v_mm_s = 125;
% SP.flow_v_mm_s = 1250;
SP.flow_dim = 3; %%%%%%%%
% new_cyl_vessel = movePoints(cyl_vessel, dim, flow_v_mm_s, frameRate, vesselDiam, startDepthMM, endDepthMM, xstart, ystart, zstart);
% [test_new_cyl_vessel, test_SP] = movePoints(cyl_vessel, SP);
% plotPoints(test_new_cyl_vessel, test_SP)

%% Define a voxel
voxel.center = [0, 0, 0]; % Center coords of the voxel
voxel.size = [100e-6, 100e-6, 100e-6]; % Define x, y, z dimensions of the voxel

% Define time steps
% SP.numFrames = 50;
SP.numFrames = 500;

% Get all the data within the voxel at frame 1
voxel.data = getDataInVoxel(cyl_vessel, voxel); % Note: voxel.data for now is just a container that is always changing
voxel.sIQ(1) = voxel_sIQ(voxel, SP);

[new_cyl_vessel, SP] = movePoints(cyl_vessel, SP);
% plotPoints(new_cyl_vessel, SP)

% Go through each frame, moving the points, and update the voxel data/sIQ
for fi = 2:SP.numFrames
    voxel.data = getDataInVoxel(new_cyl_vessel, voxel);
    voxel.sIQ(fi) = voxel_sIQ(voxel, SP);
    % temp_sIQ = voxel_sIQ(voxel, SP);
    % voxel.sIQ(fi) = temp_sIQ - mean(temp_sIQ);
    % voxel.sIQ(fi) = temp_sIQ;

    % plotPoints(voxel.data, SP)
    % plotPoints(new_cyl_vessel, SP)

    [new_cyl_vessel, SP] = movePoints(new_cyl_vessel, SP); % Update points after moving
end
clearvars fi temp_sIQ
voxel.sIQ = voxel.sIQ - mean(voxel.sIQ); % ZERO MEAN THE sIQ

%% Plot for testing
tau = 0:1/SP.frameRate:(SP.numFrames-1)/SP.frameRate;

% plotPoints(new_cyl_vessel, SP)
% figure; plot(tau, abs(voxel.sIQ))
% figure; plot(tau, abs(voxel.sIQ - mean(voxel.sIQ)))
% figure; plot(tau, real(voxel.sIQ - mean(voxel.sIQ)))
% figure; plot(tau, imag(voxel.sIQ - mean(voxel.sIQ)))
voxel.g1 = sim_g1T(voxel.sIQ);
% voxel.g1 = sim_g1T(voxel.sIQ - mean(voxel.sIQ));
figure; plot(tau .* 1e3, abs(voxel.g1)); xlabel('tau [ms]'); ylabel("|g_1|")
figure; plot(real(voxel.g1), imag(voxel.g1), '-o')
% figure; scatter3(voxel.data(:, 1), voxel.data(:, 2), voxel.data(:, 3), '.'); axis square


%% FFT to visualize the flow speed's effect
fD = -2 * SP.f * (SP.flow_v_mm_s/1e3)/SP.c;
F = fftshift(fft(voxel.sIQ));
f = linspace(-SP.frameRate/2, SP.frameRate/2, length(F));
figure; plot(f, abs(F)); xlabel('f [Hz]'); hold on
xline(abs(fD), 'r-', 'LineWidth', 2)

%% Nonlinear Least Squares fitting of |g1T|

% Define a maximum tau to fit to (improve accuracy)
tau_max = 50 / 1e3; % [s]
tau_mask = tau < tau_max;
tau_range = tau(tau_mask);

% Define a function handle for |g1T|
% Let x = [Ns, vx, vy, vz]
Rs = SP.scatterReflectivity;
% Rs = mean(voxel.data(:, 4)); % Average reflection coefficient of scatterer
% Re = SP.snr; % "Noise level of imaging system"
Re = 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re = 1e-8
abs_g1T_model = @(x, tau) x(1)*Rs*(pi^1.5)*SP.sigma(1)*SP.sigma(2)*SP.sigma(3) ...
                          / ( x(1)*Rs*(pi^1.5)*SP.sigma(1)*SP.sigma(2)*SP.sigma(3) + Re ) ...
                          .* exp( -(x(2) .* tau).^2./(4.*SP.sigma(1)^2) -(x(3) .* tau).^2./(4.*SP.sigma(2)^2) -(x(4) .* tau).^2./(4.*SP.sigma(3)^2) );


% Fit the unknown parameters (x) with the simulation and compare to the ground truth
% TESTING: INITIAL GUESS AND UPPER/LOWER BOUNDS
testx = [size(voxel.data, 1), 0, 0, SP.flow_v_mm_s/1e3];
x0 = testx %%%%%%%% testing
lb = [0, 0, 0, 0]
ub = [Inf, 0, 0, 1000e-3]
% x_fit = lsqcurvefit(abs_g1T_model, x0, tau, abs(voxel.g1))
% x_fit = lsqcurvefit(abs_g1T_model, x0, tau, abs(voxel.g1), lb, ub)
x_fit = lsqcurvefit(abs_g1T_model, x0, tau_range, abs(voxel.g1(tau_mask)), lb, ub)

%% Input the ground truth parameters to see what the |g1T| model looks like

test_abs_g1T_model = abs_g1T_model(testx, tau_range);
figure; plot(tau_range*1e3, test_abs_g1T_model); xlabel('tau [ms]'); ylabel("|g_1|")
hold on; plot(tau_range .* 1e3, abs(voxel.g1(tau_mask)));

% Plot the fit
abs_g1T_fit = abs_g1T_model(x_fit, tau_range);
plot(tau_range .* 1e3, abs_g1T_fit);
% test2 = abs_g1T_model([7849, 0, 0, .0445], tau);
% test2 = abs_g1T_model([7849, 0, 0, .0845], tau);
% plot(tau .* 1e3, test2);
hold off
legend("|g_1| model with ground truth parameters", "|g_1| simulation", "Fit")
%%

% %% Nonlinear Least Squares fitting of FULL g1T
% % Define a function handle for g1T
% % Let x = [Ns, vx, vy, vz]
% Rs = SP.scatterReflectivity;
% % Rs = mean(voxel.data(:, 4)); % Average reflection coefficient of scatterer
% % Re = SP.snr; % "Noise level of imaging system"
% Re = 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Re = 1e-8
% k = 2*pi/SP.wl;
% g1T_model = @(x, tau) x(1)*Rs*(pi^1.5)*SP.sigma(1)*SP.sigma(2)*SP.sigma(3) ...
%                           / ( x(1)*Rs*(pi^1.5)*SP.sigma(1)*SP.sigma(2)*SP.sigma(3) + Re ) ...
%                           .* exp( -(x(2) .* tau).^2./(4.*SP.sigma(1)^2) -(x(3) .* tau).^2./(4.*SP.sigma(2)^2) -(x(4) .* tau).^2./(4.*SP.sigma(3)^2) ...
%                           .* exp( 2.*1i.*(k).*x(4).*tau ));
% 
% 
% % Fit the unknown parameters (x) with the simulation and compare to the ground truth
% x0 = testx %%%%%%%% testing
% x_fit = lsqcurvefit(g1T_model, x0, tau, voxel.g1)
% 
% %% Input the ground truth parameters to see what the g1T model looks like
% testx = [size(voxel.data, 1), 0, 0, SP.flow_v_mm_s/1e3];
% test_abs_g1T_model = g1T_model(testx, tau);
% figure; plot(tau*1e3, test_abs_g1T_model); xlabel('tau [ms]'); ylabel("|g_1|")
% hold on; plot(tau .* 1e3, abs(voxel.g1));
% legend("g_1 model with ground truth parameters", "g_1 simulation")
% % Plot the fit
% abs_g1T_fit = g1T_model(x_fit, tau);
% plot(tau .* 1e3, abs_g1T_fit);
% % test2 = abs_g1T_model([7849, 0, 0, .0445], tau);
% test2 = g1T_model([7849, 0, 0, .0845], tau);
% plot(tau .* 1e3, test2);
% hold off
