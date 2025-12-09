clearvars
%% Define Simulation Parameter struct
SP.endDepthMM = 5; % End depth [mm]
SP.startDepthMM = 0; % Start depth [mm]
SP.wl = 1540 ./ 13.8889 ./ 1e6; % Wavelength [m]
SP.frameRate = 2500; % Frame rate [Hz]
% vesselX = 100e-6;    % x dimension
% vesselY = 100e-6;    % y dimension
% vesselZ = endDepthMM/1e3;  % z dimension
SP.scatterReflectivity = 1.0;
SP.sigma = [300e-6, 300e-6, 150e-6]; %%%% PSF testing %%%%


% SP.vesselDiam = 50e-6; % Vessel diameter [m]
SP.vesselDiam = 100e-6; % Vessel diameter [m]

SP.vesselLength = (SP.endDepthMM - SP.startDepthMM)/1e3;  % Vessel length [m]

% Define the center of the vessel  
SP.xstart = 0;
SP.ystart = 0;
SP.zstart = 0;

% Get points in a cylindrical "vessel"
% cyl_vessel = genRandomPts3D_cyl(vesselDiam, vesselLength, startDepthMM/1e3, xstart, ystart, zstart);
[cyl_vessel, SP] = genRandomPts3D_cyl(SP);
plotPoints(cyl_vessel, SP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define a rotation matrix for final manipulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SP.dim = 3;

SP.flow_v_mm_s = 30;
% SP.flow_v_mm_s = 125;
SP.flow_dim = 3; %%%%%%%%
% new_cyl_vessel = movePoints(cyl_vessel, dim, flow_v_mm_s, frameRate, vesselDiam, startDepthMM, endDepthMM, xstart, ystart, zstart);
[test_new_cyl_vessel, test_SP] = movePoints(cyl_vessel, SP);
plotPoints(test_new_cyl_vessel, test_SP)

%% Define a voxel
voxel.center = [0, 0, 0]; % Center coords of the voxel
voxel.size = [100e-6, 100e-6, 100e-6]; % Define x, y, z dimensions of the voxel

% Define time steps
SP.numFrames = 50;

% Get all the data within the voxel at frame 1
voxel.data = getDataInVoxel(cyl_vessel, voxel); % Note: voxel.data for now is just a container that is always changing
voxel.sIQ(1) = voxel_sIQ(voxel, SP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% need to add noise to the sIQ %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

new_cyl_vessel = movePoints(cyl_vessel, SP);
% Go through each frame, moving the points, and update the voxel data/sIQ
for fi = 2:SP.numFrames
    voxel.data = getDataInVoxel(new_cyl_vessel, voxel);
    voxel.sIQ(fi) = voxel_sIQ(voxel, SP);

    % plotPoints(voxel.data, SP)
    plotPoints(new_cyl_vessel, SP)

    [new_cyl_vessel, SP] = movePoints(new_cyl_vessel, SP); % Update points after moving
end

%% Plot for testing
tau = 0:1/SP.frameRate:(SP.numFrames-1)/SP.frameRate;

plotPoints(new_cyl_vessel, SP)
figure; plot(tau, abs(voxel.sIQ))
voxel.g1 = sim_g1T(voxel.sIQ);
figure; plot(tau, abs(voxel.g1))
figure; plot(real(voxel.g1), imag(voxel.g1), '-o')
% figure; scatter3(voxel.data(:, 1), voxel.data(:, 2), voxel.data(:, 3), '.'); axis square



% test = autocorr(abs(voxel.sIQ));
% test = sim_g1T(abs(voxel.sIQ));
% figure; plot(abs(test))
test = autocorr(abs(voxel.g1), NumLags=length(voxel.g1)-1);
figure; plot(tau.*1e3, test, '-o'); xlabel('Tau [ms]')


%% Test
% t = 0:0.1:2*pi * 4;
% y1 = sin(t);
% % y2 = cos(t);
% figure;
% hold on
% % plot(t, y1, t, y2);
% % plot(t, y1);
% plot(y1);
% 
% 
% % y1 = ones(size(t));
% 
% test = sim_g1T(y1);
% autocorr_test = autocorr(y1, "NumLags", length(y1) - 1);
% % plot(t, test)
% plot(test)
% plot(autocorr_test)
% hold off