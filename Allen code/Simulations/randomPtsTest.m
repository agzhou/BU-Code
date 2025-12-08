clearvars
%% Define Simulation Parameter struct
SP.endDepthMM = 1; % End depth [mm]
SP.startDepthMM = 0; % Start depth [mm]
SP.wl = 1540 ./ 13.8889 ./ 1e6; % Wavelength [m]
SP.frameRate = 2500; % Frame rate [Hz]
% vesselX = 100e-6;    % x dimension
% vesselY = 100e-6;    % y dimension
% vesselZ = endDepthMM/1e3;  % z dimension
SP.scatterReflectivity = 1.0;
SP.sigma = [300e-6, 300e-6, 150e-6]; %%%% PSF testing %%%%


SP.vesselDiam = 50e-6; % Vessel diameter [m]
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

% SP.flow_v_mm_s = 30;
SP.flow_v_mm_s = 125;
SP.flow_dim = 3; %%%%%%%%
% new_cyl_vessel = movePoints(cyl_vessel, dim, flow_v_mm_s, frameRate, vesselDiam, startDepthMM, endDepthMM, xstart, ystart, zstart);
[test_new_cyl_vessel, test_SP] = movePoints(cyl_vessel, SP);
plotPoints(test_new_cyl_vessel, test_SP)

%% Define a voxel
voxel.center = [0, 0, 0]; % Center coords of the voxel
voxel.size = [100e-6, 100e-6, 100e-6]; % Define x, y, z dimensions of the voxel

% Define time steps
SP.numFrames = 100;

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

    [new_cyl_vessel, SP] = movePoints(new_cyl_vessel, SP);
end

%% Plot for testing
plotPoints(new_cyl_vessel, SP)
figure; plot(abs(voxel.sIQ))
voxel.g1 = sim_g1T(voxel.sIQ);
figure; plot(abs(voxel.g1))
figure; plot(real(voxel.g1), imag(voxel.g1), '-o')
% figure; scatter3(voxel.data(:, 1), voxel.data(:, 2), voxel.data(:, 3), '.'); axis square


%%%% something is wrong with moving the points - there are many being added
%%%% over time, but the total # should be more or less constant