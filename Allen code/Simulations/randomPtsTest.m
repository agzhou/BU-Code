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

SP.vesselDiam = 100e-6; % Vessel diameter [m]
SP.vesselLength = SP.endDepthMM/1e3;  % Vessel length [m]

% Define the center of the vessel  
SP.xstart = 0;
SP.ystart = 0;
SP.zstart = 0;

% Get points in a cylindrical "vessel"
% cyl_vessel = genRandomPts3D_cyl(vesselDiam, vesselLength, startDepthMM/1e3, xstart, ystart, zstart);
[cyl_vessel, SP] = genRandomPts3D_cyl(SP);
figure; scatter3(cyl_vessel(:, 1), cyl_vessel(:, 2), cyl_vessel(:, 3), '.'); axis square



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define a rotation matrix for final manipulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SP.dim = 3;

SP.flow_v_mm_s = 30;
% flow_v_mm_s = 3000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test to see if replacement works
SP.flow_dim = 3
% new_cyl_vessel = movePoints(cyl_vessel, dim, flow_v_mm_s, frameRate, vesselDiam, startDepthMM, endDepthMM, xstart, ystart, zstart);
new_cyl_vessel = movePoints(cyl_vessel, SP);
figure; scatter3(new_cyl_vessel(:, 1), new_cyl_vessel(:, 2), new_cyl_vessel(:, 3), '.'); axis square

%%