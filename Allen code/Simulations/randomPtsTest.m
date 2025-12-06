clearvars
%%
% endDepthMM = 6;
endDepthMM = 1; % End depth [mm]
startDepthMM = 0; % Start depth [mm]
wl = 1540 ./ 13.8889 ./ 1e6; % Wavelength [m]
frameRate = 2500; % Frame rate [Hz]
% vesselX = 100e-6;    % x dimension
% vesselY = 100e-6;    % y dimension
% vesselZ = endDepthMM/1e3;  % z dimension

vesselDiam = 100e-6; % Vessel diameter [m]
vesselLength = endDepthMM/1e3;  % Vessel length [m]

% Define the center of the vessel  
xstart = 0;
ystart = 0;
zstart = 0;

% Get points in a cylindrical "vessel"
cyl_vessel = genRandomPts3D_cyl(vesselDiam, vesselLength, startDepthMM/1e3, xstart, ystart, zstart);
figure; scatter3(cyl_vessel(:, 1), cyl_vessel(:, 2), cyl_vessel(:, 3), '.'); axis square



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define a rotation matrix for final manipulation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    flow_v_mm_s = 30;
    % flow_v_mm_s = 3000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test to see if replacement works
    dim = 3
new_cyl_vessel = movePoints(cyl_vessel, dim, flow_v_mm_s, frameRate, vesselDiam, startDepthMM, endDepthMM, xstart, ystart, zstart);
figure; scatter3(new_cyl_vessel(:, 1), new_cyl_vessel(:, 2), new_cyl_vessel(:, 3), '.'); axis square

% Move all the params into a struct..................