endDepthMM = 6;
startDepth = 0;
wl = 1540 ./ 13.8889 ./ 1e6;

vesselX = 100e-6;    % x dimension
vesselY = 100e-6;    % y dimension
vesselZ = endDepthMM/1e3;  % z dimension
%  
xstart = 0;
ystart = 0;
zstart = 0;
test = randomPts3D_func(vesselX, vesselY, vesselZ, wl, startDepth, xstart, ystart, zstart);

 figure; scatter3(test(:, 1), test(:, 2), test(:, 3), '.')