% Generate a cylindrical amount of points
% Inputs:
%   vesselDiam: diameter of the vessel [m]
%   vesselLength = length of the vessel [m]

% Output is in meters

function pts = genRandomPts3D_cyl(vesselDiam, vesselLength, startDepth, xstart, ystart, zstart)
    % x, y, start are measured with x, y = 0 at center of transducer
    % z start is distance starting from startDepth
    rng('shuffle')

    % Start by defining a rectangular vessel
    vesselWidth = vesselDiam;

    cellDensity = 1 / 100 / 1e-18; % 1/(100 um^3) from Bingxue's paper
%     cellDensity = 1 / 1e-12; %%%%%%%%%%%% test %%%%%%%%%%%%
    % cellDensity = 1 / 1e-15; %%%%%%%%%%%% test %%%%%%%%%%%%
    numCells = uint16(vesselWidth * vesselWidth * vesselLength * cellDensity)
    
    % MIGHT want randn so the distribution of cells is in the Gaussian profile?
    xr = rand([1, numCells])';
    yr = rand([1, numCells])';
    zr = rand([1, numCells])';
    scatterReflectivity = 1.0;
    
%     for n = 1:numCells
    %     Media.MP(n, :) = [ -vesselWidth / wl / 2 + xr(n) * vesselWidth / wl, 0, 50 + zr(n) * vesselHeight / wl, scatterReflectivity];
        pts = [ xstart - vesselWidth ./ 2 + xr .* vesselWidth, ...
                ystart - vesselWidth ./ 2 + yr .* vesselWidth, ...
                zstart + startDepth + zr .* vesselLength, ...
                repmat(scatterReflectivity, length(xr), 1)];

%         pts = [ -vesselX / wl / 2 + xr(n) * vesselX / wl, 0, 0 + zr(n) * vesselZ / wl, scatterReflectivity];
        
        % figure; scatter3(pts(:, 1), pts(:, 2), pts(:, 3), '.'); axis square

%     end
    mask = sqrt(pts(:, 1).^2 + pts(:, 2).^2) > vesselDiam/2; % Remove points outside of the vessel radius
    pts(mask, :) = [];
return