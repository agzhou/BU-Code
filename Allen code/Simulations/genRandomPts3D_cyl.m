% Generate a cylindrical amount of points
% Inputs: struct SP with the parameters:
%   vesselDiam: diameter of the vessel [m]
%   vesselLength = length of the vessel [m]

% Output is in meters

% function pts = genRandomPts3D_cyl(vesselDiam, vesselLength, startDepth, xstart, ystart, zstart)
function [pts, SP] = genRandomPts3D_cyl(SP)
    % x, y, start are measured with x, y = 0 at center of transducer
    % z start is distance starting from startDepth
    rng('shuffle')

    % Start by defining a rectangular vessel
    SP.vesselWidth = SP.vesselDiam;

    % if ~exist('SP.')
    % SP.cellDensity = 1 / 100 / 1e-18; % 1/(100 um^3) from Bingxue's paper; here it's defined as [cells/m^3]
    % SP.cellDensity = 1 / 1e-12; %%%%%%%%%%%% test %%%%%%%%%%%%
    SP.cellDensity = 1 / 1e-15; %%%%%%%%%%%% test %%%%%%%%%%%%
%     SP.cellDensity = 1 / 1e-14; %%%%%%%%%%%% test %%%%%%%%%%%%
    numCells = round(double(SP.vesselWidth * SP.vesselWidth * SP.vesselLength * SP.cellDensity)); % Round in case there are less than 1 cells...
    
    % MIGHT want randn so the distribution of cells is in the Gaussian profile?
    xr = rand([1, numCells])';
    yr = rand([1, numCells])';
    zr = rand([1, numCells])';
    
    
%     for n = 1:numCells
    %     Media.MP(n, :) = [ -vesselWidth / wl / 2 + xr(n) * vesselWidth / wl, 0, 50 + zr(n) * vesselHeight / wl, scatterReflectivity];
        pts = [ SP.xstart - SP.vesselWidth ./ 2 + xr .* SP.vesselWidth, ...
                SP.ystart - SP.vesselWidth ./ 2 + yr .* SP.vesselWidth, ...
                SP.zstart - SP.vesselLength ./ 2 + zr .* SP.vesselLength, ... % SP.zstart + SP.startDepthMM/1e3 + zr .* SP.vesselLength, ...
                repmat(SP.scatterReflectivity, length(xr), 1)];

%         pts = [ -vesselX / wl / 2 + xr(n) * vesselX / wl, 0, 0 + zr(n) * vesselZ / wl, scatterReflectivity];
        
        % figure; scatter3(pts(:, 1), pts(:, 2), pts(:, 3), '.'); axis square

%     end
    mask = sqrt(pts(:, 1).^2 + pts(:, 2).^2) > SP.vesselDiam/2; % Remove points outside of the vessel radius
    pts(mask, :) = [];
return