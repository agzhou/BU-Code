function pts = randomPts2D_func(vesselX, vesselZ, wl, startDepth, xstart, zstart)
    % zstart is # of wavelengths below startDepth to start top surface of the vessel
% function pts = randomPts2D_func(vesselX, vesselZ, wl, startDepth)
    rng('shuffle')
    persistent cellDensity;
    persistent numCells;

%     % horizontal
%     vesselX = Trans.numelements*Trans.spacing * wl;
%     vesselZ = 30e-6; % say 30 um for now
    
%     % vertical
%     vesselX = 30e-6;    % x dimension
%     vesselZ = endDepthMM * 1e-3; % z dimension
    
    cellDensity = 1 / 20 / 1e-12; % 1/(20 um^2) from Bingxue's paper
    numCells = uint16(vesselX * vesselZ * cellDensity);
    
    % MIGHT want randn so the distribution of cells is in the Gaussian profile?
    xr = rand([1, numCells])';
    zr = rand([1, numCells])';
    scatterReflectivity = 1.0;
    
%     for n = 1:numCells
    %     Media.MP(n, :) = [ -vesselWidth / wl / 2 + xr(n) * vesselWidth / wl, 0, 50 + zr(n) * vesselHeight / wl, scatterReflectivity];
        pts = [ xstart - vesselX ./ wl ./ 2 + xr .* vesselX ./ wl, ...
                zeros(length(xr), 1), ...
                zstart + startDepth + zr .* vesselZ ./ wl, ...
                repmat(scatterReflectivity, length(xr), 1)];

%         pts = [ -vesselX / wl / 2 + xr(n) * vesselX / wl, 0, 0 + zr(n) * vesselZ / wl, scatterReflectivity];
    
%     end
return