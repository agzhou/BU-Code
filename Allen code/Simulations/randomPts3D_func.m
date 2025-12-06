function pts = randomPts3D_func(vesselX, vesselY, vesselZ, wl, startDepth, xstart, ystart, zstart)
    % x, y, start are measured with x, y = 0 at center of transducer
    % z start is distance starting from startDepth
    rng('shuffle')
    persistent cellDensity;
    persistent numCells;

%     % horizontal
%     vesselX = Trans.numelements*Trans.spacing * wl;
%     vesselZ = 30e-6; % say 30 um for now
    
%     % vertical
%     vesselX = 30e-6;    % x dimension
%     vesselZ = endDepthMM * 1e-3; % z dimension
    
    % cellDensity = 1 / 100 / 1e-18; % 1/(100 um^3) from Bingxue's paper
%     cellDensity = 1 / 1e-12; %%%%%%%%%%%% test %%%%%%%%%%%%
    cellDensity = 1 / 1e-15; %%%%%%%%%%%% test %%%%%%%%%%%%
    numCells = uint16(vesselX * vesselY * vesselZ * cellDensity)
    
    % MIGHT want randn so the distribution of cells is in the Gaussian profile?
    xr = rand([1, numCells])';
    yr = rand([1, numCells])';
    zr = rand([1, numCells])';
    scatterReflectivity = 1.0;
    
%     for n = 1:numCells
    %     Media.MP(n, :) = [ -vesselWidth / wl / 2 + xr(n) * vesselWidth / wl, 0, 50 + zr(n) * vesselHeight / wl, scatterReflectivity];
        pts = [ xstart - vesselX ./ wl ./ 2 + xr .* vesselX ./ wl, ...
                ystart - vesselY ./ wl ./ 2 + yr .* vesselY ./ wl, ...
                zstart + startDepth + zr .* vesselZ ./ wl, ...
                repmat(scatterReflectivity, length(xr), 1)];

%         pts = [ -vesselX / wl / 2 + xr(n) * vesselX / wl, 0, 0 + zr(n) * vesselZ / wl, scatterReflectivity];
    
%     end
return