% Description: 
%   Create a mask based on g1-based SNR for pixels/voxels to keep
%   (This is for 2D or 3D vUS)

% Input:
%   g1: complex g1 data (2D: [nz, nx, nTau]; 3D: [nx, ny, nz, nTau]). nTau must be at least 3. Note: my
%       code outputs g1 starting with tau = 0, so tau1 corresponds to index 2.
%   PP: Processing Parameters struct with at least the fields:
%       PP.fDim: dimension corresponding to frequency (or time)
%       PP.xDim: dimension corresponding to x in space
%       PP.yDim: dimension corresponding to y in space (3D only)
%       PP.zDim: dimension corresponding to z in space
%       PP.dimensionality: 2 for 2D, 3 for 3D
%   type: 'full' or 'half' --> dictates some constants used, depending on
%         if we are using the whole frequency spectrum or just part of it (the
%         directional filtering)

% Output:
%   g1SNR: g1-based SNR values (2D: [nz, nx, nTau]; 3D: [nx, ny, nz, nTau])
%   g1SNR_mask: g1-based SNR mask (2D: [nz, nx]; 3D: [nx, ny, nz])
%   g1SNR_express_mask: g1-based SNR mask, based on |g1(tau1)| (2D: [nz, nx]; 3D: [nx, ny, nz])

function [g1SNR, g1SNR_mask, g1SNR_express_mask] = g1BasedSNR(g1, PP, type)

    % Adjust some constants based on the 'type' input
    switch type
        case 'full'
            const = 0.4;
            floor = 0.08;
            express_threshold = 0.4;
        case 'half'
            const = 0.4;
            floor = 0.25;
            express_threshold = 0.6;
    end
    
    switch PP.dimensionality
        case 2 % 2D
            % g1SNR = mean(real(g1(:, :, 2:3)), fDim); % Average real(g1) over the first two time lags --> proxy for the SNR, based on the g1
            g1SNR = mean(abs(g1(:, :, 2:3)), PP.fDim); % Average abs(g1) over the first two time lags --> proxy for the SNR, based on the g1
            g1SNR_pixel_avg = mean(g1SNR, [PP.zDim, PP.xDim]); % Pixel/Voxel dimensions should be 1-2
            g1SNR_std = std(g1SNR, 0, [PP.zDim, PP.xDim]); % Pixel/Voxel dimensions should be 1-2
            g1SNR_mask = g1SNR > max( (g1SNR_pixel_avg - const*g1SNR_std), floor );
            g1SNR_express_mask = abs(g1(:, :, 2)) > express_threshold; % Immediate pass mask according to |g1(tau1)| > threshold

        case 3 % 3D
            % g1SNR = mean(real(g1(:, :, :, 2:3)), fDim); % Average real(g1) over the first two time lags --> proxy for the SNR, based on the g1
            g1SNR = mean(abs(g1(:, :, :, 2:3)), PP.fDim); % Average abs(g1) over the first two time lags --> proxy for the SNR, based on the g1
            g1SNR_pixel_avg = mean(g1SNR, [PP.xDim, PP.yDim, PP.zDim]); % Pixel/Voxel dimensions should be 1-3
            g1SNR_std = std(g1SNR, 0, [PP.xDim, PP.yDim, PP.zDim]); % Pixel/Voxel dimensions should be 1-3
            g1SNR_mask = g1SNR > max( (g1SNR_pixel_avg - const*g1SNR_std), floor );
            g1SNR_express_mask = abs(g1(:, :, :, 2)) > express_threshold; % Immediate pass mask according to |g1(tau1)| > threshold

    end
    
end