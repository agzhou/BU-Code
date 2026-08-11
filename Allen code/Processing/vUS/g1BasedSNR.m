% Description: 
%   Create a mask based on g1-based SNR for pixels to keep
%   (This is for 2D vUS)
% Input:
%   g1: complex g1 data [nz, nx, nTau]. nTau must be at least 3. Note: my
%       code outputs g1 starting with tau = 0, so tau1 corresponds to index 2.
%   PP: Processing Parameters struct with at least the fields:
%       PP.fDim: dimension corresponding to frequency (or time)
%       PP.xDim: dimension corresponding to x in space
%       PP.zDim: dimension corresponding to z in space
%   type: 'full' or 'half' --> dictates some constants used, depending on
%         if we are using the whole frequency spectrum or just part of it (the
%         directional filtering)


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

    % g1SNR = mean(real(g1(:, :, 2:3)), fDim); % Average real(g1) over the first two time lags --> proxy for the SNR, based on the g1
    g1SNR = mean(abs(g1(:, :, 2:3)), PP.fDim); % Average abs(g1) over the first two time lags --> proxy for the SNR, based on the g1
    g1SNR_pixel_avg = mean(g1SNR, [PP.zDim, PP.xDim]);
    g1SNR_std = std(g1SNR, 0, [1, 2]);
    g1SNR_mask = g1SNR > max( (g1SNR_pixel_avg - const*g1SNR_std), floor );
    g1SNR_express_mask = abs(g1(:, :, 2)) > express_threshold; % Immediate pass mask according to |g1(tau1)| > threshold

end