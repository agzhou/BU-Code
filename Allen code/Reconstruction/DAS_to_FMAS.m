%% Apply FMAS to delay and summed IQ data (pre-compounding)
% for row column array
% Output: FMAS IQ volumes (x, y, z, frames)
% Input: IQ volumes (x, y, z, 2*# angles, frames)

function [IQ_FMAS] = DAS_to_FMAS(IQ)
    %% Initialize variabes and separate R-C and C-R volumes
    [xp, yp, zp, nacq, nf] = size(IQ);
%     na = P.na;
    na = nacq/2;
    
    IQ_CR = IQ(:, :, :, 1:na);          % column row volumes
    IQ_RC = IQ(:, :, :, na + 1:2*na);   % row column volumes
    
    IQ_FMAS = zeros(xp, yp, zp, nf);        % initialize final FMAS image
    %% Do the FMAS
    for cri = 1:na       % column row index
        for rci = 1:na   % row column index
            temp = IQ_CR(:, :, :, cri, :) .* IQ_RC(:, :, :, rci, :);
            mag = sqrt(abs(temp));  % magnitude but sqrt to maintain units
            s = sign(temp);         % phase
            IQ_FMAS = IQ_FMAS + s .* mag;   % update the summed variabe
    
        end
    end
end