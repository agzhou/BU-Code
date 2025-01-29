%% Create rolling re-sampled frames on delay and summed IQ data (pre-compounding)
% for linear array

% Input: IQ data (z pixels, x pixels, # angles, # frames)
% Output: rolled IQ data (z pixels, x pixels, # angles * # frames - # angles + 1)

function [IQr] = LA_rollingFrames(IQ)
    [zp, xp, nacq, nf] = size(IQ);

    IQr = zeros(zp, xp, nacq * nf - nacq + 1);        % initiaize rolling IQ variable
    
    for rfi = 1:nacq * nf - nacq + 1                % rolling frame index
        m = mod(rfi - 1, nacq) + 1;
        f = floor((rfi - 1)/nacq) + 1; % the "original frame" index
        if m == 1
            temp = IQ(:, :, m:nacq, f);
        else
            temp = cat(3, IQ(:, :, m:nacq, f), IQ(:, :, 1 : m - 1, f + 1));
        end
        IQr(:, :, rfi) = sum(temp, 3);
    
    end
end