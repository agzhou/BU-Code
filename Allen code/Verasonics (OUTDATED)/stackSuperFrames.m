function [RcvDataStacked, P_new] = stackSuperFrames(RcvData, P, numSupFramesToUse)
    [s, nch, ~] = size(RcvData);
%     RcvDataStacked = zeros(s*P.numSupFrames, nch);
    RcvDataStacked = zeros(s*numSupFramesToUse, nch);

%     for supf = 1:P.numSupFrames
    for supf = 1:numSupFramesToUse
        RcvDataStacked((supf - 1) * s + 1 : supf * s, :) = RcvData(:, :, supf);

    end

    P_new = P;
%     P_new.numSubFrames = P.numSupFrames * P.numSubFrames;
    P_new.numSubFrames = P.numSubFrames * numSupFramesToUse;

    P_new.numSupFrames = 1;
end