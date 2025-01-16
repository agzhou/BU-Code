% Take the data with each subframe stacked and turn that into another
% dimension



%% Updates
% 11/6/24: should work fine with a single superframe now
%%
function [shapedRcvData] = reshapeRcvData_ULM(RcvData, P)
%     numSupFrames = P.numSupFrames;
    numSubFrames = P.numSubFrames;
    s = size(RcvData);
    nzspa = s(1) ./ numSubFrames; % # z samples per acquisition

    if length(s) == 2 % if only 1 frame and there's no 3rd dimension
        shapedRcvData = zeros(nzspa, s(2), numSubFrames, 1);
    else
        shapedRcvData = zeros(nzspa, s(2), numSubFrames, s(3));

    end
    
%     for nsupf = 1:numSubFrames
        for nsubf = 1:numSubFrames
%             shapedRcvData(:, :, nsubf, nsupf) = RcvData((nsubf - 1) * nzspa + 1 : nsubf*nzspa, :, nsupf);
            shapedRcvData(:, :, nsubf, :) = RcvData((nsubf - 1) * nzspa + 1 : nsubf*nzspa, :, :);
        end
%     end
end
%% old 11/6/24
% function [shapedRcvData] = reshapeRcvData(RcvData, P)
% %     numSupFrames = P.numSupFrames;
%     numSubFrames = P.numSubFrames;
%     s = size(RcvData);
%     
% %     shapedRcvData = reshape(RcvData, s(1) ./ numSubFrames, s(2), numSubFrames, s(3));
%     nzspa = s(1) ./ numSubFrames; % # z samples per acquisition
%     shapedRcvData = zeros(nzspa, s(2), numSubFrames, s(3));
% %     for nsupf = 1:numSubFrames
%         for nsubf = 1:numSubFrames
% %             shapedRcvData(:, :, nsubf, nsupf) = RcvData((nsubf - 1) * nzspa + 1 : nsubf*nzspa, :, nsupf);
%             shapedRcvData(:, :, nsubf, :) = RcvData((nsubf - 1) * nzspa + 1 : nsubf*nzspa, :, :);
%         end
% %     end
% end
