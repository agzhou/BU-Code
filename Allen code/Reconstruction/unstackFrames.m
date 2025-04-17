function [RcvData_unstacked, P_unstacked] = unstackFrames(RcvData, P)
    %% Use parallel processing for speed
    % https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox
    
%     pp = gcp('nocreate');
%     if isempty(pp)
%         % There is no parallel pool
%         parpool LocalProfile1
%     
%     end

    %%
    try rcvChunkSize = P.rcvChunkSize;
    catch
    end

    if exist('rcvChunkSize', 'var')
%         nsf = rcvChunkSize; % number of stacked frames
        nsf = P.numFramesPerBuffer;
    else
        nsf = P.numFramesPerBuffer;
    end

    nspa = P.Receive(1).endSample; % # samples per acquisition
    nspf = nspa * P.na * 2; % # samples per frame
    
    RcvData_unstacked = zeros(nspf, P.Resource.Parameters.numRcvChannels, nsf, 'int16');

    % Unstack the frames
%     RcvData_unstacked = reshape(RcvData, [nspf, P.Resource.Parameters.numRcvChannels, nsf]);
    
    if length(size(RcvData)) == 2 % all frames are stacked into one "frame" in the buffer
        for n = 1:nsf
            RcvData_unstacked(:, :, n) = RcvData((n-1) * nspf + 1 : n * nspf, :);
        end
    else % if only some frames are stacked, RcvData should have 3 dimensions
        for n = 1:rcvChunkSize:nsf
            for c = 0:rcvChunkSize - 1
%                 RcvData_unstacked(:, :, n + c) = RcvData((n-1) * nspf + 1 : n * nspf, :, floor(n/rcvChunkSize) + 1);
                RcvData_unstacked(:, :, n + c) = RcvData(c * nspf + 1 : (c + 1) * nspf, :, floor(n/rcvChunkSize) + 1);
            end
    
        end
    end
    
    % Remap according to the element - channel map (comment the line below
    % if using Verasonics recon)
%     RcvData_unstacked = RcvData_unstacked(:, P.Trans.Connector, :);

end

% figure; imagesc(abs(squeeze(double(RcvData_unstacked(:, :, 1))) .^ 0.5))
