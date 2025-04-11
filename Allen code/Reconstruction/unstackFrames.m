% function [RcvData_unstacked] = unstackFrames(RcvData, P)
    %% Use parallel processing for speed
    % https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox
    
    pp = gcp('nocreate');
    if isempty(pp)
        % There is no parallel pool
        parpool LocalProfile1
    
    end

    %%
    if exist('P.rcvChunkSize', 'var')
        nsf = P.rcvChunkSize; % number of stacked frames
    else
        nsf = P.numFramesPerBuffer;
    end

    nspa = P.Receive(1).endSample; % # samples per acquisition
    nspf = nspa * P.na * 2; % # samples per frame
    
    RcvData_unstacked = zeros(nspf, P.Resource.Parameters.numRcvChannels, nsf);

    % ... unstack
%     RcvData_unstacked = reshape(RcvData, [nspf, P.Resource.Parameters.numRcvChannels, nsf]);
    
    for n = 1:nsf
        RcvData_unstacked(:, :, n) = RcvData((n-1) * nspf + 1 : n * nspf, :);


    end
    
    % Remap according to the element - channel map
    RcvData_unstacked = RcvData_unstacked(:, P.Trans.Connector, :);


% end

% figure; imagesc(abs(squeeze(double(RcvData_unstacked(:, :, 1))) .^ 0.5))