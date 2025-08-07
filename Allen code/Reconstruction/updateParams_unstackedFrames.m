function [P_unstacked] = updateParams_unstackedFrames(P)
    
    % Derive some parameters
    try rcvChunkSize = P.rcvChunkSize;
    catch
    end

    if exist('rcvChunkSize', 'var')
        nsf = rcvChunkSize; % number of stacked frames
    else
        nsf = P.numFramesPerBuffer;
    end

    % Change 8/7/25
    if mod(nsf, 2) ~= 0
        warning('# of frames is odd, cutting the last frame so VSX can handle it')
        nsf = nsf - 1;
    end

    nspa = P.Receive(1).endSample; % # samples per acquisition

    % Set up flags according to which type of probe is used
    if isequal(P.Trans.name, 'RC15gV')
        probeType = 'RCA';
    elseif isequal(P.Trans.name, 'L22-14v')
        probeType = 'LA';
    else
        error('Probe name not coded into the script')
    end
    
    % Update param structures
    switch probeType
        case 'RCA'
            nspf = nspa * P.na * 2; % # samples per frame
            
            % Redefine the parameters structure with the unstacked parameters
            if exist('rcvChunkSize', 'var')
                P_unstacked = P;
                P_unstacked.numFramesPerBuffer = P.numFramesPerBuffer;
                P_unstacked.Resource.RcvBuffer.numFrames = P_unstacked.numFramesPerBuffer;
                P_unstacked.Resource.RcvBuffer.rowsPerFrame = nspf;
                P_unstacked.Receive = updateReceiveStructure_RCA(P_unstacked);
            else
                P_unstacked = P;
                P_unstacked.numFramesPerBuffer = nsf;
                P_unstacked.Resource.RcvBuffer.numFrames = nsf;
                P_unstacked.Resource.RcvBuffer.rowsPerFrame = nspf;
                P_unstacked.Receive = updateReceiveStructure_RCA(P_unstacked);
            end

        case 'LA'
            nspf = nspa * P.na; % # samples per frame
    
            % Redefine the parameters structure with the unstacked parameters
            if exist('rcvChunkSize', 'var')
                P_unstacked = P;
                P_unstacked.numFramesPerBuffer = P.numFramesPerBuffer;
                P_unstacked.Resource.RcvBuffer.numFrames = P_unstacked.numFramesPerBuffer;
                P_unstacked.Resource.RcvBuffer.rowsPerFrame = nspf;
                P_unstacked.Receive = updateReceiveStructure_LA(P_unstacked);
            else
                P_unstacked = P;
                P_unstacked.numFramesPerBuffer = nsf;
                P_unstacked.Resource.RcvBuffer.numFrames = nsf;
                P_unstacked.Resource.RcvBuffer.rowsPerFrame = nspf;
                P_unstacked.Receive = updateReceiveStructure_LA(P_unstacked);
            end
    end
    
end


%% Helper function for setting the Receive structure
% P_unstacked.Receive(1).
function [Receive] = updateReceiveStructure_RCA(P_unstacked)
    pair = 2;
    Receive = repmat(struct('Apod', zeros(1, P_unstacked.Trans.numelements), ... 
                            'startDepth', P_unstacked.Receive(1).startDepth, ...
                            'endDepth', P_unstacked.Receive(1).endDepth, ...
                            'TGC', P_unstacked.Receive(1).TGC, ...
                            'bufnum', P_unstacked.Receive(1).bufnum, ...
                            'framenum', P_unstacked.Receive(1).framenum, ...
                            'acqNum', P_unstacked.Receive(1).acqNum, ...
                            'sampleMode', P_unstacked.Receive(1).sampleMode, ...
                            'mode', P_unstacked.Receive(1).mode, ...
                            'callMediaFunc', P_unstacked.Receive(1).callMediaFunc, ...
                            'LowPassCoef', P_unstacked.Receive(1).LowPassCoef, ...
                            'InputFilter', P_unstacked.Receive(1).InputFilter), 1, pair * P_unstacked.na * P_unstacked.numFramesPerBuffer);
    j = 1;
    % an = 0;
    for nf = 1:P_unstacked.numFramesPerBuffer
        % Move points after all the acquisitions for one frame
%         Receive(j).callMediaFunc = movePointsOrNot;
        an = 0;

        for n = 1:P_unstacked.na
            an = an + 1;
            Receive(j).framenum = nf;
            Receive(j).acqNum = an;
            Receive(j).Apod(P_unstacked.Trans.numelements/2 + 1 : end) = ones(1, P_unstacked.Trans.numelements/2);
            j = j + 1;
        end
    
        for n = 1:P_unstacked.na
            an = an + 1;
            Receive(j).framenum = nf;
            Receive(j).acqNum = an;
            Receive(j).Apod(1:P_unstacked.Trans.numelements/2) = ones(1, P_unstacked.Trans.numelements/2);
            j = j + 1;
        end
        
    end
end

function [Receive] = updateReceiveStructure_LA(P_unstacked)
    Receive = repmat(struct('Apod', ones(1, P_unstacked.Trans.numelements), ... 
                            'startDepth', P_unstacked.Receive(1).startDepth, ...
                            'endDepth', P_unstacked.Receive(1).endDepth, ...
                            'TGC', P_unstacked.Receive(1).TGC, ...
                            'bufnum', P_unstacked.Receive(1).bufnum, ...
                            'framenum', P_unstacked.Receive(1).framenum, ...
                            'acqNum', P_unstacked.Receive(1).acqNum, ...
                            'sampleMode', P_unstacked.Receive(1).sampleMode, ...
                            'mode', P_unstacked.Receive(1).mode, ...
                            'callMediaFunc', P_unstacked.Receive(1).callMediaFunc, ...
                            'LowPassCoef', P_unstacked.Receive(1).LowPassCoef, ...
                            'InputFilter', P_unstacked.Receive(1).InputFilter), 1, P_unstacked.na * P_unstacked.numFramesPerBuffer);
    j = 1;
    % an = 0;
    for nf = 1:P_unstacked.numFramesPerBuffer
        % Move points after all the acquisitions for one frame
%         Receive(j).callMediaFunc = movePointsOrNot;
        an = 0;

        for n = 1:P_unstacked.na
            an = an + 1;
            Receive(j).framenum = nf;
            Receive(j).acqNum = an;
            j = j + 1;
        end
        
    end
end