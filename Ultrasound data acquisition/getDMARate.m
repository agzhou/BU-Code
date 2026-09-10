% Helper function to get the DMA rate [GB/s] of some Verasonics connector
% plate types.

function [DMARate] = getDMARate(connectorPlate)
    
    switch connectorPlate
        case 'UTA-260D'
            % DMARate = 3.3; % [GB/s]
            DMARate = 6.6; % TESTING
        case 'UTA-260S'
            DMARate = 6.6; % [GB/s]
        case 'UTA-408GE'
            DMARate = 6.6; % [GB/s]
    end
end