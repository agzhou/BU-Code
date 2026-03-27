% Description: convert a list of angles in 1 dimension to a full list of
% angles across x and y. This assumes that for TX, the angles are all the C-R
% ([theta_x, 0]) and then all the R-C ([0, theta_y]), and the reverse for
% RX.

% Inputs:
%   - anglesList: list of angles in 1 dimension [# angles, 1]
%   - config: 'TX' or 'RX' options for transmit and receive, respectfully

function angles = listToAnglesRCA(anglesList, config)
    anglesList = anglesList(:); % Make into a column vector
    switch config
        case 'TX'
            angles = [[anglesList; zeros(size(anglesList))], [zeros(size(anglesList)); anglesList]]; % List of all transmit angles. Dimensions [ntaTX, 2 (x and y angle)]
        case 'RX'
            angles = [[zeros(size(anglesList)); anglesList], [anglesList; zeros(size(anglesList))]]; % List of all receive angles. Dimensions [ntaRX, 2 (x and y angle)]
    end
end