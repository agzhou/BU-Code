% Description: calculates the per-frame displacement of some simulated
%              points, assuming the vessel is cylindrical along the axis
%              specified by axisDim
function [distPerFrame] = tangelderFlowProfile(pts, v_gp, frameRate, k, a, axisDim, vesselDiam, vesselCenter)
    
    % axisDim = SP.dim;
    R = vesselDiam; % Vessel diameter
    axisCenterCoord = vesselCenter(axisDim); % Coordinate of the dimension SP.dim corresponding to the vessel center
    otherDims = [1, 2, 3]; otherDims(axisDim) = []; % Get the non-axis dims
    % r = sqrt( (pts(:, axisDim).^2 + axisCenterCoord.^2) ); % Calculate each point's radial position
    r = sqrt( (pts(:, otherDims(1)) - vesselCenter(otherDims(1))).^2 + (pts(:, otherDims(2)) - vesselCenter(otherDims(2))).^2 ); % Calculate each point's radial position

    v_max = v_gp ./ (1 - 2/(2 + k) * a^k);

    distPerFrame = v_max.* (1 - (a.*(r./R)).^k) ./ frameRate;

end