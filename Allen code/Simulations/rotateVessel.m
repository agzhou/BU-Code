% Vessel is scatterer data. Columns are x, y, z, reflection coefficient.
% Rows are each scatterer.

% Rotate the vessel about its center.

function [vessel_rotated] = rotateVessel(vessel, x_angle, y_angle, z_angle, SP)
    % Zero mean the scatterers before rotating
    vessel_zeromeaned = vessel;
    vessel_zeromeaned(:, 1) = vessel_zeromeaned(:, 1) - SP.xstart;
    vessel_zeromeaned(:, 2) = vessel_zeromeaned(:, 2) - SP.ystart;
    vessel_zeromeaned(:, 3) = vessel_zeromeaned(:, 3) - SP.zstart;

    % Perform the roation on the centered data
    R = rotationMatrix(x_angle, y_angle, z_angle);
    vessel_rotated = vessel_zeromeaned;
    vessel_rotated(:, 1:3) = (R * vessel_rotated(:, 1:3)')';

    % Re-translate the center in space
    vessel_rotated(:, 1) = vessel_rotated(:, 1) + SP.xstart;
    vessel_rotated(:, 2) = vessel_rotated(:, 2) + SP.ystart;
    vessel_rotated(:, 3) = vessel_rotated(:, 3) + SP.zstart;
        
end