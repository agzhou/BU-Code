% Description: extract unit direction vector for a plane wave given the
% angles used

% Input:
%   - angles: [n x 2] list of x and y angles in radians

function [u] = getUnitVector(angles)
    % y rotation and then x rotation
    u = [sin(theta(2)), -sin(theta(1))*cos(theta(2)), cos(theta(1))*cos(theta(2))];

    % x rotation and then y rotation
    % u = [cos(theta(1)) * sin(theta(2)), -sin(theta(1)), cos(theta(1))*cos(theta(2))];

end