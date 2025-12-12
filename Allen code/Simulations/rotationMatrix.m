% Generate a rotation matrix R in 3D
function R = rotationMatrix(x_angle, y_angle, z_angle)

    Rx = rotx(x_angle); % x rotation matrix (input is in degrees)
    Ry = roty(y_angle); % y rotation matrix (input is in degrees)
    Rz = rotz(z_angle); % z rotation matrix (input is in degrees)
    R = Rx * Ry * Rz;   % Combined rotation matrix
end