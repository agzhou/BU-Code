%% Description: plot direction vector of transmit or receive angles
function plotAngles(angles)
    origin = [0, 0, 0];
    figure; hold on
    for ai = 1:size(angles, 1)
        theta = angles(ai, :); % [theta_x, theta_y, theta_z]
        % u = [sin(theta(2)), -sin(theta(1)), cos(theta(1))*cos(theta(2))];
        % u = [cos(theta(1))*sin(theta(2)), -sin(theta(1))*cos(theta(2)), cos(theta(1))*cos(theta(2))];
        u = [sin(theta(2)), -sin(theta(1))*cos(theta(2)), cos(theta(1))*cos(theta(2))];
        plot3([origin(1), origin(1) + u(1)], [origin(2), origin(2) + u(2)], [origin(3), origin(3) + u(3)] ,  'o-')
        % disp(norm(u))
        % disp(sqrt(u(1)^2 + u(2)^2 + u(3)^2))
    end
    axis equal

end