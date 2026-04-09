% Description: plot the delta angles for RCA RX-TX

function plotDeltaAngles_RCA(anglesTX, anglesRX)
    naTX = size(anglesTX, 1)/2;
    naRX = size(anglesRX, 1)/2;
    delta_angles_deg = rad2deg(calcDeltaThetaRCA(anglesTX, anglesRX));

    L = size(delta_angles_deg, 1)/2;
    delta_angles_deg_CR = delta_angles_deg(1:L, :);
    delta_angles_deg_RC = delta_angles_deg(L + 1:2*L, :);

    figure; hold on
    % plot(delta_angles_deg(:, 1), delta_angles_deg(:, 2), 'o')
    plot(delta_angles_deg_CR(:, 1), delta_angles_deg_CR(:, 2), 'o', 'MarkerSize', 8, 'LineWidth', 2)
    plot(delta_angles_deg_RC(:, 1), delta_angles_deg_RC(:, 2), 'x', 'MarkerSize', 8, 'LineWidth', 2)

    axis image; title('Delta angles'); xlabel('x angle [deg]'); ylabel('y angle [deg]'); fontsize(20, 'points')
    hold off
end