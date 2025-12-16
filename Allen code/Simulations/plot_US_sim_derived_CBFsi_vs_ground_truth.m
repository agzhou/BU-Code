


%% Plot the ultrasound-derived CBFsi (vessel avg) vs. v_set
figure
plot(v_set, CBFsi_va, 'o', 'MarkerSize', 11, 'LineWidth', 4)
hold on
%% Least squares linear fit
A = [v_set', ones(length(v_set), 1)]; % [y = slope * v_actual + intercept]
x_lsq = (A'*A) \ (A'*CBFsi_va');

%% Calculate R^2

SS_res = sum( (CBFsi_va - ( x_lsq(1) .* v_set + x_lsq(2) )).^2 );
SS_tot = sum( (CBFsi_va - mean(CBFsi_va)).^2 );

R2 = 1 - SS_res/SS_tot

%% Add least squares fit to the plot
v_grid = linspace(min(v_set), max(v_set), 100);
v_lsq = x_lsq(1) .* v_grid + x_lsq(2);
plot(v_grid, v_lsq, '-', 'LineWidth', 2)
hold off



xlabel("v_{set} [mm/s]")
ylabel("CBF_{speed} index")
title("Set flow velocity vs. ultrasound CBF_{speed} index")
legend("Ultrasound simulation", "Least-squares linear fit: CBF_{speed} index = " + num2str(x_lsq(1)) + " * v_{set} + " + num2str(x_lsq(2)) + " (R^2 = " + num2str(R2) + ")")
fontsize(20, "points")

xlim([0, max(v_set) + 10])
% ylim([0, max(v_actual) + 10])
axis square
grid on

