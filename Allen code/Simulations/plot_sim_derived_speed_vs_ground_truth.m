v_actual = [5, 10, 20, 30, 40, 50]'; % [mm/s]
v_sim_fit = [4.6566, 10.174, 20.7116, 30.0592, 42.377, 51.1006]';
figure
plot(v_actual, v_sim_fit, 'o', 'MarkerSize', 11, 'LineWidth', 4)
hold on
%% Least squares linear fit
A = [v_actual, ones(length(v_actual), 1)]; % [y = slope * v_actual + intercept]
x_lsq = (A'*A) \ (A'*v_sim_fit);

%% Calculate R^2

SS_res = sum( (v_sim_fit - ( x_lsq(1) .* v_actual + x_lsq(2) )).^2 );
SS_tot = sum( (v_sim_fit - mean(v_sim_fit)).^2 );

R2 = 1 - SS_res/SS_tot

%% Add least squares fit to the plot
v_grid = linspace(min(v_actual), max(v_actual), 100);
v_lsq = x_lsq(1) .* v_grid + x_lsq(2);
plot(v_grid, v_lsq, '-', 'LineWidth', 2)
hold off



xlabel("v_{set} [mm/s]")
ylabel("v_{fit} [mm/s]")
title("Set vs. fitted flow velocity")
legend("Simulation", "Least-squares linear fit: v_{fit} = " + num2str(x_lsq(1)) + " * v_{set} + " + num2str(x_lsq(2)) + " (R^2 = " + num2str(R2) + ")")
fontsize(20, "points")

xlim([0, max(v_actual) + 10])
ylim([0, max(v_actual) + 10])
axis square
grid on

