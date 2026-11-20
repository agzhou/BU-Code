% Plotting helper function for vUS numerical simulations. 

% Inputs:
%   sIQ: vector of sIQ values for each time step (complex)
%   g1: vector of g1 values for each time step (complex)
%   tau: time lag vector
%   (Optional) tau_units: character vector of the units for the time lags
%   (Optional) v_vector: numeric vector of the velocities used
%   (Optional) v_units: character vector of the units for the velocities
function plot_g1(sIQ, g1, tau, tau_units, v_vector, v_units)
    lw = 2; % Line width
    fs = 20; % Font size (points)

    if nargin > 3
        tau_units = [' [', tau_units, ']']; % If units are specified, add a space and brackets
    else
        tau_units = '';
    end

    if nargin > 4 % if a velocity vector is specified (e.g., for a sweep)
        v_labels = {};
        if nargin > 5
            v_units = " " + v_units;
        else
            v_units = "";
        end
        for vi = 1:length(v_vector)
            v_labels{vi} = "v_{gp} = " + num2str(v_vector(vi)) + v_units;
        end
    end

    % figure; plot(tau, real(sIQ - mean(sIQ))); hold on
    % plot(tau, imag(sIQ - mean(sIQ))); hold off
    % legend('real[sIQ - mean(sIQ)]', 'imag[sIQ - mean(sIQ)]')
    % xlabel(['tau', tau_units])


    figure; plot(tau, abs(g1), 'LineWidth', lw); xlabel(['tau', tau_units]); ylabel("|g_1|"); title('|g1|'); fontsize(fs, 'points'); legend(v_labels)
    figure; plot(real(g1), imag(g1), '-o', 'LineWidth', lw); xlim([-1, 1]); ylim([-1, 1]); axis square; title('Complex g1'); fontsize(fs, 'points'); legend(v_labels)

end