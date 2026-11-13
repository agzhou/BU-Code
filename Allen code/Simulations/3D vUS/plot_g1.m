% Plotting helper function for vUS numerical simulations. 

% Inputs:
%   sIQ: vector of sIQ values for each time step (complex)
%   g1: vector of g1 values for each time step (complex)
%   tau: time lag vector
%   (Optional) tau_units: character vector of the units for the time lags
function plot_g1(sIQ, g1, tau, tau_units)
    if nargin > 3
        tau_units = [' [', tau_units, ']']; % If units are specified, add a space and brackets
    else
        tau_units = '';
    end

    % figure; plot(tau, real(sIQ - mean(sIQ))); hold on
    % plot(tau, imag(sIQ - mean(sIQ))); hold off
    % legend('real[sIQ - mean(sIQ)]', 'imag[sIQ - mean(sIQ)]')
    % xlabel(['tau', tau_units])


    figure; plot(tau, abs(g1)); xlabel(['tau', tau_units]); ylabel("|g_1|"); title('|g1|')
    figure; plot(real(g1), imag(g1), '-o'); xlim([-1, 1]); ylim([-1, 1]); axis square; title('Complex g1')

end