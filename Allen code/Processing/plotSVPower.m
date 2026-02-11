% Description: Plot the power of singular values for every time point in a
% measurement.

% Input:
%   - SVsallPoints is a 2D matrix; each column is the singular value
%     magnitude across all frames of the "superframe" or block
%   - sv_threshold_lower: lower singular value threshold to keep
%   - sv_threshold_upper: upper singular value threshold to keep

function SVPowers = plotSVPower(SVsallPoints, sv_threshold_lower, sv_threshold_upper)
    SVPowers = squeeze(sum(SVsallPoints(sv_threshold_lower:sv_threshold_upper, :), 1));
    figure
    plot(SVPowers, 'LineWidth', 2)
    xlabel('Superframe/Block number')
    ylabel('Singular value power')

end