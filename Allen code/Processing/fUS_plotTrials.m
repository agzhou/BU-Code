% Plot each trial of a hemodynamic parameter.
% 
% Inputs:
%       data: a cell array with each cell a separate trial with the same # of sample points
%       point: a vector representing a point to plot the timecourses at

% Output:

function [] = fUS_plotTrials(data, point)
    figure; legend;
    hold on
    for trial = 1:length(data)
        plot(squeeze(data{trial}(point(1), point(2), point(3), :)))
    end
    hold off
end