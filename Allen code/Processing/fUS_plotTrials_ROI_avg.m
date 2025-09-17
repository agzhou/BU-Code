% Plot each trial of a hemodynamic parameter.
% 
% Inputs:
%       data: a cell array with each cell a separate trial with the same # of sample points
%       point: a vector representing a point to plot the timecourses at

% Output:

function [] = fUS_plotTrials_ROI_avg(data, roi_indices)
    figure; legend;
    hold on
    for trial = 1:length(data)
        roi_data = calc_ROI_avg(data{trial}, roi_indices);
        plot(squeeze(roi_data))
    end
    hold off
end