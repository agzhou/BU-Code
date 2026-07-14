% Trial average [the relative change in] a hemodynamic parameter (assumed
% to be a cell array with each cell a separate trial with the same # of sample points)

% Output:
%       data_trial_average: a matrix of the trial averaged data
function [data_trial_average] = fUS_trialAverage(data)
    data_trial_average = data{1};
    if length(data) > 1 % If there are multiple trials to average
        for trial = 2:length(data)
            data_trial_average = data_trial_average + data{trial};
        end
    end
    data_trial_average = data_trial_average ./ length(data);
end