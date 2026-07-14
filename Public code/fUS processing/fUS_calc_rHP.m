% Calculate r(Hemodynamic parameter) -- relative change

% Inputs:
%       data: a (# trials x 1) cell array of CBVi/PDI/etc. across time (2D
%             space + time or 3D space + time). Time is assumed to be the
%             last dimension.
% Outputs:
%       data_baseline: baseline average
%       data_relative_change: Relative change of the data relative to the
%                             baseline average

function [data_baseline, data_relative_change] = fUS_calc_rHP(data, P, interp_factor)
    
    data_numDims = length(size(data{1})); % # of dimensions of the input data

    data_baseline = cell(size(data)); % Baseline average for each trial
    data_relative_change = cell(size(data)); % Relative change for each trial
    
    switch data_numDims
        case 3                  % 3D: 2D space + time
            for trial = 1:length(data)
                data_baseline{trial} = mean(data{trial}(:, :, 1 : P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor), data_numDims);
                % data_relative_change{trial} = (data{trial} - data_baseline{trial}) ./ data_baseline{trial};
                data_relative_change{trial} = (data{trial}) ./ data_baseline{trial};
            end
        case 4                  % 4D: 3D space + time
            for trial = 1:length(data)
                data_baseline{trial} = mean(data{trial}(:, :, :, 1 : P.Mcr_fcp.apis.delay_time_ms/1000 * P.daqrate / interp_factor), data_numDims);
                % data_relative_change{trial} = (data{trial} - data_baseline{trial}) ./ data_baseline{trial};
                data_relative_change{trial} = (data{trial}) ./ data_baseline{trial};
            end
    end
end