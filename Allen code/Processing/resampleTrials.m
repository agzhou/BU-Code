% Description:
%       Resample fUS data within a trial (and interpolate with splines)
% Inputs:
%       data: CBVi/PDI/etc. across all (super)frames in the experiment (2D
%             space + time or 3D space + time)
%       trial_sf: a (# trials x 1) cell array. Each cell contains the
%                 superframe indices within that trial.
%       trial_windows: a (# trials x 1) cell array. Each cell contains the
%                      DAQ sample indices within that trial.
%       sfStarts: a vector containing the DAQ sample indices corresponding
%                 to the start of each superframe acquisition
%       P: a structure with acquisition parameters (see
%          makeParameterStructure.m)
%       interp_factor: a scalar (natural number) >= 1 to divide the DAQ
%                      sampling rate by when interpolating the data
% Outputs:
%       data_resampled: data that is aligned for each trial and
%                       resampled + interpolated

function [data_resampled] = resampleTrials(data, trial_sf, trial_windows, sfStarts, P, interp_factor, varargin)

    % Add an optional input to choose the interpolation method (update
    % 8/25/25)
    interp_type = 'makima'; % Default interpolation method
    if nargin > 6
        interp_type = varargin{1};
    end

    if interp_factor < 1
        error('Interpolation factor needs to be at least 1')
    end
    if interp_factor ~= floor(interp_factor)
        interp_factor = round(interp_factor);
        warning('Interpolation factor is not a natural number, rounding')
    end

    % Resample and interpolate
    data_resampled = cell(size(trial_sf)); % Store each resampled trial individually
%     interp_factor = 100; % Factor by which to "decimate" the daq rate 
%     for interpolation timepoints
    
    interp_times = 1:interp_factor:P.daqrate * P.Mcr_fcp.apis.seq_length_s; % Time points at which we calculate an interpolated value

    data_numDims = length(size(data)); % # of dimensions of the input data
    
    switch data_numDims
        case 3                  % 3D: 2D space + time
            for trial = 1:length(trial_windows)
                disp("Resampling trial " + num2str(trial))
                temp_indices = sfStarts(trial_sf{trial});
                temp_indices_shifted = temp_indices - trial_windows{trial}(1) + 1; % Shift the indices so they correspond to a trial start at 1
                % data_resampled{trial} = spline(temp_indices_shifted, data(:, :, trial_sf{trial}), interp_times);
                try
                    switch interp_type
                        case 'spline'
                            data_resampled{trial} = spline(temp_indices_shifted, data(:, :, trial_sf{trial}), interp_times);
                        case 'makima'
                            data_resampled{trial} = makima(temp_indices_shifted, data(:, :, trial_sf{trial}), interp_times);
                        case 'pchip'
                            data_resampled{trial} = pchip(temp_indices_shifted, data(:, :, trial_sf{trial}), interp_times);
                    end
%                 data_resampled{trial} = makima(temp_indices_shifted, data(:, :, trial_sf{trial}), interp_times);
%                 data_resampled{trial} = pchip(temp_indices_shifted, data(:, :, trial_sf{trial}), interp_times);
                catch % Case where the masked points being all NaN causes an error with the chckxy internal helper function
                    temp_data = data(:, :, trial_sf{trial});
                    temp_data(isnan(temp_data)) = 0;
                    switch interp_type
                        case 'spline'
                            data_resampled{trial} = spline(temp_indices_shifted, temp_data, interp_times);
                        case 'makima'
                            data_resampled{trial} = makima(temp_indices_shifted, temp_data, interp_times);
                        case 'pchip'
                            data_resampled{trial} = pchip(temp_indices_shifted, temp_data, interp_times);
                    end
                    warning('Make sure the [masked] voxels are set to a value compatible with chckxy, e.g., not NaN')
                end
            end
        case 4                  % 4D: 3D space + time
            for trial = 1:length(trial_windows)
                disp("Resampling trial " + num2str(trial))
                temp_indices = sfStarts(trial_sf{trial});
                temp_indices_shifted = temp_indices - trial_windows{trial}(1) + 1; % Shift the indices so they correspond to a trial start at 1
                try
                    data_resampled{trial} = spline(temp_indices_shifted, data(:, :, :, trial_sf{trial}), interp_times);
%                 data_resampled{trial} = makima(temp_indices_shifted, data(:, :, :, trial_sf{trial}), interp_times);
%                 data_resampled{trial} = pchip(temp_indices_shifted, data(:, :, :, trial_sf{trial}), interp_times);
                catch % Case where the masked points being all NaN causes an error with the chckxy internal helper function
                    temp_data = data(:, :, :, trial_sf{trial});
                    temp_data(isnan(temp_data)) = 0;
                    data_resampled{trial} = spline(temp_indices_shifted, temp_data, interp_times);
                    warning('Make sure the [masked] voxels are set to a value compatible with chckxy, e.g., not NaN')
                end
            end
    end
end