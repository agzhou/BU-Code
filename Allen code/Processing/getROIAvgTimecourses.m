% Description: recover ROI-averaged hemodynamic timecourses

% Inputs:
%   - data: 4D (3D space + time) hemodynamic data that is registered to the
%           atlas that 'roi' is based on
%   - roi: a struct that has all the atlas/ROI information

function [data_ROI_timecourses, data_ROI_hemis_timecourses] = getROIAvgTimecourses(data, roi)
    % (old incorrect title: Correlation without resampling PDI in time (with sliding window))
    % figure; plot(sfTimeTags)
    % figure; plot(diff(sfTimeTags))
    num_tps = size(data, 4); % # of time points
    
    % PDI_ROI_timecourses = cell(num_regions, num_sf); % Store average ROI PDI timecourses in a cell array (each cell is an average timecourse)
    data_ROI_timecourses = cell(roi.num_regions, 1); % Store average ROI PDI timecourses in a cell array (each cell is an average timecourse)
    data_ROI_hemis_timecourses = cell(roi.num_regions, 2); % Store average ROI PDI (hemisphere-separated) timecourses in a cell array (each cell is an average timecourse)
    
    % -- Calculate PDI [ROI average] timecourses -- %
    tic
    for ti = 1:num_tps % "time" index -- go through each superframe
    % for ti = 11
        disp(ti)
        data_ti_temp = squeeze(data(:, :, :, ti)); % Registered volume at "time" index ti
    
        for ri = 1:roi.num_regions % region/ROI index -- loop through each region
    
            % Normal
            ROI_mask_temp = roi.masks_50um{ri}; % ROI #ri mask
            data_ri_masked_temp = data_ti_temp(ROI_mask_temp); % Vectorized voxels of the registered PDI at "time" index ti
            data_ROI_timecourses{ri}(ti) = mean(data_ri_masked_temp);
            
            % Hemisphere-separated
            ROI_mask_temp_left = roi.masks_50um_hemis{ri, 1}; % ROI #ri mask (left)
            ROI_mask_temp_right = roi.masks_50um_hemis{ri, 2}; % ROI #ri mask (right)
            data_ri_masked_temp_left = data_ti_temp(ROI_mask_temp_left); % Vectorized voxels of the registered PDI at "time" index ti
            data_ri_masked_temp_right = data_ti_temp(ROI_mask_temp_right); % Vectorized voxels of the registered PDI at "time" index ti
            data_ROI_hemis_timecourses{ri, 1}(ti) = mean(data_ri_masked_temp_left);
            data_ROI_hemis_timecourses{ri, 2}(ti) = mean(data_ri_masked_temp_right);
    
        end
    end
    % clearvars ti ri PDIallSF_reg_ti_temp ROI_mask_temp PDI_ri_masked_temp
    toc
    
    % Make any row timecourses into column vectors
    for ri = 1:roi.num_regions
        data_ROI_timecourses{ri} = squeeze(data_ROI_timecourses{ri}');
    
        data_ROI_hemis_timecourses{ri, 1} = squeeze(data_ROI_hemis_timecourses{ri, 1}');
        data_ROI_hemis_timecourses{ri, 2} = squeeze(data_ROI_hemis_timecourses{ri, 2}');
    
    end
    % clearvars ri ROI_mask_temp PDI_ri_masked_temp ROI_mask_temp_left ROI_mask_temp_right PDI_ri_masked_temp_left PDI_ri_masked_temp_right
end