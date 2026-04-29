% Description: Store the ROI-averaged timecourses in matrix form (from a
%              cell array)
function [data_ROI_timecourses_mat, data_ROI_hemis_timecourses_mat] = ROIAvgTimecourses2mat(data_ROI_timecourses, data_ROI_hemis_timecourses, roi)
    data_ROI_timecourses_mat = zeros(length(data_ROI_timecourses{1}), roi.num_regions); % Still ROI-averaged PDI timecourses, but in matrix form (each column is a separate ROI timecourse). Dimensions: [# time points, # ROIs]
    for ri = 1:roi.num_regions % region/ROI index -- loop through each region
        data_ROI_timecourses_mat(:, ri) = data_ROI_timecourses{ri};
    end
    
    % Hemisphere-separated version
    data_ROI_hemis_timecourses_mat = zeros(length(data_ROI_timecourses{1}), roi.num_regions*2); % Still ROI-averaged PDI timecourses, but in matrix form (each column is a separate ROI timecourse). Dimensions: [# time points, # ROIs]
    for ri = 1:roi.num_regions % region/ROI index -- loop through each region
        data_ROI_hemis_timecourses_mat(:, (ri - 1)*2 + 1) = data_ROI_hemis_timecourses{ri, 1};
        data_ROI_hemis_timecourses_mat(:, ri*2) = data_ROI_hemis_timecourses{ri, 2};
    end