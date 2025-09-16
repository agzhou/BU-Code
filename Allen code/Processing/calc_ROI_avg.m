function roi_data_TA = calc_ROI_avg(data_TA, roi_region_inds)
    %     roi_data_TA = zeros(size(am_data));
    for ti = 1:size(data_TA, 4) % for 3D + time data
%     for ti = 1
         temp_data_TA = data_TA(:, :, :, ti);
         temp_roi_data_avg = mean(temp_data_TA(roi_region_inds));
         roi_data_TA(ti) = temp_roi_data_avg;
    end

end