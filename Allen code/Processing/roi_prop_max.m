function [roi_region_inds] = roi_prop_max(am_data, fraction)

    [am_max, am_max_ind] = max(am_data, [], 'all'); % Find the maximum point of the activation map
    am_halfmax = am_max * fraction;

%     [i, j, k] = ind2sub(size(am_data), am_max_ind); % get the yxz coordinates of the maximum point
%     figure; plot(squeeze(data_TA(i, j, k, :)))

    roi = am_data;
%     roi(roi < am_halfmax) = 0; % Remove the values less than the half max
    roi_bw = roi > am_halfmax; % Create a binary volume where voxels of the activation map > halfmax are true, and otherwise false
    CC = bwconncomp(roi_bw); % Get the connected component information

    % Go through each connected region and see which one contains the max
    % point
    CC_region_ind = [];
    for ccri = 1:length(CC.PixelIdxList)% connected component region index
        if any(CC.PixelIdxList{ccri} == am_max_ind)
            CC_region_ind = ccri;
            break
        end
    end

    roi_region_inds = CC.PixelIdxList{CC_region_ind}; % linear indices of the roi

end