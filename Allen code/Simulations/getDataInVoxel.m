
% Get the points in the "vessel"/data within a voxel (a struct with the
% center, size defined)

function data_in_voxel = getDataInVoxel(data, voxel)
    temp_voxel_mask = abs(data(:, 1) - voxel.center(1)) < voxel.size(1)/2 & abs(data(:, 2) - voxel.center(2)) < voxel.size(2)/2 & abs(data(:, 3) - voxel.center(3)) < voxel.size(3)/2;
    data_in_voxel = data(temp_voxel_mask, :);



end