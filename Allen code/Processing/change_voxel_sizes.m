function [data_avg_rs, prereg_params] = change_voxel_sizes(data, P, PData)
    %% Prepare template(s) for atlas registration
    % Create templates for each hemodynamic parameter, averaging across superframes
    data_avg = mean(data, 4);
    
    voxel_size = PData.PDelta .* P.wl; % Voxel size (y, x, z) in meters
    fUS_volume_dimensions_m = [P.Trans.numelements/2 * P.Trans.spacingMm / 1e3, P.Trans.numelements/2 * P.Trans.spacingMm / 1e3, (P.endDepthMM - P.startDepthMM)/1e3]; % Volume size in meters
    fUS_volume_dimensions_voxels = PData.Size; % Volume size in voxels (from the recon PData)
    
    % Adjust the sizes based on the pre-SVD/clutter filtering cropping
    fUS_cropped_volume_dimensions_voxels = size(data_avg);
    fUS_cropped_volume_dimensions_m = fUS_cropped_volume_dimensions_voxels ./ fUS_volume_dimensions_voxels .* fUS_volume_dimensions_m;
    
    % User input for target voxel size (post-interpolation)
    targetVoxelSizePrompt = {'y Target Voxel Size [um]', 'x Target Voxel Size [um]', 'z Target Voxel Size [um]'};
    % targetVoxelSizeDefaults = {'10', '10', '10'};
    targetVoxelSizeDefaults = {'50', '50', '50'};
    targetVoxelSizeUserInput = inputdlg(targetVoxelSizePrompt, 'Input Target Voxel Size', 1, targetVoxelSizeDefaults);
    
    % Store target voxel size inputs and convert to meters
    target_voxel_size(1) = str2double(targetVoxelSizeUserInput{1}) ./ 1e6;
    target_voxel_size(2) = str2double(targetVoxelSizeUserInput{2}) ./ 1e6;
    target_voxel_size(3) = str2double(targetVoxelSizeUserInput{3}) ./ 1e6;
    
    prereg_interp_factor = voxel_size ./ target_voxel_size;
    
    % Resample hemodynamic parameter template maps to the desired voxel size
    data_avg_rs = imresize3(data_avg, 'Scale', prereg_interp_factor, 'Method', 'cubic');
    
    % Store pre-registration parameters
    prereg_params.orig_voxel_size = voxel_size;
    prereg_params.fUS_volume_dimensions_m = fUS_volume_dimensions_m;
    prereg_params.fUS_volume_dimensions_voxels = fUS_volume_dimensions_voxels;
    prereg_params.fUS_cropped_volume_dimensions_voxels = fUS_cropped_volume_dimensions_voxels;
    prereg_params.fUS_cropped_volume_dimensions_m = fUS_cropped_volume_dimensions_m;
    prereg_params.target_voxel_size = target_voxel_size;
    prereg_params.prereg_interp_factor = prereg_interp_factor;
    % prereg_params. = 
    
    % save([savepath, 'prereg_PDI_params_50um.mat'], "data_avg_rs", "prereg_params")
end