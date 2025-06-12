
%% Description:
%       Calculate the correlation coefficient for 3D + time data with some
%       stim pattern. Then, process it to get the activation maps

% Processing steps:
%   1. Median filter
%   2. 

% Inputs:
%       volumeData: (y, x, z, time) data
%       stim: (time x 1) data
%       zt: z-score threshold

% Outputs:
%       r: Pearson's correlation coefficient (y, x, z)
%       z: z-score with Fisher's transform

%%
function [r, z, activationMap] = activationMap3D(volumeData, stim, zt)
    
    % Connected region size
    crs = 9;

    vs = size(volumeData); % volume data's size
    ntp = vs(4); % # of time points
    stim_4D = repmat( permute(stim, [4, 3, 2, 1]), [vs(1), vs(2), vs(3), 1] ); % Replicate the stim pattern for each voxel

    % Calculate the Pearson's correlation coefficient for every voxel
    r = sum( (volumeData - mean(volumeData, 4)) .* (stim_4D - mean(stim_4D, 4)) , 4) ...
        ./ sqrt ( sum( (volumeData - mean(volumeData, 4)) .^ 2, 4) ) ...
        ./ sqrt ( sum( (stim_4D - mean(stim_4D, 4)) .^ 2, 4) );

    % z score using Fisher's transform
    z = fisherTransform3D(r, ntp);

    % Median filter
    mf_kernel_size = [3, 3, 3];
    r_mf = medfilt3(r, mf_kernel_size);
    % generateTiffStack_multi({testr}, [8.8, 8.8, 8], 'hot', 5) % test
    
    % z scores for the median-filtered correlation coefficient map
    z_mf = fisherTransform3D(r_mf, ntp);

    % Create a mask where the z-score exceeds a threshold, and where the
    % resulting mask has connectivity above some value "crs"
    z_mf_mask = z_mf > zt;
    bwconn_z_mf_mask = bwconncomp(z_mf_mask);
    z_mf_mask_thresholded = z_mf_mask;
    for i = 1:length(bwconn_z_mf_mask.PixelIdxList)
        tempPixList = bwconn_z_mf_mask.PixelIdxList{i};
        if length(tempPixList) < crs
            z_mf_mask_thresholded(tempPixList) = 0;
        end
    end

    activationMap = r_mf .* z_mf_mask_thresholded;

    % Convolve the activation map
    conv_kernel = ones(7, 7, 7);
    conv_kernel(4, 4, 4) = 49;
    conv_kernel = conv_kernel ./ 49;
    activationMap_covn = convn(activationMap, conv_kernel, 'same');
    
end


%% Helper functions
function [z] = fisherTransform3D(r, ntp)
    z = sqrt(ntp - 3)/2 .* log((1 + r) ./ (1 - r));
end