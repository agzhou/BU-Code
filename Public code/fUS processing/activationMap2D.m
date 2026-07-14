
%% Description:
%       Calculate the correlation coefficient for 2D + time data with some
%       stim pattern. Then, process it to get the activation maps

% Processing steps:
%   1. Median filter
%   2. 

% Inputs:
%       volumeData: (z, x, time) data
%       stim: (time x 1) data
%       zt: z-score threshold

% Outputs:
%       r: Pearson's correlation coefficient (z, x)
%       z: z-score with Fisher's transform

%%
function [r, z, activationMap] = activationMap2D(planarData, stim, zt)
    
    % Connected region size
    crs = 9;

    timeDim = 3; % Time is on the 3rd dimension of our 2D + time data

    ps = size(planarData); % planar data's size
    ntp = ps(timeDim); % # of time points
    stim_3D = repmat( permute(stim, [3, 2, 1]), [ps(1), ps(2), 1] ); % Replicate the stim pattern for each pixel

    % Calculate the Pearson's correlation coefficient for every voxel
    r = sum( (planarData - mean(planarData, timeDim)) .* (stim_3D - mean(stim_3D, timeDim)) , timeDim) ...
        ./ sqrt ( sum( abs((planarData - mean(planarData, timeDim))) .^ 2, timeDim) ) ...
        ./ sqrt ( sum( abs((stim_3D - mean(stim_3D, timeDim))) .^ 2, timeDim) );

    % z score using Fisher's transform
    z = fisherTransform(r, ntp);

    % Median filter
    mf_kernel_size = [3, 3];
    r_mf = medfilt2(r, mf_kernel_size);
    % generateTiffStack_multi({testr}, [8.8, 8.8, 8], 'hot', 5) % test
    
    % z scores for the median-filtered correlation coefficient map
    z_mf = fisherTransform(r_mf, ntp);

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

    activationMap_nonconv = r_mf .* z_mf_mask_thresholded;

    % Convolve the activation map for smoothing.
    % These parameters are what Bingxue and Jianbo previously used
    conv_kernel = ones(7, 7);
    conv_kernel(4, 4) = 49;
    conv_kernel = conv_kernel ./ 49;
    activationMap = convn(activationMap_nonconv, conv_kernel, 'same');
    
end


%% Helper functions
function [z] = fisherTransform(r, ntp)
    z = sqrt(ntp - 3)/2 .* log((1 + r) ./ (1 - r));
end