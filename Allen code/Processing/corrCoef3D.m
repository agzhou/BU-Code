%% Description:
%       Calculate the correlation coefficient for 3D + time data with some
%       stim pattern

% Inputs:
%       volumeData: (y, x, z, time) data
%       stim: (time x 1) data

% Outputs:
%       r: Pearson's correlation coefficient (y, x, z)
%       z: z-score with Fisher's transform

%%
function [r, z] = corrCoef3D(volumeData, stim)
    vs = size(volumeData); % volume data's size
    stim_4D = repmat( permute(stim, [4, 3, 2, 1]), [vs(1), vs(2), vs(3), 1] ); % Replicate the stim pattern for each voxel

    % Calculate the Pearson's correlation coefficient for every voxel
    r = sum( (volumeData - mean(volumeData, 4)) .* (stim_4D - mean(stim_4D, 4)) , 4) ...
        ./ sqrt ( sum( abs((volumeData - mean(volumeData, 4))) .^ 2, 4) ) ...
        ./ sqrt ( sum( abs((stim_4D - mean(stim_4D, 4))) .^ 2, 4) );

    % z score using Fisher's transform
    z = sqrt(size(volumeData, 4) - 3)/2 .* log((1 + r) ./ (1 - r));
end