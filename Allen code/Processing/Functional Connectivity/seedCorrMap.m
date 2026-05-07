%% Description: calculate a seed correlation map for 4D (3D space + time) data

% Inputs:
%   - seed: a seed timecourse (# time points x 1)
%   - volumeData: volumetric data over time (# x voxels, # y voxels, # z voxels, # time points) 

function [r, z] = seedCorrMap(seed, volumeData)
    [r, z] = corrCoef3D(volumeData, seed);

end