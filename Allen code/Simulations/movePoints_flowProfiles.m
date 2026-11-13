% Move points in dimension dim, according to some flow profile
% Distance moved (per frame) depends on the inputs:
%   pts: (# points, 4) matrix. First 3 columns are coordinates, 4th is reflectivity.
%   SP.flow_vgp_mm_s: group flow speed [mm/s] in dimension dim
%   flowProfile: function that returns the displacement of a scatterer matrix
%                   per frame. Inputs: pts
%   SP.frameRate: frame rate
function [newPts, SP] = movePoints_flowProfiles(pts, flowProfile, SP)

    % Calculate per-frame displacement for each point
    % SP.dist_per_frame_mm = SP.flow_vgp_mm_s/SP.frameRate; % move v mm/s, which is (v/fps_target) mm / frame
    % SP.dist_per_frame_m = SP.dist_per_frame_mm/1e3;
    frame_disp = flowProfile(pts);

    dim = SP.dim;        % dimension to change (x, y, z) -> (1, 2, 3) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Sort points to determine which are past the boundary after moving them
    % temp_pts = sortPts(pts, dim);
    % temp_pts(:, dim) = temp_pts(:, dim) + frame_disp(order///);

    temp_pts = pts;
    temp_pts(:, dim) = temp_pts(:, dim) + frame_disp;
    temp_pts = sortPts(temp_pts, dim);

    % figure; scatter3(temp_pts(:, 1), temp_pts(:, 2), temp_pts(:, 3), '.'); axis square
    % figure; scatter3(pts(:, 1), pts(:, 2), pts(:, 3), '.'); axis square

    %
    % SP.bound = [SP.vesselDiam, ...
    %          SP.vesselDiam, ...
    %          SP.endDepthMM / 1e3]; % x, y, z boundary in m
    SP.bound = [SP.vesselCenter(1) + SP.vesselDiam/2, ...
             SP.vesselCenter(2) + SP.vesselDiam/2, ...
             SP.vesselCenter(3) + SP.vesselLength/2]; % x, y, z (positive) boundary in m --> assumes points flow in positive direction in dim
    mask_past_boundary = temp_pts(:, dim) > SP.bound(dim);
    
    % % Test
    % test = temp_pts(mask_past_boundary, :);
    % figure; scatter3(test(:, 1), test(:, 2), test(:, 3), '.'); axis square
    

    % % Replacement Params
    % RP = SP;
    % % RP.vesselLength = max(temp_pts(:, dim)) - SP.bound(dim); % Get the approximate width of points that need replacing
    % RP.vesselLength = SP.dist_per_frame_m;
    % 
    % which_dims_need_extra_shift = zeros(3, 1);
    % which_dims_need_extra_shift(dim) = 1;
    % 
    % %%%% This only works for z movement I'm pretty sure %%%%
    % % Update the center of the replacement points to be in the right place
    % RP.vesselCenter(1) = SP.vesselCenter(1) - which_dims_need_extra_shift(1) * (SP.vesselDiam/2 - RP.vesselLength/2);
    % RP.vesselCenter(2) = SP.vesselCenter(2) - which_dims_need_extra_shift(2) * (SP.vesselDiam/2 - RP.vesselLength/2);
    % RP.vesselCenter(3) = SP.vesselCenter(3) - which_dims_need_extra_shift(3) * (SP.vesselLength/2 - RP.vesselLength/2);
    % % [replacePoints, ~] = genRandomPts3D_cyl(replaceDiam, replaceLength, startDepthMM/1e3, xstart, ystart, zstart);
    % [replacePoints, ~] = genRandomPts3D_cyl(RP);
    % 
    replacePoints = temp_pts(mask_past_boundary, :);
    replacePoints(:, dim) = replacePoints(:, dim) - SP.vesselLength;
    % figure; scatter3(replacePoints(:, 1), replacePoints(:, 2), replacePoints(:, 3), '.'); axis square
    
    newPts = [replacePoints; temp_pts(~mask_past_boundary, :)];
    % figure; scatter3(newPts(:, 1), newPts(:, 2), newPts(:, 3), '.'); axis square
end