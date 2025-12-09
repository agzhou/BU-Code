% Move points in dimension dim
% Distance moved (per frame) depends on the inputs:
%   flow_v_mm_s: flow speed [mm/s] in dimension dim
%   frameRate: frame rate
% function newPts = movePoints(pts, dim, flow_v_mm_s, frameRate, vesselDiam, startDepthMM, endDepthMM, xstart, ystart, zstart)
function [newPts, SP] = movePoints(pts, SP)

    SP.dist_per_frame_mm = SP.flow_v_mm_s/SP.frameRate; % move v mm/s, which is (v/fps_target) mm / frame
    SP.dist_per_frame_m = SP.dist_per_frame_mm/1e3;
    
    % dim = 3;        % dimension to change (x, y, z) -> (1, 2, 3) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dim = 2;
    
    % Sort points to determine which are past the boundary after moving them
    temp_pts = sortPts(pts, SP.dim);
    temp_pts(:, SP.dim) = temp_pts(:, SP.dim) + SP.dist_per_frame_m;
    % figure; scatter3(temp_pts(:, 1), temp_pts(:, 2), temp_pts(:, 3), '.'); axis square
    
    %
    % SP.bound = [SP.vesselDiam, ...
    %          SP.vesselDiam, ...
    %          SP.endDepthMM / 1e3]; % x, y, z boundary in m
    SP.bound = [SP.xstart + SP.vesselDiam/2, ...
             SP.ystart + SP.vesselDiam/2, ...
             SP.zstart + SP.vesselLength/2]; % x, y, z (positive) boundary in m --> assumes points flow in positive direction in dim
    mask_past_boundary = temp_pts(:, SP.dim) > SP.bound(SP.dim);
    
    % % Test
    % test = temp_pts(mask_past_boundary, :);
    % figure; scatter3(test(:, 1), test(:, 2), test(:, 3), '.'); axis square
    
    % Replacement Params
    RP = SP;
    % RP.vesselLength = max(temp_pts(:, SP.dim)) - SP.bound(SP.dim); % Get the approximate width of points that need replacing
    RP.vesselLength = SP.dist_per_frame_m;
    
    which_dims_need_extra_shift = zeros(3, 1);
    which_dims_need_extra_shift(SP.dim) = 1;

    %%%% This only works for z movement I'm pretty sure %%%%
    % Update the center of the replacement points to be in the right place
    RP.xstart = SP.xstart - which_dims_need_extra_shift(1) * (SP.vesselDiam/2 - RP.vesselLength/2);
    RP.ystart = SP.ystart - which_dims_need_extra_shift(2) * (SP.vesselDiam/2 - RP.vesselLength/2);
    RP.zstart = SP.zstart - which_dims_need_extra_shift(3) * (SP.vesselLength/2 - RP.vesselLength/2);
    % [replacePoints, ~] = genRandomPts3D_cyl(replaceDiam, replaceLength, startDepthMM/1e3, xstart, ystart, zstart);
    [replacePoints, ~] = genRandomPts3D_cyl(RP);
    
    % figure; scatter3(replacePoints(:, 1), replacePoints(:, 2), replacePoints(:, 3), '.'); axis square
    
    newPts = [replacePoints; temp_pts(~mask_past_boundary, :)];
    % figure; scatter3(newPts(:, 1), newPts(:, 2), newPts(:, 3), '.'); axis square
end