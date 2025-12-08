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
    
    temp_pts = sortPts(pts, SP.dim);
    temp_pts(:, SP.dim) = temp_pts(:, SP.dim) + SP.dist_per_frame_m;
    % figure; scatter3(temp_pts(:, 1), temp_pts(:, 2), temp_pts(:, 3), '.'); axis square
    
    %
    SP.bound = [SP.vesselDiam, ...
             SP.vesselDiam, ...
             SP.endDepthMM / 1e3]; % x, y, z boundary in m
    mask_past_boundary = temp_pts(:, SP.dim) > SP.bound(SP.dim);
    
    % % Test
    % test = temp_pts(mask_past_boundary, :);
    % figure; scatter3(test(:, 1), test(:, 2), test(:, 3), '.'); axis square
    
    % Vertical replacement %%%%%%%%
    % replacement params
    RP = SP;
    RP.vesselLength = max(temp_pts(:, SP.dim)) - SP.bound(SP.dim);
    % [replacePoints, ~] = genRandomPts3D_cyl(replaceDiam, replaceLength, startDepthMM/1e3, xstart, ystart, zstart);
    [replacePoints, ~] = genRandomPts3D_cyl(RP);
    
    % figure; scatter3(replacePoints(:, 1), replacePoints(:, 2), replacePoints(:, 3))
    
    newPts = [replacePoints; temp_pts(~mask_past_boundary, :)];
    % figure; scatter3(newPts(:, 1), newPts(:, 2), newPts(:, 3), '.'); axis square
end