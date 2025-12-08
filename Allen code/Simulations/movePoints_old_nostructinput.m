% Move points in dimension dim
% Distance moved (per frame) depends on the inputs:
%   flow_v_mm_s: flow speed [mm/s] in dimension dim
%   frameRate: frame rate
function newPts = movePoints(pts, dim, flow_v_mm_s, frameRate, vesselDiam, startDepthMM, endDepthMM, xstart, ystart, zstart)
    dist_per_frame_mm = flow_v_mm_s/frameRate; % move v mm/s, which is (v/fps_target) mm / frame
    dist_per_frame_m = dist_per_frame_mm/1e3;
    
    dim = 3;        % dimension to change (x, y, z) -> (1, 2, 3) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dim = 2;
    
    temp_pts = sortPts(pts, dim);
    temp_pts(:, dim) = temp_pts(:, dim) + dist_per_frame_m;
    % figure; scatter3(temp_pts(:, 1), temp_pts(:, 2), temp_pts(:, 3), '.'); axis square
    
    %
    bound = [vesselDiam, ...
             vesselDiam, ...
             endDepthMM / 1e3]; % x, y, z boundary in m
    mask_past_boundary = temp_pts(:, dim) > bound(dim);
    
    % % Test
    % test = temp_pts(mask_past_boundary, :);
    % figure; scatter3(test(:, 1), test(:, 2), test(:, 3), '.'); axis square
    
    % Vertical replacement %%%%%%%%
    replaceDiam = vesselDiam;
    replaceLength = max(temp_pts(:, dim)) - bound(dim);
    replacePoints = genRandomPts3D_cyl(replaceDiam, replaceLength, startDepthMM/1e3, xstart, ystart, zstart);
    
    % figure; scatter3(replacePoints(:, 1), replacePoints(:, 2), replacePoints(:, 3))
    
    newPts = [replacePoints; temp_pts(~mask_past_boundary, :)];
    % figure; scatter3(newPts(:, 1), newPts(:, 2), newPts(:, 3), '.'); axis square
end