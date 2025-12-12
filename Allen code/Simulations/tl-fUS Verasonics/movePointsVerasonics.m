% Move points in dimension dim --> for use in Verasonics simulations

% Distance moved (per frame) depends on the inputs:
%   flow_v_mm_s: flow speed [mm/s] in dimension dim
%   frameRate: frame rate
function movePointsVerasonics

    % Get the data from the workspace
    if evalin('base','exist(''Media'',''var'')')
        Media = evalin('base', 'Media');
    else
        disp('Media object not found in workplace.');
        return
    end

    if evalin('base','exist(''Mcr_SP'',''var'')')
        Mcr_SP = evalin('base', 'Mcr_SP');
        Mcr_SP.testflag = true;
    else
        disp('Mcr_SP not found in workplace.');
        return
    end

    if evalin('base','exist(''vrp'',''var'')')
        vrp = evalin('base', 'vrp');
    else
        disp('vrp not found in workplace.');
        return
    end

    % Load the un-rotated vessel
    if evalin('base','exist(''cyl_vessel'',''var'')')
        cyl_vessel = evalin('base', 'cyl_vessel');
    else
        disp('cyl_vessel not found in workplace.');
        return
    end
    

    Mcr_SP.dist_per_frame_mm = Mcr_SP.flow_v_mm_s/Mcr_SP.frameRate; % move v mm/s, which is (v/fps_target) mm / frame
    Mcr_SP.dist_per_frame_m = Mcr_SP.dist_per_frame_mm/1e3;
    
    % Sort points to determine which are past the boundary after moving them
    pts = cyl_vessel; % Start with the un-rotated vessel with point locations in meters
    temp_pts = sortPts(pts, Mcr_SP.dim);
    temp_pts(:, Mcr_SP.dim) = temp_pts(:, Mcr_SP.dim) + Mcr_SP.dist_per_frame_m;
    % figure; scatter3(temp_pts(:, 1), temp_pts(:, 2), temp_pts(:, 3), '.'); axis square
    
    %
    % SP.bound = [SP.vesselDiam, ...
    %          SP.vesselDiam, ...
    %          SP.endDepthMM / 1e3]; % x, y, z boundary in m
    Mcr_SP.bound = [Mcr_SP.xstart + Mcr_SP.vesselDiam/2, ...
             Mcr_SP.ystart + Mcr_SP.vesselDiam/2, ...
             Mcr_SP.zstart + Mcr_SP.vesselLength/2]; % x, y, z (positive) boundary in m --> assumes points flow in positive direction in dim
    mask_past_boundary = temp_pts(:, Mcr_SP.dim) > Mcr_SP.bound(Mcr_SP.dim);
    
    % % Test
    % test = temp_pts(mask_past_boundary, :);
    % figure; scatter3(test(:, 1), test(:, 2), test(:, 3), '.'); axis square
    
    % Replacement Params
    RP = Mcr_SP;
    % RP.vesselLength = max(temp_pts(:, SP.dim)) - SP.bound(SP.dim); % Get the approximate width of points that need replacing
    RP.vesselLength = Mcr_SP.dist_per_frame_m;
    
    which_dims_need_extra_shift = zeros(3, 1);
    which_dims_need_extra_shift(Mcr_SP.dim) = 1;

    %%%% This only works for z movement I'm pretty sure %%%%
    % Update the center of the replacement points to be in the right place
    RP.xstart = Mcr_SP.xstart - which_dims_need_extra_shift(1) * (Mcr_SP.vesselDiam/2 - RP.vesselLength/2);
    RP.ystart = Mcr_SP.ystart - which_dims_need_extra_shift(2) * (Mcr_SP.vesselDiam/2 - RP.vesselLength/2);
    RP.zstart = Mcr_SP.zstart - which_dims_need_extra_shift(3) * (Mcr_SP.vesselLength/2 - RP.vesselLength/2);
    % [replacePoints, ~] = genRandomPts3D_cyl(replaceDiam, replaceLength, startDepthMM/1e3, xstart, ystart, zstart);
    [replacePoints, ~] = genRandomPts3D_cyl(RP);
    
    % figure; scatter3(replacePoints(:, 1), replacePoints(:, 2), replacePoints(:, 3), '.'); axis square
    
    newPts = [replacePoints; temp_pts(~mask_past_boundary, :)];
    % figure; scatter3(newPts(:, 1), newPts(:, 2), newPts(:, 3), '.'); axis square



    % Re-rotate
    Media.MP = rotateVessel(newPts, vrp.xa, vrp.ya, vrp.za, Mcr_SP);
    % Media.MP = [[0, 0, 3] .* 1e-3, Mcr_SP.scatterReflectivity]; %%%%%%%%%%%%%%%
    Media.MP(:, 1:3) = Media.MP(:, 1:3) ./ Mcr_SP.wl; % CONVERT TO WAVELENGTHS!!!!!!

    % Just for testing
%     Media.MP = repmat(Media.MP, 2, 1);

    assignin('base', 'Media', Media); % put the updated Media into the workspace
    
    assignin('base', 'Mcr_SP', Mcr_SP);

    % Update the non-rotated points for the next time
    cyl_vessel = newPts;
    assignin('base', 'cyl_vessel', cyl_vessel)
end