function movePointsZ3D_noreplace

    if evalin('base','exist(''Media'',''var'')')
        Media = evalin('base', 'Media');
    else
        disp('Media object not found in workplace.');
        return
    end

    if evalin('base','exist(''fps_target'',''var'')')
        fps_target = evalin('base', 'fps_target');
    else
        disp('fps_target not found in workplace.');
        return
    end

    if evalin('base','exist(''wl'',''var'')')
        wl = evalin('base', 'wl');
    else
        disp('wl not found in workplace.');
        return
    end

    if evalin('base','exist(''Trans'',''var'')')
        Trans = evalin('base', 'Trans');
    else
        disp('Trans not found in workplace.');
        return
    end

%     evalin('base', '');

%     disp('moving')
    
    
    %% Sort Media points based on one dimension

%     flow_v_mm_s = 30;
%     flow_v_mm_s = 3000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test to see if replacement works

%     assignin('base', 'flow_v_mm_s', flow_v_mm_s); 
    flow_v_mm_s = evalin('base', 'flow_v_mm_s'); % change 10/14/24
    dist_per_frame_mm = flow_v_mm_s/fps_target; % move v mm/s, which is (v/fps_target) mm / frame
%     wl_mm = 1540 / 13.8889 / 1e3; %
    wl_mm = wl * 1e3;
    dist_per_frame_wl = dist_per_frame_mm/wl_mm;

    dim = 3;        % dimension to change (x, y, z) -> (1, 2, 3)
    % dim = 2;

%     [temp_Media_sorted_dim, ind_Media_sorted] = sort(Media.MP(:, dim), 1);
%     temp_Media_sorted = Media.MP(ind_Media_sorted, :);
%     %%
%     temp_Media = temp_Media_sorted;
%     %%
%     temp_Media(:, dim) = temp_Media(:, dim) + dist_per_frame_wl;
%     disp('test1')
    temp_Media = sortMedia(Media, dim);
    temp_Media(:, dim) = temp_Media(:, dim) + dist_per_frame_wl;
    %%
    % bound = [(Trans.numelements / 2) / 2 * Trans.spacing, ...
    %          (Trans.numelements / 2) / 2 * Trans.spacing, ...
    %          endDepth]; % x, y, z boundary in wavelengths
    % mask_past_boundary = temp_Media(:, dim) > bound(dim);

%     figure; scatter(Media.MP(:, 1), Media.MP(:, 3), '.') %%

    % replacePoints = [];
    Media.MP = temp_Media;
    % Media.MP = [replacePoints; temp_Media(~mask_past_boundary, :)];

%     figure; scatter(Media.MP(:, 1), Media.MP(:, 3), '.') %%

%     Media.MP(:, dim) =  % Modify position of all media points
%     disp(num2str(max(Media.MP(:, 3))))
%%
%     disp(Media.MP(1, :))
    assignin('base', 'Media', Media); % put the updated Media into the workspace
% end

return
