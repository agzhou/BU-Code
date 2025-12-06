function movePointsZ2D

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

    endDepth = evalin('base', 'endDepth');
    vesselX = evalin('base', 'vesselX');
    vesselZ = evalin('base', 'vesselZ');

%     evalin('base', '');

%     disp('moving')
    
    
    %% Sort Media points based on one dimension

    flow_v_mm_s = 30;
%     flow_v_mm_s = 3000; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test to see if replacement works

    assignin('base', 'flow_v_mm_s', flow_v_mm_s); 
    dist_per_frame_mm = flow_v_mm_s/fps_target; % move v mm/s, which is (v/fps_target) mm / frame
    wl_mm = 1540 / 15.625 / 1e3; % for 1D array L22-14v
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

    %%%%%%%%%%% temp for test
%     endDepthMM = 15;
%     endDepth = endDepthMM / 1e3 / wl;
%     % vertical
%     vesselX = 30e-6;    % x dimension
%     vesselZ = endDepthMM * 1e-3; % z dimension
%     Trans.numelements = 128;
%     Trans.spacing = 0.1 / 1e3 / wl;

    %%%%%%%%


    bound = [Trans.numelements*Trans.spacing / 2, 0, endDepth]; % x, y, z boundary in wavelengths
    mask_past_boundary = temp_Media(:, dim) > bound(dim);

    % Vertical
    replaceX = vesselX;
    replaceZ = max(temp_Media(:, dim)) * wl - bound(dim)*wl;
%     disp('test2')
    replacePoints = randomPts2D_func(replaceX, replaceZ, wl); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     replacePoints = [];
%     disp('test3')


%     figure; scatter(Media.MP(:, 1), Media.MP(:, 3), '.') %%


    Media.MP = [replacePoints; temp_Media(~mask_past_boundary, :)];

%     figure; scatter(Media.MP(:, 1), Media.MP(:, 3), '.') %%

%     Media.MP(:, dim) =  % Modify position of all media points
%     disp(num2str(max(Media.MP(:, 3))))
%%
%     disp(Media.MP(1, :))
    assignin('base', 'Media', Media); % put the updated Media into the workspace
% end

return
