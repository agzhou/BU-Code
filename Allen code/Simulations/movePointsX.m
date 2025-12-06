function movePointsX

    if evalin('base','exist(''Media'',''var'')')
        Media = evalin('base','Media');
    else
        disp('Media object not found in workplace.');
        return
    end
    
    dif_mm = 30/fps_target; % move 30 mm/s, which is .03 mm / ms (1 ms per frame for 1 kHz)
    wl_mm = 1540 / 15.625 / 1e3; % for 1D array L22-14v
    dif = dif_mm/wl_mm;
    
    
    dim = 1;
    % dim = 2;
        Media.MP(:, dim) = Media.MP(:, dim) + dif; % Modify x position of all media points
    
    %     Media.MP
        assignin('base', 'Media', Media);
    % end

return
