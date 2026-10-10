%% Description:
%       Overlay a functional (activation) map on a power Doppler background
%       using two stacked axes: the bottom one shows the background, the top
%       one shows the activation map, transparent wherever the mask is false.

% Inputs:
%       bg: (z, x) background image (e.g. PDI.^0.5)
%       amap: (z, x) activation map (e.g. correlation coefficient map)
%       x_mm, z_mm: x and z coordinates [mm] of the pixels
%
% Optional name-value inputs:
%       'Mask': logical (z, x) mask of the pixels to show. Default is amap ~= 0 & ~isnan(amap)
%       'Alpha': opacity of the overlay inside the mask [0, 1]. Default 0.7
%       'BGColormap': colormap of the background. Default 'gray'
%       'AmapColormap': colormap of the activation map. Default 'jet'
%       'BGClim': color limits of the background. Default [] (automatic)
%       'AmapClim': color limits of the activation map. Default [] (automatic, from the masked values)
%       'ColorbarLabel': label for the activation map's colorbar. Default 'r'
%       'Parent': figure (or tab, panel...) to draw in. Default is a new figure
%       'Title': title for the plot. Default ''
%
% Outputs:
%       ax_bg, ax_amap: axes handles of the background and activation map
%       h_amap: image handle of the activation map

%%
function [ax_bg, ax_amap, h_amap] = overlayActivationMap2D(bg, amap, x_mm, z_mm, varargin)

    p = inputParser;
    addParameter(p, 'Mask', [])
    addParameter(p, 'Alpha', 0.7)
    addParameter(p, 'BGColormap', 'gray')
    addParameter(p, 'AmapColormap', 'jet')
    addParameter(p, 'BGClim', [])
    addParameter(p, 'AmapClim', [])
    addParameter(p, 'ColorbarLabel', 'r')
    addParameter(p, 'Parent', [])
    addParameter(p, 'Title', '')
    parse(p, varargin{:})
    o = p.Results;

    assert(isequal(size(bg), size(amap)), 'The background (%s) and activation map (%s) must be the same size', mat2str(size(bg)), mat2str(size(amap)))

    if isempty(o.Mask)
        o.Mask = amap ~= 0 & ~isnan(amap);
    end
    if isempty(o.Parent)
        o.Parent = figure;
    end
    if isempty(o.AmapClim)
        vals = amap(o.Mask);
        if isempty(vals)
            o.AmapClim = [0, 1]; % Nothing to show; avoid an invalid (empty) color limit
        else
            o.AmapClim = [min(vals), max(vals)];
            if o.AmapClim(1) == o.AmapClim(2) % Constant map; clim needs increasing limits
                o.AmapClim = o.AmapClim + [-0.5, 0.5];
            end
        end
    end

    % Bottom axes: background
    ax_bg = axes(o.Parent);
    imagesc(ax_bg, x_mm, z_mm, bg);
    colormap(ax_bg, o.BGColormap)
    if ~isempty(o.BGClim)
        clim(ax_bg, o.BGClim)
    end
    axis(ax_bg, 'image')
    xlabel(ax_bg, 'x [mm]'); ylabel(ax_bg, 'z [mm]')
    title(ax_bg, o.Title)

    % Top axes: activation map, transparent outside the mask
    ax_amap = axes(o.Parent);
    h_amap = imagesc(ax_amap, x_mm, z_mm, amap);
    h_amap.AlphaData = o.Alpha .* double(o.Mask);
    colormap(ax_amap, o.AmapColormap)
    clim(ax_amap, o.AmapClim)
    axis(ax_amap, 'image')
    ax_amap.Color = 'none'; % Transparent axes background
    ax_amap.XTick = []; ax_amap.YTick = []; % The bottom axes show the ticks
    ax_amap.Box = 'off';

    % Colorbar for the activation map. Adding it shrinks ax_amap, so link the positions afterwards
    cb = colorbar(ax_amap, 'Location', 'eastoutside');
    cb.Label.String = o.ColorbarLabel;
    ax_bg.Position = ax_amap.Position;
    linkaxes([ax_bg, ax_amap]) % Keep zoom and pan in sync
    lp = linkprop([ax_bg, ax_amap], {'Position', 'InnerPosition'}); % Keep the axes on top of each other when the figure is resized
    setappdata(ax_bg, 'PositionLink', lp) % The link is deleted if its handle is garbage-collected, so store it with the axes
end
