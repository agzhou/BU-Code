% Description:
%       Overlay 3D activation maps onto a template (like CBV index)

%%
% figure; imagesc(squeeze(max(CBVi_allSF_avg(:, :, :), [], 1) .^ 0.5)'); colormap gray

test = CBVi_allSF_avg .^ 0.5;
mask = test < 0.25;

test(mask) = 0;

figure; imagesc(squeeze(max(test(:, :, :), [], 1) .^ 0.5)'); colormap gray

%%

cmap = colormap_ULM;

% compareVolumes(CBVi_allSF_avg, am_rPDI, 'gray', cmap)
compareVolumes(CBVi_allSF_avg .^ 0.5, am_rPDI)
% compareVolumes(CBVi_allSF_avg .^ 1, am_rPDI)
% compareVolumes(test, am_rPDI)

% compareVolumes(CBVi_allSF_avg, am_rCBV)
% compareVolumes(CBVi_allSF_avg, am_rCBFspeed)

% grid_alphamap = linspace(0, 1, 256)';
% % ramp_alphamap = grid_alphamap - round(length(grid_alphamap) ./ 4);
% ramp_alphamap = grid_alphamap - max(grid_alphamap) ./ 4;
% ramp_alphamap(ramp_alphamap < 0) = 0;
%%
function compareVolumes(vol1, vol2, varargin) % Can change this so it has a cell array input and goes through more than 2 volumes

    % Set default colormaps
    cmap1 = gray;
    cmap2 = colormap_ULM;
    if nargin > 2
        cmap1 = varargin{1};
        if nargin > 3
            cmap2 = varargin{2};
        end
    end

    viewerThresholded = viewer3d(BackgroundColor = "white", ...
                                 BackgroundGradient="off", ...
                                 RenderingQuality = "high", ...
                                 Denoising = "off", ...
                                 SpatialUnits = "pixels", ...
                                 Box = "on");
    % volshow(vol1 .^ 1, 'Colormap', cmap1, ...
    %         Parent = viewerThresholded, ...
    %         RenderingStyle = "MaximumIntensityProjection", ...
    %         Alphamap = "linear");   

    grid_alphamap = linspace(0, 1, 256)';
    % ramp_alphamap = grid_alphamap - round(length(grid_alphamap) ./ 4);
    ramp_alphamap = grid_alphamap - max(grid_alphamap) .* 0.1;
    ramp_alphamap(ramp_alphamap < 0) = 0;

    volshow(vol1 .^ 1, 'Colormap', cmap1, ...
            Parent = viewerThresholded, ...
            RenderingStyle = "VolumeRendering", ...
            SpecularReflectance = 0.8, ...
            Alphamap = ramp_alphamap);
            % Alphamap = "linear");   

    % volshow(vol1 .^ 1, 'Colormap', cmap1, ...
    %         Parent = viewerThresholded, ...
    %         RenderingStyle = "Isosurface", ...
    %         IsosurfaceValue = 0.10, ...
    %         IsosurfaceAlpha = 0.8, ...
    %         Alphamap = "linear");   

    volshow(vol2 .^ 1, Parent = viewerThresholded, ...
            RenderingStyle = "MaximumIntensityProjection", ...
            Colormap = cmap2, ...
            Alphamap = "linear");

end