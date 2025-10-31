% Description:
%   Overlay two volumes:
%       1) Ultrasound template map (e.g., registered PDI averaged across superframes)
%       2) Region masks of interest from an annotated atlas

function compareUStoAtlasROIs(US_template, region_masks)


    % US_template = imresize3(US_template, 1/2, 'Method', 'cubic');
    % region_masks = imresize3(region_masks, 1/2, 'Method', 'cubic');

    US_gamma = 0.5;
    
    % % Go through all the ROI masks and union them together
    % allMasks = region_masks{1};
    % for ind = 2:length(region_masks)
    %     allMasks = allMasks | region_masks{ind};
    % end

    % Go through all the ROI masks and assign a "color" to each of them
    % allMasks = uint8(region_masks{1});
    % for ind = 2:length(region_masks)
    %     allMasks = allMasks + uint8(region_masks{ind}) .* ind;
    % end

    % allMasks = allMasks + 10;

    viewerThresholded = viewer3d(BackgroundColor = "white", ...
                                 BackgroundGradient="off", ...
                                 RenderingQuality = "high", ...
                                 Denoising = "off", ...
                                 SpatialUnits = "pixels", ...
                                 Box = "on");

    grid_alphamap = linspace(0, 1, 256)';
    shift_factor = 0.1;
    % ramp_alphamap = grid_alphamap - round(length(grid_alphamap) ./ 4);
    ramp_alphamap = grid_alphamap - max(grid_alphamap) .* shift_factor;
    ramp_alphamap(ramp_alphamap < 0) = 0;

    % normalize
    max_factor = 1;
    ramp_alphamap(ramp_alphamap > 0) = ramp_alphamap(ramp_alphamap > 0) ./ max(ramp_alphamap) .* max_factor;

    % Plot the ultrasound template
    cmap1 = gray;
    volshow(US_template .^ US_gamma, 'Colormap', cmap1, ...
            Parent = viewerThresholded, ...
            RenderingStyle = "MaximumIntensityProjection", ...
            SpecularReflectance = 0.8, ...
            Alphamap = ramp_alphamap);
            % Alphamap = "linear");

    % % Mask alphamap: make everything 1 except when the input voxel is 0
    % mask_alphamap = ones(size(grid_alphamap));
    % mask_alphamap(1) = 0;
    % volshow(allMasks .^ 1, Parent = viewerThresholded, RenderingStyle = "VolumeRendering", ...
    %     Colormap=jet, Alphamap = "linear");

    for ind = 1:length(region_masks)
        volshow(region_masks{ind} .^ 1, Parent = viewerThresholded, RenderingStyle = "Isosurface", ...
                Colormap=rand(1, 3), Alphamap = 1);
    end
    % volshow(US_template .^ US_gamma, Parent=viewerThresholded, RenderingStyle = "Isosurface", ...
    %     Colormap=[1 0 1], Alphamap=1);
    % volshow(allMasks .^ 1, Parent = viewerThresholded, RenderingStyle = "Isosurface", ...
    %     Colormap=[1 0 0], Alphamap = 0.5);

end