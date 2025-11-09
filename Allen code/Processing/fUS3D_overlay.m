% Description:
%       Overlay 3D activation maps onto a template (like CBV index)

%% Load data (manual for now)
% load('J:\08-28-2025 fUS activation map plotting\Trials 1 to 13\activation_maps.mat')
% load('J:\08-28-2025 fUS activation map plotting\Trials 1 to 13 diff proc\activation_maps.mat')
load('H:\Ultrasound data from 7-21-2025\08-28-2025 AZ01 RCA fUS air puff\run 1\fUS proc 11-08-2025 larger region\fUS_avg_templates.mat')
load('H:\Ultrasound data from 7-21-2025\08-28-2025 AZ01 RCA fUS air puff\run 1\fUS proc 11-08-2025 larger region\fUS_activationmaps.mat')
%%
% figure; imagesc(squeeze(max(CBVi_allSF_avg(:, :, :), [], 1) .^ 0.5)'); colormap gray

test = CBVi_allSF_avg .^ 0.5;
mask = test < 0.25;

test(mask) = 0;

figure; imagesc(squeeze(max(test(:, :, :), [], 1) .^ 0.5)'); colormap gray

%%

cmap = colormap_ULM;

% compareVolumes(CBVi_allSF_avg, am_rPDI, 'gray', cmap)


% compareVolumes(CBVi_allSF_avg_rs .^ 0.5, am_rPDI_rs)
compareVolumes(CBVi_allSF_avg_rs .^ 0.5, am_rCBV_rs)
% compareVolumes(CBVi_allSF_avg_rs .^ 0.5, am_rCBFspeed_rs .^ 1)

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
    shift_factor = 0.1;
    % ramp_alphamap = grid_alphamap - round(length(grid_alphamap) ./ 4);
    ramp_alphamap = grid_alphamap - max(grid_alphamap) .* shift_factor;
    ramp_alphamap(ramp_alphamap < 0) = 0;

    % normalize
    max_factor = 1;
    ramp_alphamap(ramp_alphamap > 0) = ramp_alphamap(ramp_alphamap > 0) ./ max(ramp_alphamap) .* max_factor;

    volshow(vol1 .^ 1, 'Colormap', cmap1, ...
            Parent = viewerThresholded, ...
            RenderingStyle = "VolumeRendering", ...
            SpecularReflectance = 0.8, ...
            Alphamap = ramp_alphamap);
            % Alphamap = "linear");   
 
    shift_factor2 = 0.4;
    ramp_alphamap2 = grid_alphamap - max(grid_alphamap) .* shift_factor2;
    ramp_alphamap2(ramp_alphamap2 < 0) = 0;

    % normalize
    max_factor2 = 1;
    ramp_alphamap2(ramp_alphamap2 > 0) = ramp_alphamap2(ramp_alphamap2 > 0) ./ max(ramp_alphamap2) .* max_factor2;

    % volshow(vol2 .^ 1, Parent = viewerThresholded, ...
    %         RenderingStyle = "MaximumIntensityProjection", ...
    %         Colormap = cmap2, ...
    %         Alphamap = "quadratic");
    volshow(vol2 .^ 1, Parent = viewerThresholded, ...
            RenderingStyle = "MaximumIntensityProjection", ...
            Colormap = cmap2, ...
            Alphamap = ramp_alphamap);
            
    % Hard coding for now...
    % viewerThresholded.CameraPosition = [91.4553 -69.5688 -11.5372];
    % viewerThresholded.CameraTarget = [40.5000 40.5000 31];
    % viewerThresholded.CameraUpVector = [0.1305 -0.2423 -0.9614];
end