% Description:
%       Overlay 3D activation maps onto a template (like CBV index)

%% Load data (manual for now)
% load('J:\08-28-2025 fUS activation map plotting\Trials 1 to 13\activation_maps.mat')
% load('J:\08-28-2025 fUS activation map plotting\Trials 1 to 13 diff proc\activation_maps.mat')
% load('H:\Ultrasound data from 7-21-2025\08-28-2025 AZ01 RCA fUS air puff\run 1\fUS proc 11-08-2025 larger region\fUS_avg_templates.mat')
% load('H:\Ultrasound data from 7-21-2025\08-28-2025 AZ01 RCA fUS air puff\run 1\fUS proc 11-08-2025 larger region\fUS_activationmaps.mat')

load('J:\08-28-2025 fUS activation map plotting\fUS proc 11-08-2025 larger region\Trials 1 to 13\fUS_avg_templates.mat')
load('J:\08-28-2025 fUS activation map plotting\fUS proc 11-08-2025 larger region\Trials 1 to 13\fUS_activationmaps.mat')

%%
% figure; imagesc(squeeze(max(CBVi_allSF_avg(:, :, :), [], 1) .^ 0.5)'); colormap gray

test = CBVi_allSF_avg .^ 0.5;
mask = test < 0.25;

test(mask) = 0;

figure; imagesc(squeeze(max(test(:, :, :), [], 1) .^ 0.5)'); colormap gray

%%

cmap = colormap_ULM;

% compareVolumes(CBVi_allSF_avg, am_rPDI, 'gray', cmap)


% v_rPDI = compareVolumes(CBVi_allSF_avg_rs .^ 0.5, am_rPDI_rs, gray, cmap, 0.5);
% v_rCBV = compareVolumes(CBVi_allSF_avg_rs .^ 0.5, am_rCBV_rs, gray, cmap, 0.5);
% v_rCBFspeed = compareVolumes(CBVi_allSF_avg_rs .^ 0.5, am_rCBFspeed_rs, gray, cmap, 0.0);

v_rPDI = compareVolumes(CBVi_allSF_avg_rs_reg .^ 0.5, am_rPDI_rs_reg, gray, cmap, 0.2, 0.5);
v_rCBV = compareVolumes(CBVi_allSF_avg_rs_reg .^ 0.5, am_rCBV_rs_reg, gray, cmap, 0.2, 0.5);
v_rCBFspeed = compareVolumes(CBVi_allSF_avg_rs_reg .^ 0.5, am_rCBFspeed_rs_reg, gray, cmap, 0.2, 0.0);



% Top view
% camPos = [100.7983 118.8956 -302.8276];
% camTarget = [109.3857 129.9169 80.6888];
% camUpVec = [-0.0092 0.7996 -0.6005];
% cz = 1.4;

% Side view
% camPos = [85.0166 -241.9843 159.0635];
% camTarget = [114.5000 132.5000 80.5000];
% camUpVec = [-0.0650 -0.7294 -0.6810];
camPos = [76.5936 -238.5990 170.6589];
camTarget = [114.5000 132.5000 80.5000];
camUpVec = [-0.0772 -0.7490 -0.6580];
cz = 1.7;

% Set the camera properties
v_rPDI.CameraZoom = cz;
v_rCBV.CameraZoom = cz;
v_rCBFspeed.CameraZoom = cz;

v_rPDI.CameraPosition = camPos;
v_rCBV.CameraPosition = camPos;
v_rCBFspeed.CameraPosition = camPos;
v_rPDI.CameraTarget = camTarget;
v_rCBV.CameraTarget = camTarget;
v_rCBFspeed.CameraTarget = camTarget;
v_rPDI.CameraUpVector = camUpVec;
v_rCBV.CameraUpVector = camUpVec;
v_rCBFspeed.CameraUpVector = camUpVec;
% compareVolumes(PDI_allSF_avg_rs .^ 0.5, am_rPDI_rs)
% compareVolumes(PDI_allSF_avg_rs .^ 0.5, am_rCBV_rs)
% compareVolumes(PDI_allSF_avg_rs .^ 0.5, am_rCBFspeed_rs .^ 1)

%% Plot the barrel region
v_barrelfield = compareVolumes(CBVi_allSF_avg_rs_reg .^ 0.5, region_masks_50um{1}, gray, cmap, 0.2, 0.0);
v_barrelfield.CameraZoom = cz;
v_barrelfield.CameraPosition = camPos;
v_barrelfield.CameraTarget = camTarget;
v_barrelfield.CameraUpVector = camUpVec;

%%
function [viewerThresholded] = compareVolumes(vol1, vol2, varargin) % Can change this so it has a cell array input and goes through more than 2 volumes

    % Set default colormaps
    cmap1 = gray;
    cmap2 = colormap_ULM;

    % Default shift factors
    shift_factor1 = 0.1;
    shift_factor2 = 0.5;

    if nargin > 2
        cmap1 = varargin{1};
        if nargin > 3
            cmap2 = varargin{2};
            if nargin > 4
                shift_factor1 = varargin{3};
                    if nargin > 5
                        shift_factor2 = varargin{4};
                    end
            end
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
    % shift_factor = 0.1;
    % shift_factor1 = 0.4;
    % ramp_alphamap = grid_alphamap - round(length(grid_alphamap) ./ 4);
    ramp_alphamap = grid_alphamap - max(grid_alphamap) .* shift_factor1;
    ramp_alphamap(ramp_alphamap < 0) = 0;

    % normalize
    max_factor1 = 1;
    % max_factor1 = 0.8;
    ramp_alphamap(ramp_alphamap > 0) = ramp_alphamap(ramp_alphamap > 0) ./ max(ramp_alphamap) .* max_factor1;

    volshow(vol1 .^ 1, 'Colormap', cmap1, ...
            Parent = viewerThresholded, ...
            RenderingStyle = "VolumeRendering", ...
            SpecularReflectance = 0.8, ...
            Alphamap = ramp_alphamap);
            % Alphamap = "linear");   
 

    ramp_alphamap2 = grid_alphamap - max(grid_alphamap) .* shift_factor2;
    ramp_alphamap2(ramp_alphamap2 < 0) = 0;

    % normalize
    max_factor2 = 1;
    ramp_alphamap2(ramp_alphamap2 > 0) = ramp_alphamap2(ramp_alphamap2 > 0) ./ max(ramp_alphamap2) .* max_factor2;
    % figure; plot(ramp_alphamap2)

    % volshow(vol2 .^ 1, Parent = viewerThresholded, ...
    %         RenderingStyle = "MaximumIntensityProjection", ...
    %         Colormap = cmap2, ...
    %         Alphamap = "cubic");
    % volshow(vol2 .^ 1, Parent = viewerThresholded, ...
    %         RenderingStyle = "MaximumIntensityProjection", ...
    %         Colormap = cmap2, ...
    %         Alphamap = "quadratic");
    volshow(vol2 .^ 1, Parent = viewerThresholded, ...
            RenderingStyle = "MaximumIntensityProjection", ...
            Colormap = cmap2, ...
            Alphamap = ramp_alphamap2);
            
    % Hard coding for now...
   
    % viewerThresholded.CameraPosition = [91.4553 -69.5688 -11.5372];
    % viewerThresholded.CameraTarget = [40.5000 40.5000 31];
    % viewerThresholded.CameraUpVector = [0.1305 -0.2423 -0.9614];
end