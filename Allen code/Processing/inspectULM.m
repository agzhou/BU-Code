%% Inspect ULM data 
%   - quantify the speeds in __ vessels
%   - Look at the stroke core quantitatively


%% BDM test
BDM = thresholdMaps(SMs_AZ02_day7.SM_SmoothedKF_counter, SMs_AZ02_day7.SM_SmoothedKF_counter, 1, 60);
figure; imagesc(squeeze(max(BDM(300:500, :, :), [], 1))'); colormap hot
%% Manually define the input data for now
% function [] = inspectULM(ULMData)
% ULMData = SMs_reg{4};
% ULMData = SMs_AZ04_day3.SM_SmoothedKF_LI_Rfn;
ULMData = SMs_AZ04_day3.SM_SmoothedKF_LI;
% ULMData = SMs_AZ02_day7.SM_SmoothedKF_LI_Rfn;
% ULMData = SMs_AZ02_day7.SM_SmoothedKF_LI;
% ULMData = thresholdMaps(SMs_AZ02_day7.SM_SmoothedKF_LI, SMs_AZ02_day7.SM_SmoothedKF_counter, 1, 300);
%%
lowspeed_lim = 5; % [mm/s]
% lowspeed_lim = 10; % [mm/s]
ULMData_lowspeed = ULMData;
ULMData_lowspeed(ULMData_lowspeed > lowspeed_lim) = 0;

vcmap = colormap_ULM;
% figure; imagesc(squeeze(max(ULMData(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
% figure; imagesc(squeeze(max(ULMData(:, :, :), [], 3))'); colormap(vcmap); clim([0, 40])
figure; imagesc(squeeze(max(ULMData_lowspeed(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, lowspeed_lim])

% generateTiffStack_multi({ULMData}, [8.8, 8.8, 8], vcmap, 50, [0, 40])
%% Define the stroke core ROI manually
% AZ04 day 3 non-registered
yr = 430:520;
xr = 500:700;
% figure; imagesc(squeeze(max(ULMData(yr, xr, :), [], 1))'); colormap(vcmap); clim([0, 40])
zr = 200:450;
% figure; imagesc(squeeze(max(ULMData(yr, xr, zr), [], 1))'); colormap(vcmap); clim([0, 40])

% AZ02 day 7 non-registered
% yr = 330:380;
% xr = 720:800;  figure; imagesc(squeeze(max(ULMData(yr, xr, :), [], 1))'); colormap(vcmap); clim([0, 40])
% zr = 420:550;

ROI_stroke = ULMData(yr, xr, zr);
figure; imagesc(squeeze(max(ROI_stroke, [], 1))'); colormap(vcmap); clim([0, 40])
ROI_stroke_lowspeed = ULMData_lowspeed(yr, xr, zr);
%% Stroke ROI histograms

% ROI_stroke_lowflow_nz = ROI_stroke(:); % Get (vectorized) values with low nonzero flow 
% ROI_stroke_lowflow_nz = ROI_stroke_lowflow_nz(ROI_stroke_lowflow_nz > 0 & ROI_stroke_lowflow_nz < lowspeed_lim);

figure; histogram(ROI_stroke, 'Normalization', 'pdf'); title('Stroke ROI - all speeds'); xlabel('Flow speed [mm/s]'); ylabel('Probability density function')
figure; histogram(ROI_stroke_lowspeed, 'Normalization', 'pdf'); title("Stroke ROI - speeds up to " + num2str(lowspeed_lim) + " mm/s"); xlabel('Flow speed [mm/s]'); ylabel('Probability density function')
% figure; histogram(ROI_stroke_lowflow_nz)

%% Stroke ROI histograms without 0 flow voxels
figure; histogram(ROI_stroke(ROI_stroke > 0), 'Normalization', 'pdf'); title('Stroke ROI - all speeds'); xlabel('Flow speed [mm/s]'); ylabel('Probability density function')
figure; histogram(ROI_stroke_lowspeed(ROI_stroke_lowspeed > 0), 'Normalization', 'pdf'); title("Stroke ROI - speeds up to " + num2str(lowspeed_lim) + " mm/s"); xlabel('Flow speed [mm/s]'); ylabel('Probability density function')

%% Look at the region corresponding to the stroke core, but flipped across the midline
midline_x = 410; % AZ04 day 3
% midline_x = 490; % AZ02 day 7
ROI_flipped = ULMData(yr, xr - midline_x, zr);
ROI_flipped_lowspeed = ULMData_lowspeed(yr, xr - midline_x, zr);
figure; imagesc(squeeze(max(ROI_flipped, [], 1))'); colormap(vcmap); clim([0, 40])

figure; histogram(ROI_flipped, 'Normalization', 'pdf'); title('Healthy ROI - all speeds'); xlabel('Flow speed [mm/s]'); ylabel('Probability density function')
figure; histogram(ROI_flipped_lowspeed, 'Normalization', 'pdf'); title("Healthy ROI - speeds up to " + num2str(lowspeed_lim) + " mm/s"); xlabel('Flow speed [mm/s]'); ylabel('Probability density function')

%% Healthy ROI histograms without 0 flow voxels
figure; histogram(ROI_flipped(ROI_flipped > 0), 'Normalization', 'pdf'); title('Healthy ROI - all speeds'); xlabel('Flow speed [mm/s]'); ylabel('Probability density function')
figure; histogram(ROI_flipped_lowspeed(ROI_flipped_lowspeed > 0), 'Normalization', 'pdf'); title("Healthy ROI - speeds up to " + num2str(lowspeed_lim) + " mm/s"); xlabel('Flow speed [mm/s]'); ylabel('Probability density function')

%% Plot stroke vs. healthy ROI histograms on top of each other (0 flow excluded)
% All speeds
figure; histogram(ROI_stroke(ROI_stroke > 0), 'Normalization', 'pdf', 'FaceAlpha', 0.7); title('All speeds'); xlabel('Flow speed [mm/s]'); ylabel('Probability density function')
hold on
histogram(ROI_flipped(ROI_flipped > 0), 'Normalization', 'pdf', 'FaceAlpha', 0.7)
legend('Stroke', 'Healthy')

% Low speeds only
figure; histogram(ROI_stroke_lowspeed(ROI_stroke_lowspeed > 0), 'Normalization', 'pdf', 'FaceAlpha', 0.7); title("Speeds up to " + num2str(lowspeed_lim) + " mm/s"); xlabel('Flow speed [mm/s]'); ylabel('Probability density function')
hold on
histogram(ROI_flipped_lowspeed(ROI_flipped_lowspeed > 0), 'Normalization', 'pdf', 'FaceAlpha', 0.7)
legend('Stroke', 'Healthy')


%% Compare persistence

SM_pers_3 = SMs_AZ02_day7_pers_3.SM_SmoothedKF_LI_Rfn;
SM_pers_5 = SMs_AZ02_day7_pers_5.SM_SmoothedKF_LI_Rfn;

BDM_pers_3 = SMs_AZ02_day7_pers_3.SM_SmoothedKF_counter;
BDM_pers_5 = SMs_AZ02_day7_pers_5.SM_SmoothedKF_counter;

SM_test = SM_pers_3 - SM_pers_5;


figure; imagesc(squeeze(max(SM_pers_3(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
figure; imagesc(squeeze(max(SM_test(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])

%% Look at the bubble density map difference
BDM_diff = BDM_pers_3 - BDM_pers_5;
BDM_diff(BDM_diff < 0) = 0;
BDM_diff_yr = 300:450; % y range for the MIP
figure; imagesc(squeeze(max(BDM_pers_3(BDM_diff_yr, :, :), [], 1))' .^ 0.5); colormap hot; % clim([0, 40])
figure; imagesc(squeeze(max(BDM_diff(BDM_diff_yr, :, :), [], 1))' .^ 0.5); colormap hot; % clim([0, 40])
























%% Convolve ULMData with a small blob to fill in missing parts of a vessel
FWHM_X = 30; % x resolution, FWHM-Amplitude, um
FWHM_Y = 30; % x resolution, FWHM-Amplitude, um
FWHM_Z = 30;  % z resolution, FWHM-Amplitude, um
Sigma_X = FWHM_X/(2*sqrt(2*log(2)));
Sigma_Y = FWHM_X/(2*sqrt(2*log(2)));
Sigma_Z = FWHM_Z/(2*sqrt(2*log(2)));
xPSF0 = -30:30; yPSF0 = xPSF0; zPSF0 = xPSF0; % pixels
[xPSF, yPSF, zPSF]=meshgrid(xPSF0, yPSF0, zPSF0);
% PRSSinfo.sysPSF=exp(-((xPSF/(Sigma_X/PRSSinfo.lPix)).^2+(zPSF/(Sigma_Z/PRSSinfo.lPix)).^2)/2);

% hard coding voxel sizes for now..
xvs = 10e-6; % [m]
yvs = 10e-6;
zvs = 10e-6;
xvoxelsPerM = 1/xvs;
yvoxelsPerM = 1/yvs;
zvoxelsPerM = 1/zvs;

smallPSF = exp(-((xPSF/(Sigma_X * (xvoxelsPerM / 1e6))).^2 + (yPSF/(Sigma_Y * (yvoxelsPerM / 1e6))).^2 + (zPSF/(Sigma_Z * (zvoxelsPerM / 1e6))).^2)/2);
figure; imagesc(squeeze(max(smallPSF, [], 1)))

% Convolve with ULMData
ULMData_conv = convn(ULMData, smallPSF, 'same');

% Plot a MIP of the convolved ULMData
figure; imagesc(squeeze(max(ULMData_conv(300:500, :, :), [], 1))'); colormap(vcmap); %clim([0, 40])
%%
flowSpeedLL = 8; % Flow speed lower limit [mm/s]
% flowSpeedUL = 50; % Flow speed upper limit [mm/s]

ULMData_lowspeed = ULMData;
% ULMData_lowspeed = ULMData_conv;
ULMData_lowspeed(ULMData_lowspeed > flowSpeedLL) = 0;

vcmap = colormap_ULM;
figure; imagesc(squeeze(max(ULMData_lowspeed(300:500, :, :), [], 1))'); colormap(vcmap)

%% Segment the ULMData into 2 parts
% % cut_dim = 2; % dimension to define the cut in (1 - y, 2 - x, 3 - z)
midline_cut = 548;
% ULMData_left = ULMData(:, 1:midline_cut, :);
% ULMData_right = ULMData(:, midline_cut:end, :);
% 
% figure; imagesc(squeeze(max(ULMData_left(:, :, :), [], 3))'); colormap(vcmap); clim([0, 40])
% figure; imagesc(squeeze(max(ULMData_right(:, :, :), [], 3))'); colormap(vcmap); clim([0, 40])

%%
figure; imagesc(squeeze(max(ULMData(200:250, :, :), [], 1))'); colormap(vcmap); clim([0, 40])

stroke_rough_ROI_y = 200:250;
stroke_rough_ROI_x = 770:840;
stroke_rough_ROI_z = 420:530;

st_in = ULMData(stroke_rough_ROI_y, stroke_rough_ROI_x, stroke_rough_ROI_z); % stroke test, inside core
st_out = ULMData(stroke_rough_ROI_y, stroke_rough_ROI_x - midline_cut, stroke_rough_ROI_z); % stroke test, outside core

figure; imagesc(squeeze(max(st_in(:, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
figure; imagesc(squeeze(max(st_in(:, :, :), [], 3))'); colormap(vcmap); clim([0, 40])

figure; imagesc(squeeze(max(st_out(:, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
figure; imagesc(squeeze(max(st_out(:, :, :), [], 3))'); colormap(vcmap); clim([0, 40])

%%
nbins = 1000;
bin_limits = [0, 10];
% figure; st_left_hist = histogram(st_left, nbins, 'BinLimits', bin_limits, 'Normalization', 'count');
figure; st_left_hist = histogram(st_in, nbins, 'BinLimits', bin_limits, 'Normalization', 'pdf');
% figure; st_right_hist = histogram(st_right, nbins, 'BinLimits', bin_limits, 'Normalization', 'count');
figure; st_right_hist = histogram(st_out, nbins, 'BinLimits', bin_limits, 'Normalization', 'pdf');










% end

%% Use parallel processing for speed
% https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox

pp = gcp('nocreate');
if isempty(pp)
    % There is no parallel pool
    parpool LocalProfile1

end

%% Blob and gradient for finding the stroke core?
% ULMData = SMs_AZ04_day3.SM_SmoothedKF_LI_Rfn;

% ULMData_blurred = imgaussfilt3(ULMData, 4);
ULMData_test = ULMData;
ULMData_test(ULMData_test < 0) = 0;
ULMData_test(ULMData_test == 0) = min(ULMData_test(ULMData_test > 0), [], 'all');
ULMData_blurred = imboxfilt3(ULMData_test, 11, 'padding', 'symmetric');

% Need a box filter of only the nonzero values within a window...
vcmap = colormap_ULM;
figure; imagesc(squeeze(max(ULMData(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
figure; imagesc(squeeze(max(ULMData_blurred(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])

%% Do a box filter but emphasize nonzero values
wss = 51; % window size scalar
window_size = [wss, wss, wss]; % should all be odd scalars [y, x, z]
ULMData_box = zeros(size(ULMData));

%%%% COULD TRY DOWNSAMPLING FIRST %%%%%

tic
ds = size(ULMData); % data size
yds = ds(1);
xds = ds(2);
zds = ds(3);
for yi = 1:yds
% for yi = 300:410
    Yll = max(1, voxelTempCoords(1) - round(window_size(1)/2));
    Yul = min(yds, voxelTempCoords(1) + floor(window_size(1)/2));
            
    for xi = 1:xds
%     for xi = 400:510
        Xll = max(1, voxelTempCoords(2) - round(window_size(2)/2));                 % X lower limit
        Xul = min(xds, voxelTempCoords(2) + floor(window_size(2)/2));  % X upper limit
        
        for zi = 1:zds
%         tic
%         for zi = 400:510
            voxelTempCoords = [yi, xi, zi];
            voxelTempInitValue = ULMData(yi, xi, zi);
                    
            Zll = max(1, voxelTempCoords(3) - round(window_size(3)/2));
            Zul = min(zds, voxelTempCoords(3) + floor(window_size(3)/2));

            windowTemp = ULMData(Yll:Yul, Xll:Xul, Zll:Zul);
%             windowTempNonzero = windowTemp(windowTemp > 0);
            if isempty(windowTempNonzero) % If there are no nonzero values in the window
                windowTempAvg = 0;
            else % If there are nonzero values in the window, store the average
%                 windowTempAvg = mean(windowTempNonzero);
                windowTempAvg = mean(windowTemp(windowTemp > 0));
            end
            ULMData_box(yi, xi, zi) = windowTempAvg;
        end
%         toc
    end
end
toc
figure; imagesc(squeeze(max(ULMData_box(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
figure; imagesc(squeeze(max(ULMData_box(300, :, :), [], 1))'); colormap(vcmap); clim([0, 40])

%% testing other filters
filter_LP = fspecial3("laplacian", 0.0, 1);
test_LP = imfilter(ULMData, filter_LP);

figure; imagesc(squeeze(max(test_LP(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])

%% Helper functions
function [Tmap] = thresholdMaps(map, counter, lowerCutoff, upperCutoff) % threshold a bubble density map or speed map to remove low and high counts (noise and/or false positives)
    Tmap = map;
    Tmap(counter <= lowerCutoff) = 0;
    Tmap(counter >= upperCutoff) = 0;
end