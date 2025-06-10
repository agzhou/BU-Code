%% Inspect ULM data 
%   - quantify the speeds in __ vessels
%   - Look at the stroke core quantitatively



%%
% function [] = inspectULM(ULMData)
ULMData = SMs_reg{4};

vcmap = colormap_ULM;
figure; imagesc(squeeze(max(ULMData(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
figure; imagesc(squeeze(max(ULMData(:, :, :), [], 3))'); colormap(vcmap); clim([0, 40])
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

tic
ds = size(ULMData); % data size
yds = ds(1);
xds = ds(2);
zds = ds(3);
% parfor yi = 1:yds
for yi = 300:410
    Yll = max(1, voxelTempCoords(1) - round(window_size(1)/2));
    Yul = min(yds, voxelTempCoords(1) + floor(window_size(1)/2));
            
%     for xi = 1:xds
    for xi = 400:510
%         for zi = 1:zds
        Xll = max(1, voxelTempCoords(2) - round(window_size(2)/2));                 % X lower limit
        Xul = min(xds, voxelTempCoords(2) + floor(window_size(2)/2));  % X upper limit
        
%         tic
        for zi = 400:510
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
