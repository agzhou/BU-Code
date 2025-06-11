dsf = 2; % Downsampling factor
volumeData = single(imresize3(ULMData, 1/dsf, "Method", "cubic"));
% volumeData = downsample(ULMData, dsf);
% volumeData = single(ULMData);
yr = [200, 570];
xr = [50, 824];
zr = [100, 700];
yr = round(yr/dsf);
zr = round(zr/dsf);
xr = round(xr/dsf);

volumeData = volumeData(yr(1):yr(2), xr(1):xr(2), zr(1):zr(2));

figure; imagesc(squeeze(max(volumeData(:, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
figure; imagesc(squeeze(max(volumeData(:, :, :), [], 3))'); colormap(vcmap); clim([0, 40])

%% Smoothing test
addpath('\\ad\eng\users\a\g\agzhou\My Documents\GitHub\BU-Code\Allen code\Processing\neuroelf')
fwhm_size = 11;
smoothVolume = smoothdata3(double(volumeData), fwhm_size);

figure; imagesc(squeeze(max(smoothVolume(:, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
figure; imagesc(squeeze(max(smoothVolume(:, :, :), [], 3))'); colormap(vcmap); clim([0, 40])

%% Do a box filter but emphasize nonzero values
% wss = 51; % window size scalar
% % wss = 5; % window size scalar
% window_size = [wss, wss, wss]; % should all be odd scalars [y, x, z]
% volumeData_box = zeros(size(volumeData), 'single');
% 
% tic
% ds = size(volumeData); % data size
% yds = ds(1);
% xds = ds(2);
% zds = ds(3);
% % for yi = 1:yds
% for yi = 300:310
% %     Yll = max(1, voxelTempCoords(1) - round(window_size(1)/2));
% %     Yul = min(yds, voxelTempCoords(1) + floor(window_size(1)/2));
%     Yll = max(1, yi - round(window_size(1)/2));
%     Yul = min(yds, yi + floor(window_size(1)/2));
% %     volumeData_sliced_temp_y = volumeData(Yll:Yul, :, :);
% 
% %     for xi = 1:xds
%     for xi = 400:510
% %         Xll = max(1, voxelTempCoords(2) - round(window_size(2)/2));                 % X lower limit
% %         Xul = min(xds, voxelTempCoords(2) + floor(window_size(2)/2));  % X upper limit
%         Xll = max(1, xi - round(window_size(2)/2));                 % X lower limit
%         Xul = min(xds, xi + floor(window_size(2)/2));  % X upper limit
% %         volumeData_sliced_temp_x = volumeData_sliced_temp_y(:, Xll:Xul, :);
% 
% %         for zi = 1:zds
% %         tic
%         for zi = 400:510
% %             voxelTempCoords = [yi, xi, zi];
% %             voxelTempInitValue = volumeData(yi, xi, zi);
%                     
% %             Zll = max(1, voxelTempCoords(3) - round(window_size(3)/2));
% %             Zul = min(zds, voxelTempCoords(3) + floor(window_size(3)/2));
%             Zll = max(1, zi - round(window_size(3)/2));
%             Zul = min(zds, zi + floor(window_size(3)/2));
%             
% %             volumeData_sliced_temp_z = volumeData_sliced_temp_x(:, :, Zll:Zul);
% 
% %             windowTemp = volumeData(Yll:Yul, Xll:Xul, Zll:Zul);
% 
% %             windowTemp = volumeData_sliced_temp_z;
% % %             windowTempNonzero = windowTemp(windowTemp > 0);
% %             if isempty(windowTempNonzero) % If there are no nonzero values in the window
% %                 windowTempAvg = 0;
% %             else % If there are nonzero values in the window, store the average
% % %                 windowTempAvg = mean(windowTempNonzero);
% %                 windowTempAvg = mean(windowTemp(windowTemp > 0));
% %             end
% %             volumeData_box(yi, xi, zi) = windowTempAvg;
% %             volumeData_box(yi, xi, zi) = mean(nonzeros(volumeData_sliced_temp_x(:, :, Zll:Zul)));
% 
% %             volumeSlice = volumeData(Yll:Yul, Xll:Xul, Zll:Zul);
% %             volumeSliceNZ = nonzeros(volumeSlice);
%             volumeSliceNZ = nonzeros(volumeData(Yll:Yul, Xll:Xul, Zll:Zul));
% %             volumeData_box(yi, xi, zi) = mean(volumeSliceNZ, 1);
%             volumeData_box(yi, xi, zi) = sum(volumeSliceNZ) / length(volumeSliceNZ);
% %             volumeData_box(yi, xi, zi) = mean(nonzeros(volumeData(Yll:Yul, Xll:Xul, Zll:Zul)));
%         end
% %         toc
%     end
% end
% toc
% % figure; imagesc(squeeze(max(volumeData_box(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
% figure; imagesc(squeeze(max(volumeData_box(:, :, :), [], 1))'); colormap(vcmap); clim([0, 40])

%% Do a box filter but emphasize nonzero values
wss = 25; % window size scalar
% wss = 5; % window size scalar
window_size = [wss, wss, wss]; % should all be odd scalars [y, x, z]
volumeData_box = zeros(size(volumeData), 'single');

tic
ds = size(volumeData); % data size
yds = ds(1);
xds = ds(2);
zds = ds(3);
for yi = 1:yds
% for yi = 300:310
%     Yll = max(1, voxelTempCoords(1) - round(window_size(1)/2));
%     Yul = min(yds, voxelTempCoords(1) + floor(window_size(1)/2));
    Yll = max(1, yi - round(window_size(1)/2));
    Yul = min(yds, yi + floor(window_size(1)/2));
    volumeData_sliced_temp_y = volumeData(Yll:Yul, :, :);

    for xi = 1:xds
%     for xi = 350:510
%         Xll = max(1, voxelTempCoords(2) - round(window_size(2)/2));                 % X lower limit
%         Xul = min(xds, voxelTempCoords(2) + floor(window_size(2)/2));  % X upper limit
        Xll = max(1, xi - round(window_size(2)/2));                 % X lower limit
        Xul = min(xds, xi + floor(window_size(2)/2));  % X upper limit
        volumeData_sliced_temp_x = volumeData_sliced_temp_y(:, Xll:Xul, :);

        for zi = 1:zds
%         tic
%         for zi = 400:510
%             voxelTempCoords = [yi, xi, zi];
%             voxelTempInitValue = volumeData(yi, xi, zi);
                    
%             Zll = max(1, voxelTempCoords(3) - round(window_size(3)/2));
%             Zul = min(zds, voxelTempCoords(3) + floor(window_size(3)/2));
            Zll = max(1, zi - round(window_size(3)/2));
            Zul = min(zds, zi + floor(window_size(3)/2));
            
%             volumeData_sliced_temp_z = volumeData_sliced_temp_x(:, :, Zll:Zul);

%             windowTemp = volumeData(Yll:Yul, Xll:Xul, Zll:Zul);

%             windowTemp = volumeData_sliced_temp_z;
% %             windowTempNonzero = windowTemp(windowTemp > 0);
%             if isempty(windowTempNonzero) % If there are no nonzero values in the window
%                 windowTempAvg = 0;
%             else % If there are nonzero values in the window, store the average
% %                 windowTempAvg = mean(windowTempNonzero);
%                 windowTempAvg = mean(windowTemp(windowTemp > 0));
%             end
%             volumeData_box(yi, xi, zi) = windowTempAvg;
%             volumeData_box(yi, xi, zi) = mean(nonzeros(volumeData_sliced_temp_x(:, :, Zll:Zul)));

%             volumeSlice = volumeData(Yll:Yul, Xll:Xul, Zll:Zul);
%             volumeSliceNZ = nonzeros(volumeSlice);

%             volumeSliceNZ = nonzeros(volumeData(Yll:Yul, Xll:Xul, Zll:Zul));
%             volumeSliceNZ = nonzeros(volumeData_sliced_temp_y(:, Xll:Xul, Zll:Zul));
            volumeSliceNZ = nonzeros(volumeData_sliced_temp_x(:, :, Zll:Zul));
            volumeData_box(yi, xi, zi) = sum(volumeSliceNZ) / length(volumeSliceNZ);
        end
%         toc
    end
end
toc
% figure; imagesc(squeeze(max(volumeData_box(300:500, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
% figure; imagesc(squeeze(max(volumeData_box(:, :, :), [], 1))'); colormap(vcmap); clim([0, 40])
figure; imagesc(squeeze(max(volumeData_box(:, :, :), [], 1))'); colormap(vcmap)
% generateTiffStack_multi({volumeData_box}, [8.8, 8.8, 8], vcmap, 20)

%% Fill test
fillVolume = volumeData_box;
% fillVolume = volumeData;
% fillVolume = fillmissing(fillVolume, 'constant', 100);
fillVolume(fillVolume == 0) = 100;
% volumeViewer(fillVolume)
figure; imagesc(squeeze(max(fillVolume(100:end, :, :), [], 1))'); colormap(vcmap)
figure; imagesc(squeeze(max(fillVolume(:, :, :), [], 3))'); colormap(vcmap)

%% Gradient test
% gradientVolume = volumeData_box;
gradientVolume = imgaussfilt3(volumeData_box, 20);
figure; imagesc(squeeze(max(gradientVolume(:, :, :), [], 1))'); colormap(vcmap)

%%
gradientMethod = 'central';
[Gmag, Gaz, Gelev] = imgradient3(gradientVolume, gradientMethod);
volumeViewer(Gmag)