%% Description: look at the raw IQ over time --> check if significant oscillations are present even before clutter filtering

%%

%%
IQ = [];
%%

% for sfi = 1:5
for sfi = 11:15
    IQtemp = single(load([IQpath, IQfilenameStructure, num2str(sfi)], 'IQ').('IQ'));
    IQ = cat(4, IQ, IQtemp);
end

%%
t = (0:size(IQ, 4)-1) ./ P.frameRate;
%% 
refPlane = squeeze(mean(IQ(50, :, :, :), 4));
figure; imagesc(abs(refPlane)'); axis image
%%
% vox = [50, 68, 12 + zstart];
vox = [50, 40, 12 + zstart];
voxData = squeeze(real(IQ(vox(1), vox(2), vox(3), :)));
% voxData_ds = decimate(voxData, P.frameRate); % Decimate/downsample to 1 Hz
figure; plot(t, voxData)

%%
% vox_oob = [50, 60, 14]; % Out of brain
% vox_oob = [50, 70, 14]; % Out of brain
% vox_oob = [50, 10, 7]; % Out of brain
vox_oob = [50, 24, 4]; % Out of brain
figure; plot(t, squeeze(abs(IQ(vox_oob(1), vox_oob(2), vox_oob(3), :))))

%%
% figure; plot(squeeze(mean(abs(IQ(40:50, 60:70, 1:15, :)), [1, 2, 3])))
figure; plot(t, squeeze(mean(abs(IQ(50:60, 30:32, 40:42, :)), [1, 2, 3])))

%%
globalmean = squeeze(mean(IQ, [1, 2, 3]));
figure; plot(t, real(globalmean))
%%
% genSliderV2(squeeze(abs(IQ(50, :, :, :))).^0.3)

generateTiffStack_acrossframes(squeeze(abs(IQ(50, :, :, 1:400))).^0.3, [8.8, 8.8, 8], 'gray', 1:1)