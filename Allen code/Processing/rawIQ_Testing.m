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



%% See if the SVD clutter filter then introduces fast fluctuations
[xp, yp, zp, nf] = size(IQtemp);
PP = reshape(IQtemp, [xp*yp*zp, nf]);
[U, S, V] = svd(PP, 'econ'); % Already sorted in decreasing order
SVs = diag(S);
% figure; semilogy(SVs, 'LineWidth', 2); xlabel('SV number'); ylabel('SV magnitude')

% -- Some adaptive thresholding stuff -- %
% Plot one SVD subspace as an image
%     subspace = 20;
%     subspace_img = reshape(U(:, subspace) * SVs(subspace) * V(:, subspace)', [xp, yp, zp, nf]);
%     figure; imagesc(squeeze(max(abs(subspace_img(:, :, :, 2)), [], 1))')
% %     volumeViewer(abs(subspace_img(:, :, :, 2)))
% 
SSM = plotSSM(U, false);
% %     SSM = plotSSM(U, true);
%     [~, a_opt, b_opt] = fitSSM(SSM, false); % Get the optimal singular value thresholds
% %     [~, a_opt, b_opt] = fitSSM(SSM, true); % Get the optimal singular value thresholds
%     
sv_threshold_lower = 5; sv_threshold_upper = size(IQtemp, length(size(IQtemp)));
% [IQf, noise] = applySVs2D(IQtemp, PP, SVs, V, sv_threshold_lower, sv_threshold_upper);
%     [IQf, noise] = applySVs2D(IQm, PP, SVs, V, a_opt, b_opt);

[PP, EVs, V_sort] = getSVs2D(IQtemp);
[IQf, noise] = applySVs2D(IQtemp, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
% disp('SVD filtered images put together')

%     volumeViewer(abs(IQf(:, :, :, 1)))
%     figure; imagesc(squeeze(abs(max(IQf(:, :, :, 1), [], 1)))'); colorbar
%     generateTiffStack_acrossframes(abs(IQf), [8.8, 8.8, 8], 'hot', 1:80)
% clearvars IQ

% Use the IQf with separated negative and positive frequency components
%     [IQf_separated, IQf_FT_separated] = separatePosNegFreqs(IQf);

%     [PDI] = calcPowerDoppler(IQf_separated);
PDI = sum(abs(IQf) .^ 2, 4) ./ size(IQf, 4);
%     [CDI] = calcColorDoppler(IQf_FT_separated, P);

%     figure; imagesc(squeeze(max(PDI, [], 1))' .^ 0.5); colormap hot; colorbar
%     figure; imagesc(squeeze(max(PDI ./ noise, [], 1))' .^ 0.5); colormap hot; colorbar
%     volumeViewer(PDI)
%     volumeViewer(PDI ./ noise)

%% Plot the in-brain voxel after SVD clutter filtering, for one superframe
% vox = [50, 68, 12 + zstart];
vox = [50, 40, 12 + zstart];
voxData = squeeze(real(IQf(vox(1), vox(2), vox(3), :)));
% voxData_ds = decimate(voxData, P.frameRate); % Decimate/downsample to 1 Hz
figure; plot(voxData)