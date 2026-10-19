%% Draw ROI
V = PDIA .^ 0.5;
figure; imagesc(squeeze(max(V, [], 1)).'); axis image; title('Coronal MIP (x vs z): click barrel cortex')
[xi, zi] = ginput(1);
figure; imagesc(squeeze(max(V, [], 2)).'); axis image; title('Sagittal MIP (y vs z): click AP position')
[yi, ~] = ginput(1);
c = round([yi, xi, zi]);               % centre voxel (dim 1, 2, 3)

vox_mm = [PData.PDelta(2), PData.PDelta(1), PData.PDelta(3)] * P.wl * 1e3;  % mm per voxel, dims 1-3
semi_mm = [0.5, 0.5, 0.4];             % ROI semi-axes [mm]
[I, J, K] = ndgrid(1:size(V,1), 1:size(V,2), 1:size(V,3));
mask = (((I-c(1))*vox_mm(1)/semi_mm(1)).^2 + ((J-c(2))*vox_mm(2)/semi_mm(2)).^2 + ...
        ((K-c(3))*vox_mm(3)/semi_mm(3)).^2) <= 1;

%%
% Brain mask from mean PDI (tune the threshold; or draw the S1BF ROI directly)
% mask = PDIA{3} > 0.3 * prctile(PDIA{3}(:), 99);
nSF  = size(PDIallSF, 4);

pdiVox = reshape(PDIallSF, [], nSF);   % voxels x superframes
roiTS  = mean(pdiVox(mask(:), :), 1);     % 1 x nSF ROI-mean PDI

[ev, tG, inf] = periStimulusAverage(roiTS, TD.sfStarts, TD.stimOnsetTimestamps, ...
    'SFWidth', sfWidth, 'Pre', 5, 'Post', P.apis.seq_length_s - P.apis.delay_time_ms/1e3, ...
    'Normalize', 'percent');
figure; errorbar(tG, ev, inf.sem); xline(0); xline(P.apis.stim_length_s)
xlabel('Time from stim onset [s]'); ylabel('\DeltaPDI [%]'); title('ROI-mean evoked response')
