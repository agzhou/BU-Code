%% Parameters
zHeightMM = P.endDepthMM - P.startDepthMM;
lateralWidthMM = P.Trans.spacingMm * P.Trans.numelements/2;
hwRatio = zHeightMM / lateralWidthMM;
%%
aPSFn = abs(PSF) ./ max(abs(PSF), [], 'all');
xzFig = figure; imagesc(squeeze(aPSFn(40, :, :))'); colormap gray
title('xz plane')
xlabel('x [mm]')
ylabel('z [mm]')
xzFig.Position(4) = xzFig.Position(3) * hwRatio;

yzFig = figure; imagesc(squeeze(aPSFn(:, 40, :))'); colormap gray
title('yz plane')
xlabel('y [mm]')
ylabel('z [mm]')
yzFig.Position(4) = yzFig.Position(3) * hwRatio;

xyTicksDefault = xticks;
xyTicklabels = cell(length(xyTicksDefault), 1);
for e = 1:numel(xyTicksDefault)
    xyTicklabels{e} = num2str(e * lateralWidthMM / xyTicksDefault(end));
end
xticklabels(xyTicklabels)
yticklabels(xyTicklabels)

%% xy
xyFig = figure; imagesc(squeeze(aPSFn(:, :, 101))); colormap gray
title('xy plane')
xlabel('x [mm]')
ylabel('y [mm]')
% xyFig.Position(4) = xyFig.Position(3);
axis square

xyTicksDefault = xticks;
xyTicklabels = cell(length(xyTicksDefault), 1);
for e = 1:numel(xyTicksDefault)
    xyTicklabels{e} = num2str(e * lateralWidthMM / xyTicksDefault(end));
end
xticklabels(xyTicklabels)
yticklabels(xyTicklabels)