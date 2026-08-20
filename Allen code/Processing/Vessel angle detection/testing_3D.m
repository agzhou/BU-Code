%%
brightondark = true;
tau = 0.5;

% spacing = [PData.PDelta(3) .* P.wl, PData.PDelta(1) .* P.wl]; % z, x pixel size [m]
% if spacing(1) ~= spacing(2)
%     warning('sigma values rely on an equal pixel spacing, which is not the case')
% end
% sigmas = [100:1000].*1e-6; % Sigmas [m]
spacing = [1, 1, 1];
sigmas = [1:0.5:5];

% % Create vectors defining x and z coordinates in actual length units
% x_mm = (1:PData.Size(2)) .* PData.PDelta(1) .* P.wl .* 1e3; % x [mm]
% z_mm = (1:PData.Size(1)) .* PData.PDelta(3) .* P.wl .* 1e3; % z [mm]
%% Plot PDI average (PDIA)
figure; imagesc(squeeze(max(PDIA, [], 1))' .^ 0.5); colormap hot; %axis square; colorbar; xlabel('x [mm]'); ylabel('z [mm]'); title('Power Doppler Average')
%%
vesselness = vesselness2D(PDIA .^ 0.5, sigmas, spacing, tau, brightondark);

figure; imagesc(vesselness)

%%
vesselnessThreshold = 0.01;
minBranchLengthPix = 1;
minSegLengthPix = 3;
[vessels, dir1, dir2, dir3, vesselness, vesselMask, segLabel, skel] = vesselAngle3D(PDIA, sigmas, spacing, tau, brightondark, vesselnessThreshold, minBranchLengthPix, minSegLengthPix);

%%
figure; imagesc(squeeze(max(vesselness, [], 1))'); colormap parula; colorbar; title('Vesselness probability')
figure; imagesc(squeeze(max(vesselMask, [], 1))'); title('Vessel mask')
%%
% figure; h = imagesc(angleMap); colormap hsv; colorbar; axis square; xlabel('x [mm]'); ylabel('z [mm]'); title('Vessel angle'); set(h, 'AlphaData', ~isnan(angleMap)) % make pixels transparent if the angle = NaN
figure; imagesc(squeeze(max(skel, [], 1))'); title('Vessel skeletons')

%%
volumeViewer()
