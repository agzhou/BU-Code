% load test data
load('F:\Allen\Data\12-16-2024 Phantom\RC15gV\run 2 11 angles -5 to 5 deg\IQ.mat')
load('F:\Allen\Data\12-16-2024 Phantom\RC15gV\run 2 11 angles -5 to 5 deg\params.mat')
%% Initialize variabes and separate R-C and C-R volumes
[xp, yp, zp, nacq, nf] = size(IQ);
na = P.na;

IQ_CR = IQ(:, :, :, 1:na);          % column row volumes
IQ_RC = IQ(:, :, :, na + 1:2*na);   % row column volumes

img_FMAS = zeros(xp, yp, zp, nf);        % initiaize fina FMAS image
%% Do the FMAS
for cri = 1:na       % column row index
    for rci = 1:na   % row column index
        temp = IQ_CR(:, :, :, cri, :) .* IQ_RC(:, :, :, rci, :);
        mag = sqrt(abs(temp));  % magnitude but sqrt to maintain units
        s = sign(temp);         % phase?
        img_FMAS = img_FMAS + s .* mag;   % update the summed variabe

    end
end
%% Plot and compare to normal compounding
abs_DAS = squeeze(abs(sum(IQ, 4)));
abs_FMAS = abs(img_FMAS); % envelope

savepath = 'F:\Allen\Data\12-16-2024 Phantom\RC15gV\run 2 11 angles -5 to 5 deg\FMAS vs DAS\';

%% xz plots
xzDAS = figure; imagesc(squeeze(abs_DAS(40, :, :))' ./ max(abs_DAS(40, :, :), [], 'all'))
title('xz DAS')
ylabel('z pixels')
xlabel('x pixels')
savefig(xzDAS, [savepath, 'xzDAS', '.fig'])
exportgraphics(xzDAS, [savepath, 'xzDAS', '.png'])

xzFMAS = figure; imagesc(squeeze(abs_FMAS(40, :, :))' ./ max(abs_FMAS(40, :, :), [], 'all'))
title('xz FMAS')
ylabel('z pixels')
xlabel('x pixels')
savefig(xzFMAS, [savepath, 'xzFMAS', '.fig'])
exportgraphics(xzFMAS, [savepath, 'xzFMAS', '.png'])

%% yz plots
yzDAS = figure; imagesc(squeeze(abs_DAS(:, 41, :))' ./ max(abs_DAS(:, 41, :), [], 'all'))
title('yz DAS')
ylabel('z pixels')
xlabel('y pixels')
savefig(yzDAS, [savepath, 'yzDAS', '.fig'])
exportgraphics(yzDAS, [savepath, 'yzDAS', '.png'])

yzFMAS = figure; imagesc(squeeze(abs_FMAS(:, 41, :))' ./ max(abs_FMAS(:, 41, :), [], 'all'))
title('yz FMAS')
ylabel('z pixels')
xlabel('y pixels')
savefig(yzFMAS, [savepath, 'yzFMAS', '.fig'])
exportgraphics(yzFMAS, [savepath, 'yzFMAS', '.png'])

%% xy plots
zl = 850;
xyDAS = figure; imagesc(squeeze(abs_DAS(:, :, zl))' ./ max(abs_DAS(:, :, zl), [], 'all'))
title('xy DAS')
axis square
ylabel('x pixels')
xlabel('y pixels')
savefig(xyDAS, [savepath, 'xyDAS', '.fig'])
exportgraphics(xyDAS, [savepath, 'xyDAS', '.png'])

xyFMAS = figure; imagesc(squeeze(abs_FMAS(:, :, zl))' ./ max(abs_FMAS(:, :, zl), [], 'all'))
title('xy FMAS')
axis square
ylabel('x pixels')
xlabel('y pixels')
savefig(xyFMAS, [savepath, 'xyFMAS', '.fig'])
exportgraphics(xyFMAS, [savepath, 'xyFMAS', '.png'])

%% xy maximum projection test
m = max(abs_FMAS, [], 3);
figure; imagesc(m)


%% 1D plots
axial_FMAS = squeeze(abs_FMAS(40, 41, :) ./ max(abs_FMAS(40, 41, :), [], 'all'));
axial_DAS = squeeze(abs_DAS(40, 41, :) ./ max(abs_DAS(40, 41, :), [], 'all'));
axial_pix = 1:size(abs_FMAS, 3);
axialFig = figure; plot(axial_pix, axial_DAS, axial_pix, axial_FMAS, 'LineWidth', 1.5)
title('Axial PSF')
legend('DAS', 'FMAS')
xlabel('z pixel')
ylabel('Normalized Intensity (au)')
savefig(axialFig, [savepath, 'axialPSF', '.fig'])
exportgraphics(axialFig, [savepath, 'axialPSF', '.png'])

lateral_FMAS = squeeze(abs_FMAS(40, :, 849) ./ max(abs_FMAS(40, :, 849), [], 'all'));
lateral_DAS = squeeze(abs_DAS(40, :, 850) ./ max(abs_DAS(40, :, 850), [], 'all'));
lateral_pix = 1:size(abs_FMAS, 2);
lateralFig = figure; plot(lateral_pix, lateral_DAS, lateral_pix, lateral_FMAS, 'LineWidth', 1.5)
title('Lateral PSF')
legend('DAS', 'FMAS')
xlabel('x pixel')
ylabel('Normalized Intensity (au)')
savefig(lateralFig, [savepath, 'lateralPSF', '.fig'])
exportgraphics(lateralFig, [savepath, 'lateralPSF', '.png'])