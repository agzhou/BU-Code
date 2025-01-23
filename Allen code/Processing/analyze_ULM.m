
%%

% IQ_coherent_sum = squeeze(sum(IQ, 3));

IQ_coherent_sum = squeeze(sum(IQstack, 3));
% I_coherent_sum = abs(IQ_coherent_sum); % intensity
% figure; imagesc(I_coherent_sum(:, :, 100))

%% SVD proc part 1
tic
[PP, EVs, V_sort] = getSVs1D(IQ_coherent_sum);
disp('SVs decomposed')
toc
%% SVD proc part 2
sv_threshold_lower = 50;
sv_threshold_upper = 2000;
tic
[IQ_f_50_2000] = applySVs1D(IQ_coherent_sum, PP, EVs, V_sort, sv_threshold_lower, sv_threshold_upper);
disp('SVD filtered images put together')
toc
%% Plot filtered data
abs_IQ_f = abs(IQ_f_50_2000);
figure; imagesc(abs_IQ_f(:, :, 1))
figure; imagesc(abs_IQ_f(:, :, 100))
figure; imagesc(abs_IQ_f(:, :, 500))
figure; imagesc(abs_IQ_f(:, :, 1000))
figure; imagesc(abs_IQ_f(:, :, 2000))

%% Plot square of the singular values
figure; plot(EVs, '-o')
title('Eigenvalues')
xlabel('EV #')
ylabel('EV magnitude')

figure; plot(log10(EVs), '-o')
title('Eigenvalues log plot')
xlabel('EV #')
ylabel('log10(EV magnitude)')

%% Make video
vo = VideoWriter('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processing\video_50_5000');
a = abs(IQ_f_50_5000);
a = a ./ max(a); % scale 0-1 so it can make a video

vo.Quality = 100;
vo.FrameRate = 60;
open(vo);
for f = 1:size(a, 3)
    writeVideo(vo, a(:, :, f));
end
close(vo);

%% Bubble tracking attempt
zrange = 40:120;
% frameRange = 1:10000; 
frameRange = 1:size(IQ_f_50_5000, 3);

IQs = IQ_f_50_5000(zrange, :, frameRange); % IQ section
% figure; imagesc(abs(IQs(:, :, 1)))
% figure; imagesc(abs(IQs(:, :, 2)))
% test = img(:, :, 2) - img(:, :, 1);
% figure; imagesc(test)

d = diff(IQs, 1, 3); % take first order difference along the frame dimension, seems to get rid of some background
d(:, :, end+1) = IQs(:, :, end) - IQs(:, :, end-2); % add another value so you get back to orig # frames
% figure; imagesc(abs(d(:, :, 1)))
% figure; imagesc(abs(d(:, :, 2)))

%% Make video
vo = VideoWriter('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processing\diff_50_5000'); % video object
va = abs(d);
va = va ./ max(va);

vo.Quality = 100;
vo.FrameRate = 10;
open(vo);
for f = 1:size(va, 3)
    writeVideo(vo, va(:, :, f));
end
close(vo);

%% binary image conversion of all frames
% a = abs(d);
a = abs(IQs);
% figure; imagesc(a)
% figure; plot(a(:))

threshold = 7e3;

mask = a > threshold;
bi = zeros(size(a)); bi(mask) = 1; % binary image with white above the threshold
% binaryImage = figure; imagesc(bi); colormap gray

%% Make video of binary image
vo = VideoWriter('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processing\bi_50_5000'); % video object
va = bi;
va = va ./ max(va);

vo.Quality = 100;
vo.FrameRate = 10;
open(vo);
for f = 1:size(va, 3)
    writeVideo(vo, va(:, :, f));
end
close(vo);
%% Centroid finding
nf = size(bi, 3);
centroids = cell(nf, 1);
parfor f = 1:nf
    cf = bi(:, :, f); % current frame
    CC = bwconncomp(cf); % connected components (connected regions of 1 in a binary image)
    s = regionprops(CC, cf, 'WeightedCentroid');
    centroidsCurrentFrame = zeros(numel(s), 2); % initialize centroid array. Dimensions: # centroids x 2 (z location, x location)

    % go through the s structure and get the .WeightedCentroid data
    for cn = 1:numel(s)
        centroidsCurrentFrame(cn, :) = s(cn).WeightedCentroid;
    end
    
    centroids{f} = centroidsCurrentFrame; % put back into overall variable

end

%% attempt to graph centroids
% centroidPlotFigure = figure;
% hold on
% for f = 1:nf
%     scatter(centroids{f}(:, 1), centroids{f}(:, 2), '.')
% end
% hold off
% 
zpts = [];
xpts = [];
for f = 1:nf
    zpts = [zpts; centroids{f}(:, 1)];
    xpts = [xpts; centroids{f}(:, 2)];
end

hPixFactor = 10; % increase the pixel count by this factor in each dimension
figure;
h = histogram2(zpts, xpts, [size(a, 1) * hPixFactor, size(a, 2) * hPixFactor], 'DisplayStyle','tile');
grid off
colormap hot

%%


% %% binary image test
% a = abs(d(:, :, 2));
% % figure; imagesc(a)
% % figure; plot(a(:))
% 
% threshold = 1e4;
% 
% mask = a > threshold;
% bi = zeros(size(a)); bi(mask) = 1; % binary image with white above the threshold
% binaryImage = figure; imagesc(bi); colormap gray
% 
% % figure; hold on
% % plot(a(:), 'o')
% % plot(bi(:), '^')
% % hold off
% 
% %%
% CC = bwconncomp(bi)
% s = regionprops(CC, bi, 'WeightedCentroid')
% centroids = zeros(numel(s), 2);
% for cn = 1:numel(s)
%     centroids(cn, :) = s(cn).WeightedCentroid;
% end
% 
% radii = ones(size(centroids, 1), 1);
% hold on
% viscircles(centroids, radii)
% hold off
% % s2 = regionprops(bi,'Area','WeightedCentroid')

%%
taustep = P.SeqControl(1).argument ./ 1e6 .* 1e3 .* P_new.na .* (0:P_new.numSubFrames - 1); % in ms

g1 = g1T_1D(IQ_f);

% pt = [40, 44, 63]; % y, x, z
% pt = [44, 40, 63];
pt = [97, 66];

figure; plot(taustep, abs(squeeze(g1(pt(1), pt(2), :))), '-o')
xlabel('tau (ms)')
ylabel('|g1|')
title(strcat("|g1| at ", num2str(pt), " pixel")) % (z, x)

mag_g1 = abs(g1);

%% Filtered Power Doppler for CBV comparison
IQ_cs_f_sq = IQ_f.^2;
I_f_PowerDoppler = abs(squeeze(sum(IQ_cs_f_sq, 3)));
figure; imagesc(I_f_PowerDoppler)
title('Filtered Power Doppler - xz plane')
xlabel('x pixels')
ylabel('z pixels')
%% Jianbo Power Doppler
[PDI]=sIQ2PDI(IQ_f);

figure; imagesc(PDI(:, :, 1))
title('Power Doppler - positive frequencies')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(PDI(:, :, 2))
title('Power Doppler - negative frequencies')
xlabel('x pixels')
ylabel('z pixels')
figure; imagesc(PDI(:, :, 3))
title('Power Doppler - all frequencies')
xlabel('x pixels')
ylabel('z pixels')

%% same thing but rescaled
dr_scaling = 0.3; % dynamic range scaling factor

% figure; imagesc(PDI(:, :, 1).^dr_scaling)
% title('Power Doppler - positive frequencies')
% xlabel('x pixels')
% ylabel('z pixels')
%
figure; imagesc(PDI(:, :, 2).^dr_scaling)
title('Power Doppler - negative frequencies')
xlabel('x pixels')
ylabel('z pixels')
%
figure; imagesc(PDI(:, :, 3).^dr_scaling)
title('Power Doppler - all frequencies')
xlabel('x pixels')
ylabel('z pixels')
%%
[CBF, CBV] = g1_to_CBi(g1, taustep, 2, 15, 4); % (g1, tau, tau1_index_CBF, tau2_index_CBF, tau1_index_CBV)
%% CBFi
figure; imagesc(CBF)
title('CBFi - xz plane')
xlabel('x pixels')
ylabel('z pixels')

%% CBVi
figure; imagesc(CBV)
title('CBVi - xz plane')
xlabel('x pixels')
ylabel('z pixels')
%%
%%
dr_scaling = 0.3; % dynamic range scaling factor

figure; imagesc(CBV.^dr_scaling)

title(strcat("CBVi - xz plane with exponential scaling factor = ", num2str(dr_scaling)))
xlabel('x pixels')
ylabel('z pixels')
colormap gray