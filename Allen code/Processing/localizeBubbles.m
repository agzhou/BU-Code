% Output: centroid coordinates
% Input:
% %         File saving path
%         Filtered IQ data
%         PSF (analytical or simulational) - assumed to have the same size
%         and pixel spacing as the filtered IQ data
%         range: {z pixel range to use,
%                 x pixel range to use,
%                 frame range to use}
%         image refinement factor [z factor, x factor]
%         Threshold on PSF x refined IQf cross-correlated data to create a
%         binary image
%         Threshold on connected component areas to use

function [centroidCoordinates] = localizeBubbles(IQf, PSF, range, imgRefinementFactor, binaryThreshold, areaThreshold)
    
    %% Use parallel processing for speed
    % https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox
    
    pp = gcp('nocreate');
    if isempty(pp)
        % There is no parallel pool
        parpool LocalProfile1
    
    end

    %% Bubble tracking attempt
%     zrange = 40:120;
%     frameRange = 1:size(IQ_f_50_5000, 3);
    zrange = range{1};
    xrange = range{2};
    framerange = range{3};
    IQs = IQf(zrange, xrange, framerange); % IQ section

%     d = diff(IQs, 1, 3); % take first order difference along the frame dimension, seems to get rid of some background
%     d(:, :, end+1) = IQs(:, :, end) - IQs(:, :, end-2); % add another value so you get back to orig # frames
    % figure; imagesc(abs(d(:, :, 1)))
    % figure; imagesc(abs(d(:, :, 2)))
    
    %% image refinement/interpolation
    rfnZ = imgRefinementFactor(1); % refinement pixel increase factor
    rfnX = imgRefinementFactor(2);
    
    refIQs = zeros(size(IQs, 1) * rfnZ, size(IQs, 2) * rfnX, size(IQs, 3)); % refined IQ section
    % I = IQs(:, :, 1);

    % go through all frames and refine
    parfor f = 1:size(IQs, 3)
        I_temp = IQs(:, :, f);
        refIQs(:, :, f) = imresize(I_temp/max(I_temp, [], 'all') .* 256 .* 5, [size(I_temp, 1) * rfnZ, size(I_temp, 2) * rfnX], 'bilinear');
    end
    
    %% and cross correlation
    % if ~exist('base', 'var', 'PSF')
%         load('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\PSF sim\IQ.mat') % Load simulated PSF
    %     PSF = IQ_coherent_sum(90:110, 55:75);                             % Crop the PSF
        refPSF = imresize(PSF, [size(PSF, 1) * rfnZ, size(PSF, 2) * rfnX], 'bilinear');
    %     PSF = refPSF(198:205, 125:132);
    %     PSF = refPSF(160:240, 100:156);

        % refined PSF section
        refPSFs = refPSF(190:210, 118:138);
    % end
    
    XC = normxcorr2(abs(refPSFs), abs(refIQs(:, :, 1))); % Cross correlate the filtered/refined images and the simulated PSF
    parfor f = 2:size(refIQs, 3)
        XC(:, :, f) = normxcorr2(abs(refPSFs), abs(refIQs(:, :, f)));
    end
    
    %% Make video
%     % vo = VideoWriter('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processing\diff_50_5000'); % video object
%     vo = VideoWriter('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processing\XC_50_5000'); % video object
%     
%     % va = abs(d);
%     va = XC - min(XC, [], 'all');
%     va = va ./ max(va, [], 'all');
%     
%     vo.Quality = 100;
%     vo.FrameRate = 100;
%     open(vo);
%     for f = 1:size(va, 3)
%         writeVideo(vo, va(:, :, f));
%     end
%     close(vo);
    
    %% binary image conversion of all frames
    
    mask = XC > binaryThreshold;
    bi = zeros(size(XC)); bi(mask) = 1; % binary image with white above the threshold
    % binaryImage = figure; imagesc(bi); colormap gray

    %% Make video of binary image
%     vo = VideoWriter('G:\Allen\Data\01-17-2025 AZ001 ULM\L22-14v\run 1 allen code left eye\Processing\biXC_50_5000'); % video object
%     va = bi;
%     va = va ./ max(va);
%     
%     vo.Quality = 100;
%     vo.FrameRate = 100;
%     open(vo);
%     for f = 1:size(va, 3)
%         writeVideo(vo, va(:, :, f));
%     end
%     close(vo);
    %% Centroid finding
    nf = size(bi, 3);        % # frames in the binary image stack
    centroids = cell(nf, 1); % initialize
%     areaThreshold = ;
    
    parfor f = 1:nf
        cf = bi(:, :, f); % current frame
        CC = bwconncomp(cf); % connected components (connected regions of 1 in a binary image)
    
        s = regionprops(CC, cf, 'Area', 'WeightedCentroid'); % Get the weighted centroids and area of each connected component
    
        % Remove connected regions without enough pixels (probably noise)
        for si = numel(s):-1:1
            if s(si).Area <= areaThreshold
                s = s(1:si - 1);
            end
        end
    
        centroidsCurrentFrame = zeros(numel(s), 2); % initialize centroid array. Dimensions: # centroids x 2 (z location, x location)
    
        % go through the s structure and get the .WeightedCentroid data
        for cn = 1:numel(s)
            centroidsCurrentFrame(cn, :) = s(cn).WeightedCentroid;
        end
        
        centroids{f} = centroidsCurrentFrame; % put back into overall variable
    
    end
    centroidCoordinates = centroids;
    
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

end