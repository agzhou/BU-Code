% Output: centroid coordinates
% Input:
%         File saving path
%         Filtered IQ data
%         refined PSF (analytical or simulational) - assumed to have the same size
%         as the refined image
%         and pixel spacing as the filtered IQ data
%         range: {z pixel range to use,
%                 x pixel range to use,
%                 frame range to use}
%         image refinement factor [z factor, x factor]
%         Threshold on PSF x refined IQf cross-correlated data to create a
%         binary image
%         Threshold on connected component areas to use

function [centroidCoordinates] = localizeBubbles(IQf, refPSF, range, imgRefinementFactor, binaryThreshold, areaThreshold)
    
    %% Use parallel processing for speed
    % https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox
    
    pp = gcp('nocreate');
    if isempty(pp)
        % There is no parallel pool
        parpool LocalProfile1
    
    end

    %% Section the data to a ROI
    zrange = range{1};
    xrange = range{2};
    framerange = range{3};
    IQs = IQf(zrange, xrange, framerange); % IQ section

    
    %% image refinement/interpolation
    rfnZ = imgRefinementFactor(1); % refinement pixel increase factor
    rfnX = imgRefinementFactor(2);
    
    refIQs = zeros(size(IQs, 1) * rfnZ, size(IQs, 2) * rfnX, size(IQs, 3)); % refined IQ section

    % go through all frames and refine
    parfor f = 1:size(IQs, 3)
        I_temp = IQs(:, :, f);
        refIQs(:, :, f) = imresize(I_temp/max(I_temp, [], 'all') .* 256 .* 5, [size(I_temp, 1) * rfnZ, size(I_temp, 2) * rfnX], 'bilinear');
    end
    
    %% and cross correlation

    XC = normxcorr2(abs(refPSF), abs(refIQs(:, :, 1))); % Cross correlate the filtered/refined images and the simulated PSF
    parfor f = 2:size(refIQs, 3)
        XC(:, :, f) = normxcorr2(abs(refPSF), abs(refIQs(:, :, f)));
    end

    %% binary image conversion of all frames
    
    mask = XC > binaryThreshold;
    bi = zeros(size(XC)); bi(mask) = 1; % binary image with white above the threshold

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

%     zpts = [];
%     xpts = [];
%     for f = 1:nf
%         zpts = [zpts; centroids{f}(:, 1)];
%         xpts = [xpts; centroids{f}(:, 2)];
%     end
%     
%     hPixFactor = 10; % increase the pixel count by this factor in each dimension
%     figure;
%     h = histogram2(zpts, xpts, [size(a, 1) * hPixFactor, size(a, 2) * hPixFactor], 'DisplayStyle','tile');
%     grid off
%     colormap hot

end