% Output: centroid coordinates
% Input:
%         File saving path
%         Filtered IQ data (x, y, z pixels)
%         refined PSF (analytical or simulational) - assumed to have the same size
%         as the refined image
%         and pixel spacing as the filtered IQ data
%         range: {x pixel range to use,
%                 y pixel range to use,
%                 z pixel range to use,
%                 frame range to use}
%         image refinement factor [x factor, y factor, z factor]
%         Threshold on PSF x refined IQf cross-correlated data to create a
%         binary image
%         Threshold on connected component areas to use

function [centroidCoordinates] = localizeBubbles3D(IQf, refPSF, range, imgRefinementFactor, binaryThreshold, volumeThreshold)
    
    %% Use parallel processing for speed
    % https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox
    
    pp = gcp('nocreate');
    if isempty(pp)
        % There is no parallel pool
        parpool LocalProfile1
    
    end

    %% Section the data to a ROI
%     zrange = range{1};
%     xrange = range{2};
%     framerange = range{3};
%     IQs = IQf(zrange, xrange, framerange); % IQ section

%     IQs = IQf;
    
    %% image refinement/interpolation
%     rfnX = imgRefinementFactor(1); % refinement pixel increase factor
%     rfnY = imgRefinementFactor(2);
%     rfnZ = imgRefinementFactor(3);
%     
%     refIQs = zeros(size(IQs, 1) * rfnX, size(IQs, 2) * rfnY, size(IQs, 3) * rfnZ, size(IQs, 4)); % refined IQ section
% 
%     % go through all frames and refine
%     parfor f = 1:size(IQs, 4)
% %     for f = 1
%         I_temp = IQs(:, :, :, f);
%         refIQs(:, :, :, f) = imresize3(I_temp/max(I_temp, [], 'all') .* 256 .* 5, [size(IQf, 1) * rfnX, size(IQf, 2) * rfnY, size(IQf, 3) * rfnZ], 'linear');
%     end
    refIQs = IQf;
    
    %% and cross correlation

    % normxcorr3 from file exchange
    % (https://www.mathworks.com/matlabcentral/fileexchange/73946-normxcorr3-fast-3d-ncc)
    XC = normxcorr3(abs(refPSF), abs(refIQs(:, :, :, 1))); % Cross correlate the filtered/refined images and the simulated PSF
    parfor f = 2:size(refIQs, 4)
        XC(:, :, :, f) = normxcorr3(abs(refPSF), abs(refIQs(:, :, :, f)));
    end

    %% binary image conversion of all frames
    
    mask = XC > binaryThreshold;
    bi = zeros(size(XC)); bi(mask) = 1; % binary image with white above the threshold

    %% Centroid finding
    nf = size(bi, 4);        % # frames in the binary image stack
    centroids = cell(nf, 1); % initialize
%     volumeThreshold = ;
    
    parfor f = 1:nf
%     for f = 1:nf
        cf = bi(:, :, :, f); % current frame
        CC = bwconncomp(cf); % connected components (connected regions of 1 in a binary image)
    
        s = regionprops3(CC, cf, 'Volume', 'WeightedCentroid'); % Get the weighted centroids and area of each connected component. It outputs as a table
    
        % Remove connected regions without enough pixels (probably noise)
        for si = size(s, 1):-1:1
            if s(si, :).Volume <= volumeThreshold
                s = s(1:si - 1, :);
            end
        end
    
        centroidsCurrentFrame = zeros(size(s, 1), 3); % initialize centroid array. Dimensions: # centroids x 2 (x location, y location, z location)
    
        % go through the s structure and get the .WeightedCentroid data
        for cn = 1:size(s, 1)
            centroidsCurrentFrame(cn, :) = s(cn, :).WeightedCentroid;
        end
        
        centroids{f} = centroidsCurrentFrame; % put back into overall variable
    
    end
    centroidCoordinates = centroids;
    
    %% attempt to graph centroids

%     xpts = [];
%     ypts = [];
%     zpts = [];
%     for f = 1:nf
%         xpts = [xpts; centroids{f}(:, 1)];
%         ypts = [ypts; centroids{f}(:, 2)];
%         zpts = [zpts; centroids{f}(:, 3)];
%     end
% 
%     figure; scatter3(xpts, ypts, zpts)
%     
%     hPixFactor = 10; % increase the pixel count by this factor in each dimension
%     figure;
%     h = histogram2(zpts, xpts, [size(a, 1) * hPixFactor, size(a, 2) * hPixFactor], 'DisplayStyle','tile');
%     grid off
%     colormap hot

end