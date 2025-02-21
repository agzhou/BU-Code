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

function [centers, refIQs, XC] = localizeBubbles2D(IQf, refPSF, range, imgRefinementFactor, XCThreshold)
    
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
        refIQs(:, :, f) = imresize(I_temp, [size(I_temp, 1) * rfnZ, size(I_temp, 2) * rfnX], 'bilinear');
    end
    
    %% and cross correlation

    XC = normxcorr2(abs(refPSF), abs(refIQs(:, :, 1))); % Cross correlate the filtered/refined images and the simulated PSF
    parfor f = 2:size(refIQs, 3)
        XC(:, :, f) = normxcorr2(abs(refPSF), abs(refIQs(:, :, f)));
    end

    %% remove data from XC beneath a threshold and find local maxima
    
    XCt = XC; % XC Thresholded
    XCt(XCt < XCThreshold) = 0;
    centers = imregionalmax(XCt, 4); % Center of each isolated blob (logical matrix)
    
    %% attempt to graph centroids
    
%     test = sum(centers, 3);
%     figure; imagesc(test)
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