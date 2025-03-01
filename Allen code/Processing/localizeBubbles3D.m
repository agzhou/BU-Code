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

function [centers, refIQs, XC, XCThresholdAdaptive] = localizeBubbles3D(IQf, refPSF, range, imgRefinementFactor, XCThresholdFactor)
    
    %% Use parallel processing for speed
    % https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox
    
    pp = gcp('nocreate');
    if isempty(pp)
        % There is no parallel pool
        parpool LocalProfile1
    
    end

    %% Section the data to a ROI
    xrange = range{1};
    yrange = range{2};
    zrange = range{3};
    framerange = range{4};
    IQs = IQf(xrange, yrange, zrange, framerange); % IQ section

%     IQs = IQf;
    
    %% image refinement/interpolation
    rfnX = imgRefinementFactor(1); % refinement pixel increase factor
    rfnY = imgRefinementFactor(2);
    rfnZ = imgRefinementFactor(3);
    
    if ~all(imgRefinementFactor == 1)
        refIQs = zeros(size(IQs, 1) * rfnX, size(IQs, 2) * rfnY, size(IQs, 3) * rfnZ, size(IQs, 4)); % refined IQ section
    
        % go through all frames and refine
        parfor f = 1:size(IQs, 4)
    %     for f = 1
            I_temp = IQs(:, :, :, f);
            refIQs(:, :, :, f) = imresize3(I_temp, [size(I_temp, 1) * rfnX, size(I_temp, 2) * rfnY, size(I_temp, 3) * rfnZ], 'linear');
        end
    else
        refIQs = IQs;
    end
    
    %% and cross correlation

    % normxcorr3 from file exchange
    % (https://www.mathworks.com/matlabcentral/fileexchange/73946-normxcorr3-fast-3d-ncc)
    XC = normxcorr3(abs(refPSF), abs(refIQs(:, :, :, 1))); % Cross correlate the filtered/refined images and the simulated PSF
    parfor f = 2:size(refIQs, 4)
        XC(:, :, :, f) = normxcorr3(abs(refPSF), abs(refIQs(:, :, :, f)));
    end

    %% remove data from XC beneath a threshold and find local maxima
    
    XCt = XC; % XC Thresholded
    XCThresholdAdaptive = XCThresholdFactor * max(XC, [], 'all');
%     XCt(XCt < XCThreshold) = 0;
    XCt(XCt < XCThresholdAdaptive) = 0;
    centers = imregionalmax(XCt, 6); % Center of each isolated blob (logical matrix)
    
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