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

function [coords, XCThresholdsAdaptive] = localizeBubbles3D_framewise_GPU(IQf, refPSF, range, imgRefinementFactor, XCThresholdFactor)
    
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
 
    rfnX = imgRefinementFactor(1); % refinement pixel increase factor
    rfnY = imgRefinementFactor(2);
    rfnZ = imgRefinementFactor(3);

    XCThresholdsAdaptive = zeros(length(framerange), 1);
    coords = cell(length(framerange), 1);

    refPSF_gpu = gpuArray(abs(refPSF));
    tic
%     parfor f = framerange
    parfor f = 1:2
%     for f = 1
        % image refinement/interpolation
        if ~all(imgRefinementFactor == 1)
            I_temp = squeeze(IQs(:, :, :, f));
            refIQs = imresize3(I_temp, [size(I_temp, 1) * rfnX, size(I_temp, 2) * rfnY, size(I_temp, 3) * rfnZ], 'linear');
        else
            refIQs = IQs;
        end
        
        refIQs_gpu = gpuArray(abs(refIQs));

        % and cross correlation
        %   normxcorr3 from file exchange
        %   (https://www.mathworks.com/matlabcentral/fileexchange/73946-normxcorr3-fast-3d-ncc)
        XC = normxcorr3(refPSF_gpu, refIQs_gpu);
%         clear refIQs
    
        % remove data from XC beneath a threshold and find local maxima
        XCt = abs(gather(XC)); % XC Thresholded
        XCThresholdAdaptive = XCThresholdFactor * max(XCt, [], 'all');
        XCThresholdsAdaptive(f) = XCThresholdAdaptive;
    %     XCt(XCt < XCThreshold) = 0;
        XCt(XCt < XCThresholdAdaptive) = 0;
%         clear XC
        centers = imregionalmax(XCt, 6); % Center of each isolated blob (logical matrix)
%         clear XCt

        % Get coordinates from the centers logical matrix
        tsl = size(centers, 1) * size(centers, 2) * size(centers, 3); % troubleshooting length to account for all ones in the centers matrix
        indTemp = find(centers);
        [xc, yc, zc] = ind2sub(size(centers), indTemp);

        if ~((length(xc) == tsl) & (length(yc) == tsl) & (length(zc) == tsl))
            coords{f} = [xc, yc, zc];
        end
%         clear xc yc zc indTemp tsl

    end
    toc
end