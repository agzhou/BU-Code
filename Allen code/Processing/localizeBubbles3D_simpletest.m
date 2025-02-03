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
    
    %% Use parallel processing for speed
    % https://www.mathworks.com/matlabcentral/answers/91744-how-can-i-check-if-matlabpool-is-running-when-using-parallel-computing-toolbox
    
    pp = gcp('nocreate');
    if isempty(pp)
        % There is no parallel pool
        parpool LocalProfile1
    
    end
    
    %% 3D cross correlation and regionprops testing
    psfXC = normxcorr3(abs(refPSF), abs(refPSF));
    mask = psfXC < 0.4; XCT = psfXC; XCT(mask) = 0;
    bi = zeros(size(psfXC)); bi(~mask) = 1;
    %%
    cf = bi;
    CC = bwconncomp(cf); % connected components (connected regions of 1 in a binary image)
    
    s = regionprops3(CC, cf, 'Volume', 'WeightedCentroid');
    %% and cross correlation

    refIQs = IQf;

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
%     areaThreshold = ;
    
    parfor f = 1:nf
        cf = bi(:, :, :, f); % current frame
        CC = bwconncomp(cf); % connected components (connected regions of 1 in a binary image)
    
        s = regionprops(CC, cf, 'Area', 'WeightedCentroid'); % Get the weighted centroids and area of each connected component
    
        % Remove connected regions without enough pixels (probably noise)
        for si = numel(s):-1:1
            if s(si).Area <= areaThreshold
                s = s(1:si - 1);
            end
        end
    
        centroidsCurrentFrame = zeros(numel(s), 3); % initialize centroid array. Dimensions: # centroids x 2 (x location, y location, z location)
    
        % go through the s structure and get the .WeightedCentroid data
        for cn = 1:numel(s)
            centroidsCurrentFrame(cn, :) = s(cn).WeightedCentroid;
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