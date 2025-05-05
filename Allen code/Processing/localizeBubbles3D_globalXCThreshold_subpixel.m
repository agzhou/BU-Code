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

function [centersRC, refIQs, XC] = localizeBubbles3D_globalXCThreshold(IQf, refPSF, range, imgRefinementFactor, XCThreshold)
    
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
%     IQs = IQf(xrange, yrange, zrange, framerange); % IQ section
    IQs = IQf(yrange, xrange, zrange, framerange); % IQ section

%     IQs = IQf;
    
    %% image refinement/interpolation
    rfnX = imgRefinementFactor(1); % refinement pixel increase factor
    rfnY = imgRefinementFactor(2);
    rfnZ = imgRefinementFactor(3);
    
    if ~all(imgRefinementFactor == 1)
        refIQs = zeros(size(IQs, 1) * rfnX, size(IQs, 2) * rfnY, size(IQs, 3) * rfnZ, size(IQs, 4)); % refined IQ section
    
        % go through all frames and refine
        parfor f = 1:size(IQs, 4)
%         for f = 1
            I_temp = IQs(:, :, :, f);
            refIQs(:, :, :, f) = imresize3(I_temp, [size(I_temp, 1) * rfnX, size(I_temp, 2) * rfnY, size(I_temp, 3) * rfnZ], 'linear');
        end
    else
        refIQs = IQs;
    end
    
    %% and cross correlation

    % normxcorr3 from file exchange
    % (https://www.mathworks.com/matlabcentral/fileexchange/73946-normxcorr3-fast-3d-ncc)
    XC = normxcorr3(abs(refPSF), abs(refIQs(:, :, :, 1)), 'same'); % Cross correlate the filtered/refined images and the simulated PSF
    parfor f = 2:size(refIQs, 4)
        XC(:, :, :, f) = normxcorr3(abs(refPSF), abs(refIQs(:, :, :, f)), 'same');
    end

    %% remove data from XC beneath a threshold and find local maxima
    
    XCt = abs(XC); % XC Thresholded
    XCt(XCt < XCThreshold) = 0;
    centers = imregionalmax(XCt, 6); % Center of each isolated blob (logical matrix)
    
%     figure; imagesc(squeeze(abs(max(XCt(:, :, :, 1), [], 1)))')
%     y_mip_range = 60:100;
%     figure; imagesc(abs(squeeze(max(XC(y_mip_range, :, :, 1), [], 1))' .^ 1))
%     figure; imagesc(abs(squeeze(max(XCt(y_mip_range, :, :, 1), [], 1))' .^ 1))

    %% Try using radial symmetry localization
    % Segment the XC volumes and go through a region around each detected center
    srX = 10; % # x voxels to use in the subregion
    srY = 10; % # x voxels to use in the subregion
    srZ = 10; % # x voxels to use in the subregion
    sr = [srY, srX, srZ];

%     testCenters = centers(:, :, :, 10); %%%%%%%%%%%%%%%%%%%%%%

%     tcc = cell(size(testCenters, 4), 1); % testCenters coords (in each frame)
    tcc = cell(size(centers, 4), 1); % centers coords (in each frame)

    % Get the coords for each bubble
    tsl = size(centers, 1) * size(centers, 2) * size(centers, 3); % troubleshooting length to account for all ones in the centers matrix
    for bfi = 1:size(centers, 4) % buffer frame index
        centersTemp = squeeze(centers(:, :, :, bfi));
        indTemp = find(centersTemp);
        [xc, yc, zc] = ind2sub(size(centersTemp), indTemp);

        if ~((length(xc) == tsl) & (length(yc) == tsl) & (length(zc) == tsl))
            % THIS BELOW LINE MIGHT ONLY WORK IF 'n' STARTS AT
            % 1!!!!!!!!!!!!!!!!!!!!!
%             centerCoords{(n - 1) * P.numFramesPerBuffer + bfi} = [xc, yc, zc];
            tcc{bfi} = [xc, yc, zc];
        end
    end
    clearvars centersTemp indTemp xc yc zc bfi

    %%
    centersRC = cell(size(tcc)); % Initialize the radial centers cell array

    region_limits = [1, size(centers, 2); 1, size(centers, 1); 1, size(centers, 3)];
    for bfi = 1:length(tcc) % Go through each frame
        centersRCbfiTemp = zeros(size(tcc{bfi}));
        for bi = 1:size(tcc{bfi}, 1) % bubble index
%         for bi = 1
            centerTemp = tcc{bfi}(bi, :); % Use the "center" we get from imregionalmax to approximate the local region around each bubble
            
            Xll = max(1, centerTemp(2) - round(srX/2));                 % X lower limit
            Xul = min(size(centers, 2), centerTemp(2) + round(srX/2));  % X upper limit
            Yll = max(1, centerTemp(1) - round(srY/2));
            Yul = min(size(centers, 1), centerTemp(1) + round(srY/2));
            Zll = max(1, centerTemp(3) - round(srZ/2));
            Zul = min(size(centers, 3), centerTemp(3) + round(srZ/2));

            srXL = Xll : Xul; % subregion x limits
            srYL = Yll : Yul; % subregion y limits
            srZL = Zll : Zul; % subregion z limits
    
            srL = {srYL, srXL, srZL};
            % Adjust the sizes so x and y are equal (radialcenter3D.m seems
            % to require this)

            XCsr = XC(srYL, srXL, srZL, bfi);
            pad_sizes = [length(srYL) - srY - 1, length(srXL) - srX - 1, length(srZL) - srZ - 1];
            pad_directions = sign(pad_sizes);
            pad_directions_string = cell(size(pad_directions));
%             for i = 1:length(pad_directions_string)
%                 switch pad_directions(i)
%                     case 1 % 
%                         pad_directions_string{i} = 
                
%             end

            % Set the padarray directions (pre, post, both)
            for dim = 1:3
                if pad_directions(dim) < 0 % If there are fewer XCsr elements in dim than we specified due to being on the edge of the volume
%                     srdim = sr{dim}; % the subregion limits for dimension dim
                    if centerTemp(dim) - round(sr(dim)/2) < region_limits(dim, 1) & centerTemp(dim) + round(sr(dim)/2) > region_limits(dim, 2) % If our subregion goes past the region limits in both directions somehow
                        pad_directions_string{dim} = "both"; % Pad array direction is "both"
                    elseif centerTemp(dim) - round(sr(dim)/2) < region_limits(dim, 1) % If we go lower than the lower limit
                        pad_directions_string{dim} = "pre"; % Pad array direction is "pre"
                    elseif centerTemp(dim) + round(sr(dim)/2) > region_limits(dim, 2) % If we go higher than the upper limit
                        pad_directions_string{dim} = "post"; % Pad array direction is "post"
                    end
                end
            end

            XCsrp = XCsr;
            abs_pad_sizes = diag(abs(pad_sizes));
            for dim = 1:3
                if pad_directions(dim) < 0 % If there are fewer XCsr elements in dim than we specified due to being on the edge of the volume
                    XCsrp = padarray(XCsrp, abs_pad_sizes(:, dim), "replicate", pad_directions_string{dim});
                end
            end
            
%             XCsrp = padarray(XCsr, [length(srYL) - srY), 0, 0], 
            [rc, sigma] = radialcenter3D(XCsrp);
            %%%%%%%%%%%%%%%%%%%% CORRECT FOR HOW IT USES THE CORNER
            %%%%%%%%%%%%%%%%%%%% COORDINATE AS A REFERENCE
            centersRCbfiTemp(bi, :) = rc' + centerTemp; % Store each center (the 'rc' output is relative to the subregion size!!)
        end
        % Remove NaNs
        nanmask = ~isnan(centersRCbfiTemp);
        centersRCbfiTemp_nonan = centersRCbfiTemp(nanmask);
        numNans = sum(nanmask, 1); numNans = size(centersRCbfiTemp, 1) - numNans(1);
        centersRCbfiTemp_nonan = reshape(centersRCbfiTemp_nonan, [size(centersRCbfiTemp, 1) - numNans, size(centersRCbfiTemp, 2)]);
        centersRC{bfi} = centersRCbfiTemp_nonan;
    end
   %% test plot of the RCs
%    fac = 10;
%    newgridsize = size(testCenters) .* fac;
%    newcenters = round(centersRCbfiTemp_nonan .* fac);
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