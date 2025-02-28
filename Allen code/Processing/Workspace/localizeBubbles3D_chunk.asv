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

% Note: it converts the output coords to uint16, so if for some reason the
% # of x, y, z pixels or frames > 65536, change that. It also converts the
% IQ data into a matrix of singles, so if for some reason the intensity
% goes above 2^32, then change that.

function [coords, img_size, XCThresholdsAdaptive] = localizeBubbles3D_chunk(IQf, refPSF, range, imgRefinementFactor, XCThresholdFactor)
    
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
 
    rfnX = imgRefinementFactor(1); % refinement pixel increase factor
    rfnY = imgRefinementFactor(2);
    rfnZ = imgRefinementFactor(3);

    XCThresholdsAdaptive = zeros(length(framerange), 1);
    coords = cell(length(framerange), 1);                   % Initialize the coordinate cell array (# frames, 1)

    % Define chunk size and each chunk
    chunkSize = 50; % Define chunk size, to help make normxcorr3 more efficient. Each chunk will stack chunkSize frames on top of each other, and then unstack after doing the cross correlation.
%     cs = framerange(1):chunkSize:framerange(end);
    chunks = {};
    numFrames = length(framerange);
    numChunks = ceil(numFrames / chunkSize);

    for chunkNum = 1:numChunks
        chunkStart = (chunkNum - 1) * chunkSize + 1;
        chunkEnd = min(chunkNum * chunkSize, numFrames);
        chunks{chunkNum} = framerange(chunkStart:chunkEnd);
    end

    % Go through each chunk and do the localization
    tic
    for chunkNum = 1:length(chunks)
        chunk = chunks{chunkNum};
        refIQs_chunk = zeros(size(IQs, 1) * rfnX, size(IQs, 2) * rfnY, size(IQs, 3) * rfnZ, length(chunk), 'single');

%         tic
        % image refinement/interpolation
        parfor f = chunk 
%             if ~all(imgRefinementFactor == 1)
                I_temp = squeeze(IQs(:, :, :, f));
                refIQs_chunk(:, :, :, f) = imresize3(I_temp, [size(I_temp, 1) * rfnX, size(I_temp, 2) * rfnY, size(I_temp, 3) * rfnZ], 'linear');
%             else
%                 refIQs(:, :, :, f) = IQs; % comment this out for now
%             end
        end
%         toc
        
        % and cross correlation
        %   normxcorr3 from file exchange
        %   (https://www.mathworks.com/matlabcentral/fileexchange/73946-normxcorr3-fast-3d-ncc)
        rics = size(refIQs_chunk); % refIQs chunk size
        refIQs_chunk_stack = abs(reshape(refIQs_chunk, [rics(1), rics(2), rics(3) * rics(4)])); % Stack the refIQs in the chunk
%         tic
        XC_chunk = normxcorr3(abs(refPSF), refIQs_chunk_stack, 'same'); % Do the cross correlation. The 'same' shape parameter means the output is trimmed to the original data size.
%         toc

%         XC_chunk_trim = XC_chunk(:, :, size(refPSF, 3)/2 - 1:end - size(refPSF, 3)/2 - 1);
%         XCcts = size(XC_chunk_trim); % XC chunk size
%         XC_chunk_rs = reshape(XC_chunk_trim, XCcts(1), XCcts(2), XCcts(3) / length(chunk), length(chunk));

        XCcs = size(XC_chunk); % XC chunk size
        XC_chunk_rs = reshape(XC_chunk, XCcs(1), XCcs(2), XCcs(3) / length(chunk), length(chunk)); % Unstack the XC in the chunk

        % store image size on chunk #1
        if chunkNum == 1
            img_size = size(XC_chunk_rs);
            img_size = img_size(1:3);
        end

        % remove data from XC beneath a threshold and find local maxima
        XC_chunk_t = XC_chunk_rs; % XC chunk, thresholded
        XCThresholdAdaptive = XCThresholdFactor * max(XC_chunk_t, [], 'all'); % Calculate the adaptive XC threshold
        XCThresholdsAdaptive(chunk) = XCThresholdAdaptive;
        XC_chunk_t(XC_chunk_t < XCThresholdAdaptive) = 0;
%         tic
        centers = imregionalmax(XC_chunk_t, 6); % Center of each isolated blob (logical matrix)
%         toc

        % Get coordinates from the centers logical matrix
        tsl = size(centers, 1) * size(centers, 2) * size(centers, 3); % troubleshooting length to account for all ones in the centers matrix
%         tic
        for fi = 1:length(chunk)
            centers_fi = centers(:, :, :, fi);
            f = chunk(fi); % The actual frame number
            indTemp = find(centers_fi); % Get linear indices of the centers' coordinates
            [xc, yc, zc] = ind2sub(size(centers_fi), indTemp); % Turn linear indices into voxel coordinates
    
            if ~((length(xc) == tsl) & (length(yc) == tsl) & (length(zc) == tsl))
                coords{f} = uint16([xc, yc, zc]); % Store the coordinates for frame f
            end
        end
%         toc

%         clear xc yc zc indTemp tsl

    end
    toc
end