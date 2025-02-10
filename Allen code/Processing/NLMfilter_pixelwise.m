% Pixelwise nonlocal means (NLM) denoising
% Is 2D for now

% function [img_dn] = NLMfilter_pixelwise(img, M, d)
%%
    imgSize = size(img);

    img_dn = zeros(imgSize);
%     img_pad = padarray(img, [d, d], 'symmetric'); % Pad the original image so we can search all local neighborhoods in search regions that include the border of the original image
    img_pad = padarray(img, [d, d]); % Pad the original image so we can search all local neighborhoods in search regions that include the border of the original image
    padImgSize = size(img_pad); % padded size

    searchRegionWidth = 2*M + 1;         % Width of the square or cubic search region
    localNeighborhoodWidth = (2*d + 1);  % Width of the square or cubic local neighborhood region

    %%%%%%%%%%%%%%%%%%%
    % for now, do this but later probably make it vectorized with the frame dimension
    imgv = img(:);
    %%%%%%%%%%%%%%%%%%%

    %%
    for p = 1:imgSize(1) * imgSize(2) % go through every pixel "x_i". p is the vectorized pixel index
%     for p = 1
        pixValue = imgv(p);

        % Coordinates of the current pixel
        [i, j] = ind2sub(imgSize, p); % i is the coordinate of the 1st dimension, j the 2nd dimension
        
        % Define bounds of the current search region
        % CURRENTLY THIS IS 2D ONLY: i -> rows, j -> columns
        imin = max([1, i - M]); % account for the edges of the original image
        imax = min([imgSize(1), i + M]);
        jmin = max([1, j - M]);
        jmax = min([imgSize(2), j + M]);

        imin_pad = max([1, i - M]); % account for the edges] of the original image
        imax_pad = min([padImgSize(1), i + M + d*2]);
        jmin_pad = max([1, j - M]);
        jmax_pad = min([padImgSize(2), j + M + d*2]);

        searchRegion = img(imin:imax, jmin:jmax);
        searchRegionPadded = img_pad(imin_pad:imax_pad, jmin_pad:jmax_pad);

        sRP_center = [ceil(size(searchRegionPadded) ./ 2)]; % indices of the center of the padded search region. Could reformulate this line to be outside the loop for speed
        nbh_p = searchRegionPadded(sRP_center(1) - d : sRP_center(1) + d, sRP_center(2) - d : sRP_center(2) + d); % local neighborhood around p ("x_i")
        h = std(nbh_p(:)) * 12; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Go through each local neighborhood within the search region for
        % pixel p
        pixValue_restored = 0; % "NL" in the paper
        weight_q = NaN(numel(searchRegion), 1);
        for q = 1:numel(searchRegion) % go through each pixel q within the search region
            [m, n] = ind2sub(size(searchRegion), q);
%             pixValue_q = searchRegion(q);
            nbh_q = searchRegionPadded(m : m + localNeighborhoodWidth - 1, n : n + localNeighborhoodWidth - 1);% local neighborhood around q ("x_j")
            weight_q(q) = weightGWED(nbh_p, nbh_q, h);
        end
        
        normalizedWeight_q = weight_q ./ sum(weight_q); % sum of normalized weights = 1
        pixValue_q = searchRegion(:); % vectorized values of pixel qs in the search region for p
        pixValue_restored = sum(pixValue_q .* normalizedWeight_q);

        img_dn(i, j) = pixValue_restored;

    end


% end