% Calculate image cross correlation for a movement metric (see Jonghwan's
% 2011? paper)

% Input: 2x IQ data of the same size (currently built for 4D data, but it could easily be
% adapted)
%    - VOLUME 2 IS THE SHIFTED ONE
%    - Optional input: cropValidVoxels: logical, choose if we do the XC on
%      only valid voxels after shifting or not (the shift code should fill
%      the missing voxels with NaNs)

function [ixc] = calcIXC_shift(vol1, vol2, varargin)
    cropValidVoxels = false;
    if nargin > 2
        cropValidVoxels = varargin{1};
    end

    if ~ cropValidVoxels
        % nf = size(vol1, length(size(vol1))); % # of frames (assumed to be the last dimension)
        % ixc = zeros(nf, 1); % Image cross correlation (to the first frame)
        
        rss_vol1 = sqrt(sum(abs(vol1).^2, 'all')); % root sum? square of volume 1
        rss_vol2 = sqrt(sum(abs(vol2).^2, 'all')); % root sum? square of volume 2
        % tic
        % for fi = 1:nf
        % for fi = 1
            % ifi = squeeze(vol1(:, :, :, fi)); % image #fi
            % rss_ifi = sqrt(sum(abs(ifi).^2, 'all')); % root sum? square of volume #fi
    
            % ixc(fi) = sum( (ifi - mean(ifi, "all")) .* conj(iref - mean(iref, "all")), "all") ./ (rss_vol1 * rss_ifi);  
            ixc = sum( vol1 .* conj(vol2), "all") ./ (rss_vol1 * rss_vol2);  
        % end
        % toc
    else
        validVoxelsMask = ~isnan(vol2);
        vol1_validVoxels = vol1(validVoxelsMask);
        vol2_validVoxels = vol2(validVoxelsMask);

        rss_vol1 = sqrt(sum(abs(vol1_validVoxels).^2, 'all')); % root sum? square of volume 1
        rss_vol2 = sqrt(sum(abs(vol2_validVoxels).^2, 'all')); % root sum? square of volume 2
        ixc = sum( vol1_validVoxels .* conj(vol2_validVoxels), "all") ./ (rss_vol1 * rss_vol2);  
    end
end