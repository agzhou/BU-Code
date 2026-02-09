% Calculate image cross correlation for a movement metric (see Jonghwan's
% 2011? paper)

% Input: IQ data (currently built for 4D data, but it could easily be
% adapted)

function [ixc] = calcIXC_simple(IQ)

    nf = size(IQ, length(size(IQ))); % # of frames (assumed to be the last dimension)
    ixc = zeros(nf, 1); % Image cross correlation (to the first frame)

    ref_frame = 1; % Starting reference frame number

    iref = squeeze(IQ(:, :, :, ref_frame)); % reference volume
    rss_iref = sqrt(sum(abs(iref).^2, 'all')); % root sum? square of the reference volume
    % tic
    for fi = 1:nf
    % for fi = 1
        ifi = squeeze(IQ(:, :, :, fi)); % image #fi
        rss_ifi = sqrt(sum(abs(ifi).^2, 'all')); % root sum? square of volume #fi

        % ixc(fi) = sum( (ifi - mean(ifi, "all")) .* conj(iref - mean(iref, "all")), "all") ./ (rss_iref * rss_ifi);  
        ixc(fi) = sum( ifi .* conj(iref), "all") ./ (rss_iref * rss_ifi);  

        % vxc(fi) = normxcorr3(vref, vfi);
    end
    % toc
end