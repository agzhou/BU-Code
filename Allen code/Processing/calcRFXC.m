% Calculate RF data cross correlation for a movement metric (see Jonghwan's
% 2011 paper)

% Input: IQ data (currently built for 4D data, but it could easily be
% adapted)

function [rfxc] = calcRFXC(RcvData)
    nf = size(RcvData, length(size(RcvData))); % # of frames (assumed to be the last dimension)
    rfxc = zeros(nf, 1); % Image cross correlation (to the first frame)
    rfref = squeeze(RcvData(:, :, 1)); % reference volume
    rss_rfref = sqrt(sum(abs(rfref).^2, 'all')); % root sum? square of the reference volume
    % tic
    for fi = 1:nf
        rffi = squeeze(RcvData(:, :, fi)); % image #fi
        rss_rffi = sqrt(sum(abs(rffi).^2, 'all')); % root sum? square of volume #fi

        rfxc(fi) = sum( (rffi - mean(rffi, "all")) .* conj(rfref - mean(rfref, "all")), "all") ./ (rss_rfref * rss_rffi);  

    end
    % toc
end