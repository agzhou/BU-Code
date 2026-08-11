% Description: 
%   Create a mask based on frequency-based SNR for pixels to keep
%   (This is for 2D vUS)
% Input:
%   data_FT: Fourier-transformed data [nz, nx, nt]. According to Jianbo's
%       code, this data should be masked to keep frequencies in some range only.
%   data_FT_ref: Fourier-transformed data [nz, nx, nt], without any masking.
%   PP: Processing Parameters struct with at least the fields:
%       PP.fDim: dimension corresponding to frequency (or time)
%       PP.xDim: dimension corresponding to x in space
%       PP.zDim: dimension corresponding to z in space
%   type: 'full' or 'half' --> dictates some constants used, depending on
%         if we are using the whole frequency spectrum or just part of it (the
%         directional filtering)


function [fbSNR, fbSNR_mask] = spectralSNR(data_FT, data_FT_ref, PP, type)
    zPix = max(floor(PP.zp * 0.1), 1):floor(PP.zp); % z pixels to consider (avoid NaNs and near-field artifacts)

    % Adjust some constants based on the 'type' input
    switch type
        case 'full'
            const = 0.9;
            inflation = 1.02;
        case 'half'
            const = 0.8;
            inflation = 1.05;
    end

    % fbSNR = sum(abs(IQf_FT_separated_masked{3}(zPix, :, :)), fDim) ./ sum(abs(IQf_FT_separated{3}(zPix, :, :)), fDim); % Frequency-based SNR
    fbSNR = sum(abs(data_FT(:, :, :)), PP.fDim) ./ sum(abs(data_FT_ref(:, :, :)), PP.fDim); % Frequency-based SNR
    fbSNR_zAvg = squeeze(mean(fbSNR, PP.xDim)) - const*(1 + ([1:PP.zp]./(5*PP.zp)).^2)' .* std(fbSNR, 0, PP.xDim); % Average the frequency-based SNR across the lateral dimension (x), to get an average value for each depth value (z)
    fbLinearFit = polyfit(zPix, fbSNR_zAvg(zPix), 1); % Linear fit for this z-averaged and std-subtracted threshold vector
    fbSNR_threshold = repmat(polyval(fbLinearFit, 1:PP.zp)', [1, PP.xp]) .* inflation; % Evaluate the linear fit at all z pixels and stretch over x, then add an extra 2%
    fbSNR_mask = fbSNR > fbSNR_threshold; % Mask for pixels to keep, according to this frequency-based SNR method
end