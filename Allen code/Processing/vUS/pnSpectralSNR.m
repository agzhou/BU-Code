% Description: 
%   Create a mask based on frequency-band-power based SNR for pixels/voxels to keep
%       --> effectively Eq. 13 in the vUS paper.
%       This is similar to spectralSNR.m, but specifically compares the
%       power of the negative or positive component to that of the whole spectrum.
%   (This is for 2D or 3D vUS)

% Input:
%   data_FT: [3 x 1] cell array of Fourier-transformed data (each cell: 2D: [nz, nx, nt]; 3D: [nx, ny, nz, nt]).
%            The cells correspond to [negative frequencies, positive frequencies, all frequencies] (see separatePosNegFreqs.m)
%            NOTE: this should be after removing system noise!!!!
%   PP: Processing Parameters struct with at least the fields:
%       PP.fDim: dimension corresponding to frequency (or time)
%       PP.xDim: dimension corresponding to x in space
%       PP.yDim: dimension corresponding to y in space (3D only)
%       PP.zDim: dimension corresponding to z in space
%       PP.dimensionality: 2 for 2D, 3 for 3D
%       PP.faxis: frequency vector for the FFT [Hz]
%       PP.freqMask: mask for which parts of the spectrum were set to 0
%   di: direction index --> index for which direction to process.
%       1 --> negative frequencies
%       2 --> positive frequencies

function [pnSNR, pnSNR_mask] = pnSpectralSNR(data_FT, PP, di)

    % Adjust some constants based on the 'type' input
    switch di
        case 1 % Negative frequencies
            pnSNR_threshold = 0.3;
        case 2 % Positive frequencies
            pnSNR_threshold = 0.2;
    end
    
    freqMask_full = ~PP.freqMask; freqMask_full = freqMask_full(:); % Ensure freqMask is a column vector before the reshaping steps below. freqMask is true (1) within the kept frequency bands, which is the opposite of PP.freqMask.
    % switch PP.dimensionality
    %     case 2
    %         freqMask_rep = repmat(permute(freqMask, [2, 3, 1]), [PP.zp, PP.xp, 1]); % Expand the frequency mask to spatial dimensions
    %     case 3
    %         freqMask_rep = repmat(permute(freqMask, [2, 3, 4, 1]), [PP.xp, PP.yp, PP.zp, 1]); % Expand the frequency mask to spatial dimensions
    % end
    data_FT_stacked = stackData(data_FT, PP);
    
    % Stacked in spatial dimensions (vectorized) version e.g., each cell contains [nVox, nFTPts] data
    noise_floor = median(abs(data_FT_stacked{3}(:, freqMask_full)), 2); % Define the noise floor as the median of the FFT'd filtered IQ data across the whole spectrum (after removing system noise)
    
    % Create a separate frequency mask for each directional component, so
    % we don't add a bunch of negative values (with half the spectrum being
    % zeros, minus the positive noise floor)
    freqMask_half = freqMask_full;
    switch di
        case 1 % Negative frequencies
            freqMask_half(PP.faxis > 0) = false;
        case 2 % Positive frequencies
            freqMask_half(PP.faxis < 0) = false;
    end
    
    numer = sum(abs(data_FT_stacked{di}(:, freqMask_half)) - noise_floor, 2);
    denom = sum(abs(data_FT_stacked{3}(:, freqMask_full)) - noise_floor, 2);
    
    % TESTING %
    % numer_test = unstackData(numer, PP);
    % denom_test = unstackData(denom, PP);
    % noise_floor_test = unstackData(noise_floor, PP);
    % figure; imagesc(noise_floor_test)
    % figure; imagesc(numer_test)
    % figure; imagesc(denom_test)
    % tp = [14, 87]; % test pt
    % figure; plot(squeeze(abs(data_FT{3}(tp(1), tp(2), :)))); yline(noise_floor_test(tp(1), tp(2)), 'LineWidth', 2.5, 'Color', 'r')
    % hold on
    % plot(squeeze(abs(data_FT{di}(tp(1), tp(2), :))))
    % hold off

    % figure; plot(squeeze(abs(data_FT{3}(tp(1), tp(2), :))) - noise_floor_test(tp(1), tp(2))); yline(0, 'LineWidth', 2.5, 'Color', 'r')
    % hold on
    % plot(squeeze(abs(data_FT{di}(tp(1), tp(2), :))) - noise_floor_test(tp(1), tp(2)))
    % hold off
    % numer_test(tp(1), tp(2))
    % denom_test(tp(1), tp(2))


    pnSNR = unstackData(numer ./ denom, PP); % Frequency-band-power based SNR

    pnSNR_mask = pnSNR > pnSNR_threshold; % Mask for pixels to keep, according to this frequency-based SNR method
    
    % % Non-stacked version
    % noise_floor = median(abs(data_FT{3}(freqMask_rep)), PP.fDim); % Define the noise floor as the median of the FFT'd filtered IQ data across the whole spectrum (after removing system noise)
    % 
    % numer = sum(abs(data_FT{di}(freqMask_rep)) - noise_floor, PP.fDim);
    % denom = sum(abs(data_FT{3}(freqMask_rep)) - noise_floor, PP.fDim);
    % pnSNR = numer ./ denom; % Frequency-band-power based SNR
    % pnSNR_mask = pnSNR > pnSNR_threshold; % Mask for pixels to keep, according to this frequency-based SNR method
    
end


