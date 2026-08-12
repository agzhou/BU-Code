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

    noise_floor = median(abs(data_FT{3}), PP.fDim); % Define the noise floor as the median of the FFT'd filtered IQ data across the whole spectrum (after removing system noise)
    
    numer = sum(abs(data_FT{di}) - noise_floor, PP.fDim);
    denom = sum(abs(data_FT{3}) - noise_floor, PP.fDim);
    pnSNR = numer ./ denom; % Frequency-band-power based SNR
    pnSNR_mask = pnSNR > pnSNR_threshold; % Mask for pixels to keep, according to this frequency-based SNR method
    
end