% g1T_fft calculates the temporal field autocorrelation g1, in 1D, 2D, or
% 3D, via the Wiener-Khinchin theorem (FFT of the power spectrum) instead
% of direct per-lag summation. Drop-in replacement for g1T.m: same
% inputs/outputs, same normalization (unbiased, divide by nf - f + 1).
%
% Input: 1. some data, likely IQ (coherently summed across angles).
%           The data should have the spatial dimensions first and then a
%           time/frame dimension last.
%        2. (Optional) number of tau values to calculate g1 at, starting
%           from 0
% Output: temporal g1
%
% Notes:
% - For a given lag m (0-indexed), the circular autocorrelation of data
%   zero-padded to length L equals the desired linear correlation
%   sum_{t=1}^{nf-m} data(t+m)*conj(data(t)) as long as L >= nf + m, so
%   padding to L >= nf + np - 1 covers every lag this function returns.

function [g1] = g1T_fft(data, varargin)
    frameDim = ndims(data);
    nf = size(data, frameDim);

    if nargin > 1 % If the # of points to calculate g1 at is specified
        np = varargin{1};
    else
        np = nf;
    end

    denom = mean(conj(data) .* data, frameDim); % temporal (frame) average

    L = 2^nextpow2(nf + np - 1); % zero-pad enough to avoid circular wraparound for lags 0..np-1
    F = fft(data, L, frameDim);
    C = ifft(F .* conj(F), L, frameDim); % C(...,m+1) = sum_{t=1}^{nf-m} data(t+m)*conj(data(t)), m = 0..L-1

    idx = repmat({':'}, 1, frameDim);
    idx{frameDim} = 1:np;
    numer = C(idx{:});

    normShape = ones(1, frameDim);
    normShape(frameDim) = np;
    normVec = reshape(nf - (0:np - 1), normShape); % (nf - f + 1) for f = 1:np

    g1 = numer ./ normVec ./ denom;
    g1(isnan(g1)) = 0; % Account for any 0/0 issues that result in NaNs
end
