function s = synthesizeCorrelatedG1Process(g1fun, frameRate, nFrames)
%% Description:
%   Synthesizes a zero-mean, unit-variance complex Gaussian random process
%   whose autocorrelation matches a prescribed g1(tau) (e.g. from
%   vUS_2D_erf.m), sampled at frameRate for nFrames samples. Standard
%   FFT-based circulant-embedding method: build the two-sided (Hermitian)
%   autocorrelation, take its power spectral density (must be real,
%   non-negative for a valid autocorrelation), and shape white complex
%   Gaussian noise in the frequency domain by its square root.
%
%   Used to simulate realistic sIQ-like time series for testing how the
%   actual g1T.m-style temporal-correlation ESTIMATOR (not just the final
%   g1(tau) curve) responds to different frame-rate / frame-count
%   combinations -- unlike adding noise directly to a clean g1(tau) curve,
%   this correctly captures how MORE frames reduces the estimator's own
%   statistical noise (more frame-pairs averaged per lag) and how a lower
%   frame rate coarsens tau sampling within the decay window.
%
%   Validated (see conversation/commit history) by averaging many
%   noise-free realizations and confirming the recovered g1(tau) matches
%   the target closely at early-to-mid lags (where the decay-curve shape
%   information that v_xgp/v_zgp are fit from actually lives); a small
%   residual mismatch (~1-2% of g1 magnitude) remains at very late lags
%   near the DC floor, attributed to the higher per-realization variance
%   of the near-zero-frequency spectral component representing a
%   persistent (non-decaying) DC term -- immaterial for lags used in
%   fitting, which stay well within the decay window.
%
% Inputs:
%   g1fun: function handle, g1fun(tau) -> complex g1 value(s), vectorized,
%       with g1fun(0) = 1 exactly (handle any 1/tau singularity in the
%       caller, e.g. vUS_2D_erf's tau=0 case)
%   frameRate: [Hz]
%   nFrames: number of time-domain samples to return
%
% Outputs:
%   s: [nFrames, 1] complex, E[|s|^2] ~= 1, autocorrelation ~= g1fun(tau)

    Npad = 2^nextpow2(8 * nFrames); % generous padding so g1 decays well within the buffer
    n = (0:Npad-1)';
    lagIdx = min(n, Npad - n); % circular distance -> two-sided lag
    tauLag = lagIdx / frameRate;

    R = g1fun(tauLag);
    isNeg = n > Npad/2;
    R(isNeg) = conj(R(isNeg)); % Hermitian symmetry: R(-tau) = conj(R(tau))

    S = real(fft(R)); % power spectral density; real up to FFT round-off
    S(S < 0) = 0; % guard tiny negative numerical noise

    whiteNoise = (randn(Npad, 1) + 1i*randn(Npad, 1)) / sqrt(2); % unit-variance complex white noise
    s_full = ifft(sqrt(S) .* whiteNoise) * sqrt(Npad);
    s = s_full(1:nFrames);
end
