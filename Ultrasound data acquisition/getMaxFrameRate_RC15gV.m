% Helper function to get the hardware-limited maximum PRF and frame rate [Hz]
% for the RC15gV (R-C and C-R pairs of acquisitions per angle).
%
% The sequencer cannot start an acquisition sooner than the minimum acquisition
% cycle time T_min after the previous one. 'timeToNextAcq' is only a minimum, so
% if T_min is longer than the programmed time the sequence silently runs slow.
% T_min is modeled as the receive window + a fixed per-acquisition overhead.
% overheadUs was fit to RF time tags of past runs (RC15gV, 47 kHz PRF, 384 samples):
%   13.6 MHz, BS67BW:   window 21.5 us, measured T_min 29.06 us
%   15.625 MHz, BS67BW: window 18.4 us, measured T_min 27.31 us
% 9 us is the conservative (larger) of the two fitted overheads. Not valid for
% the L22-14v, which met a 22.2 us acquisition period with a 16.4 us window.
%
% Inputs:
%   PRF - Pulse repetition frequency [Hz]
%   na - Number of angles
%   numRcvSamples - Samples per acquisition (rounded up to a multiple of 128)
%   probeFreqMHz - Probe center frequency [MHz]
%   samplesPerWave - Nominal samples per wavelength for the ADC sample mode
%   overheadUs - (optional) Per-acquisition overhead [us]

function [lim] = getMaxFrameRate_RC15gV(PRF, na, numRcvSamples, probeFreqMHz, samplesPerWave, overheadUs)

    if nargin < 6
        overheadUs = 9;
    end

    % Decimated sample rate: nearest realizable ADC rate (250 MHz/N) / decimation factor (1-8)
    % to the nominal rate. Reproduces the rates in past runs' Receive(1).decimSampleRate
    % (17.857 MHz at 13.6 MHz and 20.833 MHz at 15.625 MHz, both BS67BW)
    ADCRates = 250 ./ [4, 5, 6, 7, 8]; % [MHz]
    candidateRates = ADCRates(:) ./ (1:8); % [MHz]
    [~, ind] = min(abs(candidateRates(:) - samplesPerWave * probeFreqMHz));
    fsDecim = candidateRates(ind);

    nAcq = 2 * na; % Acquisitions per frame
    T_window = numRcvSamples / fsDecim; % [us]
    T_min = T_window + overheadUs; % [us]

    lim.nAcq = nAcq;
    lim.fsDecim = fsDecim; % [MHz]
    lim.T_window = T_window; % [us]
    lim.T_min = T_min; % [us]
    lim.PRFmax = 1 / (T_min * 1e-6); % [Hz]
    lim.frameRateMax = min(PRF, lim.PRFmax) / nAcq; % [Hz] Keep the gap after the last acquisition >= 1 PRF period
    lim.frameRateMax_atPRFmax = lim.PRFmax / nAcq; % [Hz] If the PRF is lowered to the hardware limit
end
