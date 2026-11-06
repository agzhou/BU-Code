%% RF Demodulation test

load('C:\Users\agzhou\GitHub\BU-Code\MUST\RFsignal@20MHz.mat')

demodFreq = 5e6; % center frequency (in Hz)
decimSampleRate = 20e6; % sampling frequency (in Hz)

RFData = RF.';
sampleDim = 1;
numRFSamples = size(RFData, 1);
sampleVec = (0:numRFSamples-1).';
tSample = sampleVec ./ decimSampleRate; % [Micro]-Time corresponding to each sample

% Demodulate (remove carrier) in the time domain
RF_demod = 2.* single(RFData) .* exp(1i.*2.*pi.* -demodFreq .* tSample); % Multiply by 2 because the original signal is real{envelope * exp(2j*pi*demodFreq*t)} = 1/2*envelope * exp(2j*pi*demodFreq*t) + 1/2*envelope * exp(-2j*pi*demodFreq*t) --> property of complex numbers. Then, demodulating gives you 1/2*envelope + 1/2*envelope*exp(-4j*pi*demodFreq*t), and you keep the first term

% FFT and lowpass the demodulated RF
RFFD = fft(RF_demod, numRFSamples, sampleDim); % RF in the frequency domain
% BW = abs(diff(P.Trans.Bandwidth)) * 1e6; % Bandwidth of the probe [Hz]
faxis = linspace(-decimSampleRate/2, decimSampleRate/2, size(RFFD, sampleDim));
LPF_mask = abs(faxis) <= decimSampleRate/4;
% figure; yyaxis left; plot(faxis, abs(RFFD(:, 1, 120, 1))); yyaxis right; plot(faxis, LPF_mask)
RFFD_LPF = fftshift(RFFD, sampleDim); RFFD_LPF(~LPF_mask, :, :, :) = 0; % Apply the low pass filter
% figure; yyaxis left; plot(faxis, abs(RFFD_LPF(:, 1, 120, 1))); yyaxis right; plot(faxis, LPF_mask)

RF_LPF = ifft(ifftshift(RFFD_LPF, sampleDim), [], sampleDim);

%% RF Demodulation test: undersampling

load('C:\Users\agzhou\GitHub\BU-Code\MUST\RFsignal@20MHz.mat')

demodFreq = 5e6; % center frequency (in Hz)
decimSampleRate = 20e6; % sampling frequency (in Hz)

RFData = RF.';

% Undersample
factor = 5; % undersampling factor
usRFData = RFData(1:factor:end);
usdecimSampleRate = decimSampleRate/factor; % (under)sampling frequency

sampleDim = 1;
numRFSamples = size(usRFData, 1);
sampleVec = (0:numRFSamples-1).';
tSample = sampleVec ./ usdecimSampleRate; % [Micro]-Time corresponding to each sample

% Demodulate (remove carrier) in the time domain
usRF_demod = 2.* single(usRFData) .* exp(1i.*2.*pi.* -demodFreq .* tSample); % Multiply by 2 because the original signal is real{envelope * exp(2j*pi*demodFreq*t)} = 1/2*envelope * exp(2j*pi*demodFreq*t) + 1/2*envelope * exp(-2j*pi*demodFreq*t) --> property of complex numbers. Then, demodulating gives you 1/2*envelope + 1/2*envelope*exp(-4j*pi*demodFreq*t), and you keep the first term

% FFT and lowpass the demodulated RF
usRFFD = fft(usRF_demod, numRFSamples, sampleDim); % RF in the frequency domain
% BW = abs(diff(P.Trans.Bandwidth)) * 1e6; % Bandwidth of the probe [Hz]
faxis = linspace(-usdecimSampleRate/2, usdecimSampleRate/2, size(usRFFD, sampleDim));
LPF_mask = abs(faxis) <= usdecimSampleRate/4;
% figure; yyaxis left; plot(faxis, abs(RFFD(:, 1, 120, 1))); yyaxis right; plot(faxis, LPF_mask)
usRFFD_LPF = fftshift(usRFFD, sampleDim); usRFFD_LPF(~LPF_mask, :, :, :) = 0; % Apply the low pass filter
% figure; yyaxis left; plot(faxis, abs(RFFD_LPF(:, 1, 120, 1))); yyaxis right; plot(faxis, LPF_mask)

usRF_LPF = ifft(ifftshift(usRFFD_LPF, sampleDim), [], sampleDim);

figure; plot(usRFData)
hold on
plot(abs(usRF_LPF))