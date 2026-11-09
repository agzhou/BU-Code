%% Description: writing my own code to call rcabeam from Matlab.
% Dependencies: Parallel Computing Toolbox, Matlab R2026b (for pycupyarray)

%% Add the rcabeam folder (and subfolders) to path
codeDir = cd;
codeDir_split = split(string(codeDir), filesep);
rcabeamCodePath = fullfile(join(codeDir_split(1:find(contains(codeDir_split, "GitHub"))), '\') + "\rcabeam\");
addpath(genpath(rcabeamCodePath))

%% Load RF and params
RFpath = uigetdir('G:\Allen\Data\', 'Select the RF data path');
RFpath = [RFpath, '\'];
if ~exist('P', 'var')
    load([RFpath, 'params.mat']) % load acquisition parameters
end

IQpath = uigetdir(RFpath, 'Select the folder to save reconstructed data to');
IQpath = [IQpath, '\'];

% Prompt for parameter user input
parameterPrompt = {'Start file number', 'End file number'};
parameterDefaults = {'1', ''};
parameterUserInput = inputdlg(parameterPrompt, 'Input Parameters', 1, parameterDefaults);

startFile = str2double(parameterUserInput{1}); % File to start reconstructing from
endFile = str2double(parameterUserInput{2});   % File to stop reconstructing on

clear parameterPrompt parameterDefaults parameterUserInput

%% Get parameters
nf = P.numFramesPerBuffer; % # of frames per superframe
na = P.na; % # of angles (for RC or CR)
nspa = P.Resource.RcvBuffer(1).rowsPerFrame / nf / (2*na); % of samples per acquisition
ntc = P.Resource.Parameters.numRcvChannels; % # of total channels
nc = length(P.Trans.Connector); % # of probe channels (row + col)
wl = P.wl;
Trans = P.Trans;
numElements = P.Trans.numelements/2; % # of rows or cols
angles = P.angles;
c = P.Resource.Parameters.speedOfSound; % Speed of sound [m/s]

% Sampling
sampleMode = P.Receive(1).sampleMode; % Verasonics ADC sample mode
spw = P.Receive(1).samplesPerWave; % Samples per wavelength
decimSampleRate = P.Receive(1).decimSampleRate * 1e6; % [Hz]
demodFreq = P.Receive(1).demodFrequency * 1e6; % [Hz]

%% Load RF file
filenum = 1; % Testing
RFData = loadBinaryRF([RFpath, P.fsFileName, sprintf('%05d', filenum), '.rf'], P);
RFData = unstackFrames(RFData, P);
RFData = reshape(RFData, [nspa, 2*na, ntc, nf]); % Reshape to dimensions (samples, angles, channels, frames)
RFData = RFData(:, :, P.Trans.Connector, :); % Map the receive channels properly

%% If RF is undersampled, do the necessary stuff

sampleDim = 1;
numRFSamples = size(RFData, 1);
sampleVec = (0:numRFSamples-1).';
tSample = sampleVec ./ decimSampleRate; % [Micro]-Time corresponding to each sample

% Demodulate (remove carrier) in the time domain
RF_demod = 2.* single(RFData) .* exp(1i.*2.*pi.* -demodFreq .* tSample); % Multiply by 2 because the original signal is real{envelope * exp(2j*pi*demodFreq*t)} = 1/2*envelope * exp(2j*pi*demodFreq*t) + 1/2*envelope * exp(-2j*pi*demodFreq*t) --> property of complex numbers. Then, demodulating gives you 1/2*envelope + 1/2*envelope*exp(-4j*pi*demodFreq*t), and you keep the first term

% FFT and lowpass the demodulated RF
RFFD = fft(RF_demod, numRFSamples, sampleDim); % RF in the frequency domain
BW = abs(diff(P.Trans.Bandwidth)) * 1e6; % Bandwidth of the probe [Hz]
faxis = linspace(-decimSampleRate/2, decimSampleRate/2, size(RFFD, sampleDim));
LPF_mask = abs(faxis) <= min(BW/2, decimSampleRate/4);
% figure; yyaxis left; plot(faxis, abs(RFFD(:, 1, 120, 1))); yyaxis right; plot(faxis, LPF_mask)
RFFD_LPF = fftshift(RFFD, sampleDim); RFFD_LPF(~LPF_mask, :, :, :) = 0; % Apply the low pass filter
% figure; yyaxis left; plot(faxis, abs(RFFD_LPF(:, 1, 120, 1))); yyaxis right; plot(faxis, LPF_mask)

RF_LPF = ifft(ifftshift(RFFD_LPF, sampleDim), [], sampleDim);

% Get the start sample time of each acquisition and shift the RF for each
% acquisition accordingly (start sample time wrt each acquisition
% corresponds to the first element that fired). Phase depends on the
% carrier frequency, demodFreq.
elementx = (wl .* Trans.ElementPos(1:numElements, 1))';
elementy = (wl .* Trans.ElementPos(numElements + 1 : numElements*2, 2))';
x_end = max(elementx); % End element position [m]
y_end = max(elementy);
d_correction_1way = Trans.lensCorrection * wl; % TX and RX each pass through the lens
pathConst = wl * P.TW.peak + 2 * d_correction_1way - 2 * P.startDepth * wl;
t0 = -(pathConst + x_end * abs(sin(angles))) / c;

IQ = RF_LPF .* repmat(exp(-1i.*2.*pi .* demodFreq .* t0), 1, 2); % Non-beamformed IQ

clear RFData RF_demod RFFD RFFD_LPF RF_LPF
%% Reconstruction region (same pixel grid as RCA_DAS)
% Extents of the region
xsize = numElements * Trans.spacingMm / 1e3;
ysize = xsize;
zsize = (P.endDepth - P.startDepth) * wl;

% Same pixel grid as RCA_DAS (pixel centers of a grid spanning the probe, y increasing)
% xpix_spacing = 0.5 * wl;
% ypix_spacing = 0.5 * wl;
xpix_spacing = Trans.spacingMm/1e3;
ypix_spacing = Trans.spacingMm/1e3;
zpix_spacing = 0.5 * wl;

nz = ceil(zsize / zpix_spacing);
nx = ceil(xsize / xpix_spacing);
ny = ceil(ysize / ypix_spacing);

grid.z = P.startDepth * wl + ((1:nz) - 1/2) / nz * zsize;
grid.x = ((1:nx) - (nx + 1)/2) / nx * xsize;
grid.y = ((1:ny) - (ny + 1)/2) / ny * ysize;

%% Set up the Python environment?
pyenv("Version", "C:\Users\agzhou\AppData\Local\Python\pythoncore-3.13-64\python.exe", ExecutionMode="OutOfProcess");
np = py.importlib.import_module('numpy');
rcb = py.importlib.import_module('rcabeam');

%% A more simple method to use rcabeam
geom = rcb.RCAGeometry(np.array(elementx), np.array(elementy), pyargs( ...
    'c', 1540.0, 'fs', decimSampleRate, 'f_demod', demodFreq, 'fnumber', 1.0));
pygrid = py.tuple({np.array(grid.x), np.array(grid.z), np.array(grid.y)});

rc = permute(IQ(:, 1:na, numElements+1:2*numElements, :), [1, 3, 2, 4]);
cr = permute(IQ(:, na + 1:2*na, 1:numElements, :), [1, 3, 2, 4]);
tStart = t0;
tic
BF = rcb.beamform_ensemble(np.array(rc), np.array(cr), np.array(angles), ...
    tStart, pygrid, geom, pyargs('method', 'opw'));
toc

IQ_BF = permute(double(BF), [1, 3, 2, 4]);

%% A more complicated but faster method to use rcabeam
mod = rcabeam_setup("C:\Users\agzhou\AppData\Local\Python\pythoncore-3.13-64\python.exe");
pyvec = @(v) py.list(num2cell(double(v(:)'))); % Function to convert a Matlab vector into a Python list

f_number = 1.0;
fs = decimSampleRate;

% Create a persistent Session that allocates GPU buffers...
sess = mod.Session(pyvec([numRFSamples, numElements, na, nf]), pyvec(angles), pyvec(t0), ...
    pyvec(grid.x), pyvec(grid.z), pyvec(grid.y), pyvec(elementx), pyvec(elementy), ...
    c, fs, demodFreq, f_number, 'float32');               % 0.788 = Verasonics-equivalent f-number

% Pass the demodulated but nonbeamformed IQ data to the GPU
IQ_RC_gpu = gpuArray(permute(IQ(:, 1:na, numElements+1:2*numElements, :), [4, 2, 3, 1]));
IQ_CR_gpu = gpuArray(permute(IQ(:, na + 1:2*na, 1:numElements, :), [4, 2, 3, 1]));
% IQ_RC_cupy = pycupyarray(IQ_RC_gpu);
% IQ_CR_cupy = pycupyarray(IQ_CR_gpu);

% IQ_RC_gpu, IQ_CR_gpu: gpuArray complex single [nT, na, nc, nsIQ]
tic
vol = rcabeam_from_py(sess.beamform_sum(pycupyarray(IQ_RC_gpu), pycupyarray(IQ_CR_gpu)), [nx ny nz nf]);
toc
% Convert the output back to a Matlab array (gpu or cpu)
IQ_BF_c = gather(vol);                               % conj = Verasonics convention
% toc
sess.release();

